#!/usr/bin/env python3
"""
化学方程式配平工具 — 纯 Python 标准库实现
算法：fractions.Fraction + 高斯-约当消元 (RREF)，电荷作为特殊元素 e⁻ 纳入矩阵。
"""

import re
import math
import os
import json
from fractions import Fraction
import tkinter as tk
from tkinter import messagebox, filedialog
from tkinter.scrolledtext import ScrolledText

# 可选依赖：Pillow 图片导出
_PIL_AVAILABLE = False
try:
    from PIL import Image, ImageDraw, ImageFont, ImageTk  # noqa: F811
    _PIL_AVAILABLE = True
except ImportError:
    pass


# ============================================================
# 1. 化学式解析
# ============================================================

def parse_atoms(formula: str) -> dict:
    """解析化学式中的原子计数。支持括号如 (NH4)2SO4。"""
    stack: list[dict] = [{}]
    i = 0
    n = len(formula)

    while i < n:
        ch = formula[i]

        if ch.isupper():
            # 新元素开始
            name = ch
            i += 1
            while i < n and formula[i].islower():
                name += formula[i]
                i += 1
            # 读取下标数字
            num_str = ''
            while i < n and formula[i].isdigit():
                num_str += formula[i]
                i += 1
            cnt = int(num_str) if num_str else 1
            stack[-1][name] = stack[-1].get(name, 0) + cnt
            continue  # i 已指向下一个非数字字符

        elif ch == '(':
            stack.append({})

        elif ch == ')':
            group = stack.pop()
            i += 1
            mult_str = ''
            while i < n and formula[i].isdigit():
                mult_str += formula[i]
                i += 1
            mult = int(mult_str) if mult_str else 1
            for elem, cnt in group.items():
                stack[-1][elem] = stack[-1].get(elem, 0) + cnt * mult
            continue  # i 已指向下一个字符

        i += 1

    # 合并栈中所有层（正常情况只剩一层）
    result: dict = {}
    for d in stack:
        for elem, cnt in d.items():
            result[elem] = result.get(elem, 0) + cnt
    return result


def disambiguate_charge(raw_formula: str):
    """
    尝试所有可能的电荷分割，用启发式规则选择最佳解析。
    返回 (base_formula, charge) 或 (raw_formula, 0) 若无电荷。
    """
    candidates: list[tuple[str, int]] = []

    # 从第1个字符开始尝试分割（保证 base 非空）
    for split_pos in range(1, len(raw_formula)):
        suffix = raw_formula[split_pos:]
        # 电荷模式: 可选数字 + 一个或多个 +/-
        m = re.fullmatch(r'(\d*)([+-]+)', suffix)
        if not m:
            continue

        digits = m.group(1)
        signs = m.group(2)

        if digits:
            val = int(digits)
        else:
            val = len(signs)  # "---" → 3

        if '-' in signs:
            val = -val

        base = raw_formula[:split_pos]
        if base:
            candidates.append((base, val))

    if not candidates:
        return raw_formula, 0  # 中性分子

    # 解析每个候选的基本原子组成，过滤不合理者
    parsed_candidates = []
    for base, chg in candidates:
        try:
            atoms = parse_atoms(base)
        except (ValueError, IndexError):
            continue
        max_subscript = max(atoms.values()) if atoms else 0
        total_atoms = sum(atoms.values())
        is_single_element = len(atoms) <= 1
        parsed_candidates.append((base, chg, atoms, total_atoms,
                                  is_single_element, max_subscript))

    # 过滤：单个元素下标不超过 20
    reasonable = [c for c in parsed_candidates if c[5] <= 20]

    if not reasonable:
        # 放宽条件，取 max_subscript 最小的
        parsed_candidates.sort(key=lambda x: x[5])
        reasonable = [parsed_candidates[0]]

    if len(reasonable) == 1:
        return reasonable[0][0], reasonable[0][1]

    # 启发式排序
    # 单元素 → 原子数少的好；多元素 → 原子数多的好
    single = [c for c in reasonable if c[4]]
    multi = [c for c in reasonable if not c[4]]

    if single and multi:
        # 混合情况：优先选单元素中原子最少的，或按整体规则
        # 取所有候选中原子数适中的
        reasonable.sort(key=lambda x: x[3])
        return reasonable[len(reasonable) // 2][0], reasonable[len(reasonable) // 2][1]

    if single:
        single.sort(key=lambda x: x[3])  # 原子少的好
        return single[0][0], single[0][1]

    # multi
    multi.sort(key=lambda x: -x[3])  # 原子多的好
    return multi[0][0], multi[0][1]


def parse_formula(raw: str):
    """
    解析一个化学式，返回 (原子字典, 电荷值)。
    1. 去除状态标记 (s)/(g)/(l)/(aq)
    2. 解析电荷
    3. 解析原子组成
    """
    # 去除状态标记
    cleaned = re.sub(r'\((s|g|l|aq)\)', '', raw)

    # 电荷解析
    base, charge = disambiguate_charge(cleaned)

    # 解析原子组成
    atoms = parse_atoms(base)

    # 基本校验
    if not atoms:
        raise ValueError(f"无法解析化学式: '{raw}'")

    return atoms, charge


# ============================================================
# 2. 方程式解析
# ============================================================

def split_terms(side: str) -> list[str]:
    """
    分割化学式项。只按后面紧跟大写字母或 '(' 的 '+' 分割，
    避免把电荷符号（如 Fe2+ 中的 +）误当作分隔符。
    """
    return [s for s in re.split(r'\+(?=[A-Z(])', side) if s]


def parse_equation(line: str):
    """解析整行方程式，返回 (左边化合物列表, 右边化合物列表)。"""
    # 按 = 或 -> 分割
    parts = re.split(r'->|=', line)
    if len(parts) != 2:
        raise ValueError("方程式格式错误：缺少 '=' 或 '->' 分隔符")

    left_raw, right_raw = parts[0], parts[1]

    # 按 + 分割各项（只分割真正的分隔符 +）
    left_items = split_terms(left_raw)
    right_items = split_terms(right_raw)

    if not left_items or not right_items:
        raise ValueError("方程式格式错误：反应物或生成物为空")

    left_compounds = [parse_formula(item) for item in left_items]
    right_compounds = [parse_formula(item) for item in right_items]

    return left_compounds, right_compounds


# ============================================================
# 3. 矩阵构建
# ============================================================

def build_matrix(left, right):
    """
    构建增广矩阵（非增广，齐次方程组）。
    左边化合物贡献正值，右边贡献负值（移项到左边）。
    """
    # 收集所有唯一元素
    elements = []
    seen = set()
    for atoms, _ch in left + right:
        for elem in atoms:
            if elem not in seen:
                seen.add(elem)
                elements.append(elem)
    # 电荷作为特殊元素
    elements.append('e-')

    n_vars = len(left) + len(right)
    n_elems = len(elements)
    matrix = [[Fraction(0) for _ in range(n_vars)] for _ in range(n_elems)]

    # 左边：正贡献
    for j, (atoms, charge) in enumerate(left):
        for i, elem in enumerate(elements[:-1]):
            matrix[i][j] = Fraction(atoms.get(elem, 0))
        matrix[-1][j] = Fraction(charge)

    # 右边：负贡献（移项）
    offset = len(left)
    for j, (atoms, charge) in enumerate(right):
        col = offset + j
        for i, elem in enumerate(elements[:-1]):
            matrix[i][col] = Fraction(-atoms.get(elem, 0))
        matrix[-1][col] = Fraction(-charge)

    return matrix, elements


# ============================================================
# 4. 高斯-约当消元 (RREF)
# ============================================================

def gauss_jordan_rref(matrix):
    """将矩阵化为行简化阶梯形 (RREF)。"""
    m = len(matrix)       # 行数
    n = len(matrix[0])    # 列数
    r = 0                 # 当前处理行

    for c in range(n):
        # 找主元
        pivot_row = None
        for i in range(r, m):
            if matrix[i][c] != Fraction(0):
                pivot_row = i
                break

        if pivot_row is None:
            continue

        # 交换到当前行
        if pivot_row != r:
            matrix[r], matrix[pivot_row] = matrix[pivot_row], matrix[r]

        # 归一化
        pivot = matrix[r][c]
        for j in range(c, n):
            matrix[r][j] /= pivot

        # 消去其他所有行中该列的值
        for i in range(m):
            if i != r and matrix[i][c] != Fraction(0):
                factor = matrix[i][c]
                for j in range(c, n):
                    matrix[i][j] -= factor * matrix[r][j]

        r += 1
        if r >= m:
            break

    return matrix


# ============================================================
# 5. 求解零空间
# ============================================================

def solve_nullspace(rref_matrix):
    """
    从 RREF 矩阵中提取零空间的一个基向量。
    返回 Fraction 列表，或 None。
    """
    m = len(rref_matrix)
    n = len(rref_matrix[0])

    # 找到主元列和对应的行
    pivot_cols: set = set()
    pivot_map: dict[int, int] = {}  # row → col
    for r in range(m):
        for c in range(n):
            if rref_matrix[r][c] == Fraction(1):
                # 确认是主元（该列其他行都是0）
                is_pivot = True
                for r2 in range(m):
                    if r2 != r and rref_matrix[r2][c] != Fraction(0):
                        is_pivot = False
                        break
                if is_pivot:
                    pivot_cols.add(c)
                    pivot_map[r] = c
                    break

    # 自由变量列
    free_cols = [c for c in range(n) if c not in pivot_cols]

    if len(free_cols) == 0:
        return None  # 只有平凡解

    # 自由度 > 1：设最后一个自由变量为 1，其余为 0
    solution = [Fraction(0)] * n
    solution[free_cols[-1]] = Fraction(1)

    # 回代求解主元变量
    # 找到每个主元行，按从下到上的顺序
    sorted_pivot_rows = sorted(pivot_map.keys(), reverse=True)
    for r in sorted_pivot_rows:
        c = pivot_map[r]
        total = Fraction(0)
        for j in range(c + 1, n):
            total += rref_matrix[r][j] * solution[j]
        solution[c] = -total

    return solution


def to_integer_coeffs(fractions):
    """将 Fraction 列表转换为最简整数列表。"""
    if not fractions:
        return []

    # 求所有分母的 LCM
    denominators = [f.denominator for f in fractions if f != 0]
    if not denominators:
        return [0] * len(fractions)
    lcm = denominators[0]
    for d in denominators[1:]:
        lcm = lcm * d // math.gcd(lcm, d)

    # 乘以 LCM
    ints = [f.numerator * (lcm // f.denominator) if f != 0 else 0
            for f in fractions]

    # 求 GCD 约简
    non_zero = [x for x in ints if x != 0]
    if non_zero:
        g = non_zero[0]
        for x in non_zero[1:]:
            g = math.gcd(g, x)
        ints = [x // g for x in ints]

    return ints


# ============================================================
# 6. 预处理 & 入口
# ============================================================

def preprocess(equation: str) -> str:
    """预处理方程式字符串。"""
    # 删除所有空白字符
    equation = re.sub(r'\s+', '', equation)

    if not equation:
        raise ValueError("输入为空")

    # 检测前置系数
    if re.search(r'(?:^|\+|=)\d+[A-Z(]', equation):
        raise ValueError(
            "请勿输入前置系数（如 2H2O），本工具仅负责配平，"
            "请从化学式本身（如 H2O）开始。"
        )

    return equation


def format_output(left_items, right_items, coeffs):
    """格式化输出方程式。"""
    n_left = len(left_items)
    left_coeffs = coeffs[:n_left]
    right_coeffs = coeffs[n_left:]

    def format_side(items, coeffs):
        parts = []
        for coeff, item in zip(coeffs, items):
            parts.append(f"{coeff}{item}")
        return ' + '.join(parts)

    left_str = format_side(left_items, left_coeffs)
    right_str = format_side(right_items, right_coeffs)
    return f"{left_str} = {right_str}"


def balance(equation: str) -> str:
    """主配平流程。"""
    # 1. 预处理
    eq = preprocess(equation)

    # 2. 解析方程式
    left_compounds, right_compounds = parse_equation(eq)

    # 保存原始化学式字符串用于输出
    parts = re.split(r'->|=', eq)
    left_raw_items = split_terms(parts[0])
    right_raw_items = split_terms(parts[1])

    # 3. 构建矩阵
    matrix, elements = build_matrix(left_compounds, right_compounds)

    # 4. RREF
    rref_matrix = gauss_jordan_rref(matrix)

    # 5. 求解零空间
    solution = solve_nullspace(rref_matrix)

    if solution is None:
        return "错误：该方程式无解或只有平凡解。"

    # 6. 检查解的有效性
    if any(s <= 0 for s in solution):
        return (
            "错误：该方程式配平系数不唯一或存在歧义，"
            "请检查反应物/生成物是否写全。"
        )

    # 7. 转换为最简整数
    int_coeffs = to_integer_coeffs(solution)

    # 8. 格式化输出
    return format_output(left_raw_items, right_raw_items, int_coeffs)


def balance_wrapper(equation_str: str) -> tuple[bool, str]:
    """
    封装 balance()，为 GUI 提供统一的 (成功标志, 消息) 接口。
    """
    try:
        result = balance(equation_str)
    except ValueError as e:
        return False, str(e)
    if result.startswith("错误"):
        return False, result
    return True, result


# ============================================================
# 7. Unicode 上下标辅助（仅 GUI 显示层）
# ============================================================

_SUPER_TABLE = str.maketrans('0123456789+-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻')
_SUB_TABLE   = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')
# 反向映射：Unicode 角标 → 纯文本
_UNICODE_TO_PLAIN = str.maketrans(
    '⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻₀₁₂₃₄₅₆₇₈₉',
    '0123456789+-0123456789'
)


def _digits_to_super(s: str) -> str:
    return s.translate(_SUPER_TABLE)


def _digits_to_sub(s: str) -> str:
    return s.translate(_SUB_TABLE)


def _signs_to_super(s: str) -> str:
    return s.translate(_SUPER_TABLE)


def _term_to_unicode(term: str) -> str:
    """
    将单个化学式项（如 Fe2+、MnO4-、H2O）转换为 Unicode 角标格式。
    前置系数（如 5Fe2+ 中的 5）保持平排，不转上下标。
    利用 disambiguate_charge 正确区分原子下标与电荷数字：
      Fe2+ → Fe²⁺, CO32- → CO₃²⁻, Cr2O72- → Cr₂O₇²⁻
    """
    # 分离前置系数
    m = re.match(r'(\d*)(.*)', term)
    coeff = m.group(1)
    formula = m.group(2)

    # 利用已有的电荷解析器获取正确的 base / charge 分割
    base, charge_val = disambiguate_charge(formula)

    # base 中的数字 → 下标
    body_unicode = _digits_to_sub(base)

    # 中性分子，直接返回
    if charge_val == 0:
        return coeff + body_unicode

    # 电荷 → 上标
    sign = '+' if charge_val > 0 else '-'
    abs_val = abs(charge_val)
    if abs_val > 1:
        charge_unicode = _digits_to_super(str(abs_val)) + _signs_to_super(sign)
    else:
        charge_unicode = _signs_to_super(sign)

    return coeff + body_unicode + charge_unicode


def format_formula_with_unicode(text: str) -> str:
    """将完整结果字符串中的每个化学式项转为 Unicode 角标格式。"""
    parts = re.split(r'( \+ | = )', text)
    result_parts = []
    for part in parts:
        stripped = part.strip()
        if stripped in ('+', '=', ''):
            result_parts.append(part)
        else:
            result_parts.append(_term_to_unicode(stripped))
    return ''.join(result_parts)


# ============================================================
# 8. 自定义分子式持久化
# ============================================================

FORMULAS_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'custom_formulas.json'
)

_DEFAULT_FORMULAS = [
    # 常见分子
    "CO2", "H2", "H2O", "HCl", "HNO3", "H2SO4", "H3PO4",
    "KMnO4", "NaCl", "NaOH", "NH3", "O2", "CH4", "C2H5OH",
    "C6H12O6", "Fe", "Cu", "Al", "Zn", "Ag",
    # 常见离子
    "H+", "Fe2+", "Fe3+", "OH-", "NH4+", "S2-",
    "SO42-", "NO3-", "CO32-", "PO43-", "MnO4-", "Cl-",
    "Cr2O72-", "Cu2+", "Al3+", "Mg2+", "Ca2+", "Na+", "K+",
    # 常见化合物
    "CuSO4", "CaCO3", "Ca(OH)2", "Fe2O3", "Al2O3",
    "SiO2", "Na2CO3", "NaHCO3", "BaSO4", "AgCl",
]


def _load_formulas() -> list[str]:
    try:
        with open(FORMULAS_FILE, 'r', encoding='utf-8') as f:
            data = json.load(f)
            if isinstance(data, list):
                return sorted(set(data))
    except (FileNotFoundError, json.JSONDecodeError):
        pass
    defaults = sorted(set(_DEFAULT_FORMULAS))
    _save_formulas(defaults)
    return defaults


def _save_formulas(formulas: list[str]) -> None:
    with open(FORMULAS_FILE, 'w', encoding='utf-8') as f:
        json.dump(sorted(set(formulas)), f, ensure_ascii=False, indent=2)


def manage_formulas(parent: tk.Tk, entry_widget: ScrolledText) -> None:
    """弹出自定义分子式管理对话框。"""
    formulas = _load_formulas()

    dialog = tk.Toplevel(parent)
    dialog.title("自定义分子式")
    dialog.resizable(False, False)
    dialog.transient(parent)

    dialog.update_idletasks()
    dw, dh = 380, 420
    px = parent.winfo_x() + (parent.winfo_width() - dw) // 2
    py = parent.winfo_y() + (parent.winfo_height() - dh) // 2
    dialog.geometry(f"{dw}x{dh}+{px}+{py}")

    # --- 新增区域 ---
    tk.Label(dialog, text="新分子式：").pack(anchor="w", padx=15, pady=(12, 4))

    add_frame = tk.Frame(dialog)
    add_frame.pack(padx=15, fill="x")

    new_entry = tk.Entry(add_frame, font=("Consolas", 11))
    new_entry.pack(side="left", fill="x", expand=True)

    def do_add() -> None:
        val = new_entry.get().strip()
        if not val:
            messagebox.showwarning("提示", "请输入分子式", parent=dialog)
            return
        if val in formulas:
            messagebox.showinfo("提示", f"'{val}' 已存在", parent=dialog)
            return
        formulas.append(val)
        formulas.sort()
        _save_formulas(formulas)
        _refresh_list()
        new_entry.delete(0, tk.END)
        new_entry.focus_set()

    btn_add = tk.Button(add_frame, text="保存", command=do_add)
    btn_add.pack(side="left", padx=(8, 0))
    new_entry.bind("<Return>", lambda _e: do_add())

    # --- 列表区域（ScrolledText 支持首字母粗体 + 角标）---
    tk.Label(dialog, text='已保存（点击选中，双击或点"插入"添加到输入框）：').pack(
        anchor="w", padx=15, pady=(10, 4))

    formula_text = ScrolledText(dialog, font=("Consolas", 11),
                                 wrap="none", height=14, width=32)
    formula_text.pack(padx=15, fill="both", expand=True)
    formula_text.bind("<Key>", lambda _e: "break")  # 禁止键盘编辑

    formula_text.tag_config("bold", font=("Consolas", 11, "bold"))
    formula_text.tag_config("highlight", background="lightblue")

    _selected_index: list[int | None] = [None]  # mutable 闭包容器

    def _refresh_list() -> None:
        formula_text.config(state="normal")
        formula_text.delete("1.0", tk.END)
        prev_initial: str = ""
        for f in formulas:
            unicode_f = format_formula_with_unicode(f)
            current_initial = unicode_f[0].upper()
            if current_initial != prev_initial:
                # 每组首字母的第一个公式 → 整行粗体
                formula_text.insert(tk.END, unicode_f + "\n", "bold")
                prev_initial = current_initial
            else:
                formula_text.insert(tk.END, unicode_f + "\n")
        formula_text.config(state="disabled")

    def _on_select(_event: tk.Event) -> None:
        formula_text.tag_remove("highlight", "1.0", tk.END)
        index = formula_text.index(f"@{_event.x},{_event.y}")
        line_num = int(index.split(".")[0])
        if 1 <= line_num <= len(formulas):
            _selected_index[0] = line_num - 1
            formula_text.tag_add(
                "highlight",
                f"{line_num}.0",
                f"{line_num}.0 lineend",
            )
        else:
            _selected_index[0] = None

    def do_insert() -> None:
        idx = _selected_index[0]
        if idx is None:
            messagebox.showwarning("提示", "请先在列表中点击选中一个分子式",
                                   parent=dialog)
            return
        formula = formulas[idx]
        col = int(entry_widget.index(tk.INSERT).split(".")[1])
        # 行首或等号后 → 不加 "+"
        text_before = entry_widget.get("insert linestart", tk.INSERT).rstrip()
        if col == 0 or text_before.endswith("="):
            entry_widget.insert(tk.INSERT, formula)
        else:
            entry_widget.insert(tk.INSERT, " + " + formula)
        entry_widget.focus_set()

    def do_delete() -> None:
        idx = _selected_index[0]
        if idx is None:
            messagebox.showwarning("提示", "请先在列表中点击选中一个分子式",
                                   parent=dialog)
            return
        formula = formulas[idx]
        if messagebox.askyesno("确认", f"确定删除 '{formula}'？",
                               parent=dialog):
            formulas.pop(idx)
            _save_formulas(formulas)
            _refresh_list()
            _selected_index[0] = None

    def do_insert_equals() -> None:
        """在光标位置插入 ' = '。"""
        entry_widget.insert(tk.INSERT, " = ")
        entry_widget.focus_set()

    def do_newline() -> None:
        """在光标位置插入换行。"""
        entry_widget.insert(tk.INSERT, "\n")
        entry_widget.focus_set()

    formula_text.bind("<Button-1>", _on_select)
    formula_text.bind("<Double-Button-1>", lambda e: (_on_select(e), do_insert()))

    _refresh_list()

    # --- 底部按钮 ---
    btn_frame = tk.Frame(dialog)
    btn_frame.pack(pady=(8, 12))

    tk.Button(btn_frame, text="插入", width=8, command=do_insert
              ).pack(side="left", padx=3)
    tk.Button(btn_frame, text="插入等号", width=8, command=do_insert_equals
              ).pack(side="left", padx=3)
    tk.Button(btn_frame, text="换行", width=6, command=do_newline
              ).pack(side="left", padx=3)
    tk.Button(btn_frame, text="删除选中", width=8, command=do_delete
              ).pack(side="left", padx=3)
    tk.Button(btn_frame, text="关闭", width=8, command=dialog.destroy
              ).pack(side="left", padx=3)


# ============================================================
# 9. Tkinter 图形界面
# ============================================================

def gui_main():
    """启动 Tkinter 图形界面。"""

    # === 内置工具函数（闭包内，访问 root / entry / output） ===

    def _insert_symbol(symbol: str) -> None:
        """根据焦点将符号插入到输入框或结果框。"""
        widget = root.focus_get()
        if widget not in (entry, output):
            widget = output
        widget.insert(tk.INSERT, symbol)
        widget.focus_set()

    def _handle_reversible() -> None:
        """可逆符号 ⇌：结果框中替换 =/->/→，或光标处插入。"""
        try:
            sel_ranges = output.tag_ranges("sel")
            if sel_ranges:
                start, end = sel_ranges[0], sel_ranges[1]
                selected = output.get(start, end)
                new_text = (selected
                            .replace('\u2192', '\u21cc')   # → → ⇌
                            .replace('->', '\u21cc')
                            .replace('=', '\u21cc'))
                output.delete(start, end)
                output.insert(start, new_text)
            else:
                output.insert(tk.INSERT, '\u21cc')
        except (tk.TclError, ValueError):
            output.insert(tk.INSERT, '\u21cc')
        output.focus_set()

    def _copy_plain_text() -> None:
        """复制结果框纯文本到系统剪贴板。"""
        text = output.get("1.0", tk.END).strip()
        if text:
            root.clipboard_clear()
            root.clipboard_append(text)
            messagebox.showinfo("提示", "已复制到剪贴板", parent=root)
        else:
            messagebox.showwarning("提示", "结果框为空", parent=root)

    # Cross‑platform Times New Roman search
    _FONT_CANDIDATES = [
        r"C:\Windows\Fonts\times.ttf",
        r"C:\Windows\Fonts\timesbd.ttf",
        "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf",
        "/Library/Fonts/Times New Roman.ttf",
    ]

    def _find_times_font() -> str | None:
        for fp in _FONT_CANDIDATES:
            if os.path.isfile(fp):
                return fp
        return None

    # 中文字体（用于条件文字，含 CJK 字符；优先宋体更显学术感）
    _CJK_FONT_CANDIDATES = [
        r"C:\Windows\Fonts\simsun.ttc",     # 宋体（serif，与 Times 同族）
        r"C:\Windows\Fonts\simsunb.ttf",    # SimSun-ExtB
        r"C:\Windows\Fonts\msyh.ttc",       # 微软雅黑
        r"C:\Windows\Fonts\msyhbd.ttc",     # 微软雅黑 Bold
        r"C:\Windows\Fonts\simhei.ttf",     # 黑体
        r"C:\Windows\Fonts\simfang.ttf",    # 仿宋
        r"C:\Windows\Fonts\simkai.ttf",     # 楷体
        "/usr/share/fonts/truetype/wqy/wqy-microhei.ttc",
        "/Library/Fonts/Songti.ttc",
    ]

    def _find_cjk_font() -> str | None:
        for fp in _CJK_FONT_CANDIDATES:
            if os.path.isfile(fp):
                return fp
        return None

    def _render_to_pil_image(lines: list[str], condition: str,
                             *, silent: bool = False) -> Image.Image:
        """将多行方程式渲染为 PIL Image 对象。
        条件文字（如"点燃"）放置在 = / ⇌ / → 的正上方。
        silent=True 时不弹出字体缺失警告。"""
        # 字体准备
        fp = _find_times_font()
        cjk_fp = _find_cjk_font()
        DPI = 300
        factor = DPI / 72.0

        max_chars = max((len(l) for l in lines), default=0)
        base_pt = 12 if max_chars > 70 else 14
        font_px = int(base_pt * factor)
        small_px = int(10 * factor)
        cond_gap_px = int(3 * factor)
        line_spacing_px = int(font_px * 1.5)

        use_default = False
        if fp:
            try:
                font = ImageFont.truetype(fp, font_px)
            except Exception:
                use_default = True
        if not fp or use_default:
            if not silent:
                messagebox.showwarning(
                    "字体警告",
                    "未找到 Times New Roman，使用 PIL 默认字体。",
                    parent=root,
                )
            font = ImageFont.load_default()

        # 条件文字优先用 CJK 字体（支持中文），无则回退到 Times
        if cjk_fp:
            try:
                small_font = ImageFont.truetype(cjk_fp, small_px)
            except Exception:
                small_font = ImageFont.truetype(fp, small_px) if fp \
                    else ImageFont.load_default()
        elif fp:
            small_font = ImageFont.truetype(fp, small_px)
        else:
            small_font = ImageFont.load_default()

        has_cond = bool(condition.strip())
        cond_text = condition.strip() if has_cond else ""

        # 检测条件中是否有 CJK 字符，若无 CJK 字体则警告
        if has_cond and not cjk_fp:
            has_cjk = any('\u4e00' <= ch <= '\u9fff' for ch in cond_text)
            if has_cjk and not silent:
                messagebox.showwarning(
                    "中文字体警告",
                    "未找到中文字体（宋体/微软雅黑等），"
                    "条件中的中文可能显示为方框。",
                    parent=root,
                )
        cond_width = 0.0
        cond_height = 0.0

        draw_tmp = ImageDraw.Draw(Image.new("RGB", (1, 1)))
        if has_cond:
            bb = draw_tmp.textbbox((0, 0), cond_text, font=small_font)
            cond_width = bb[2] - bb[0]
            cond_height = bb[3] - bb[1]

        _SEP_RE = re.compile(r'\s*(⇌|→|=)\s*')
        _SEP_WIDE = {'=': '\u2550\u2550\u2550',
                     '\u2192': '\u27f6',
                     '\u21cc': '\u21cc'}

        class _LineInfo:
            left: str = ''
            sep: str = ''
            right: str = ''
            left_w: float = 0.0
            right_w: float = 0.0
            sep_w: float = 0.0
            total_w: float = 0.0
            eq_height: float = 0.0

        parsed: list[_LineInfo] = []
        for line in lines:
            # 导出图片时省略系数 "1"（例如 1O₂ → O₂），图片更简洁
            line = re.sub(r'\b1(?=[A-Z])', '', line)
            info = _LineInfo()
            m = _SEP_RE.search(line)
            if m:
                info.sep = _SEP_WIDE.get(m.group(1), m.group(1))
                info.left = line[:m.start()]
                info.right = line[m.end():]
                info.left_w = float(
                    draw_tmp.textbbox((0, 0), info.left, font=font)[2])
                sep_str = ' ' + info.sep + ' '
                info.sep_w = float(
                    draw_tmp.textbbox((0, 0), sep_str, font=font)[2])
                info.right_w = float(
                    draw_tmp.textbbox((0, 0), info.right, font=font)[2])
                info.total_w = info.left_w + info.sep_w + info.right_w
                full = info.left + sep_str + info.right
            else:
                info.left = line
                bb = draw_tmp.textbbox((0, 0), line, font=font)
                info.left_w = bb[2] - bb[0]
                info.total_w = info.left_w
                full = line

            bb_eq = draw_tmp.textbbox((0, 0), full, font=font)
            info.eq_height = bb_eq[3] - bb_eq[1]
            parsed.append(info)

        max_line_w = max((p.total_w for p in parsed), default=200.0)
        MARGIN = 40
        img_w = int(max_line_w + 2 * MARGIN)

        total_h = MARGIN
        for i, p in enumerate(parsed):
            if has_cond and p.sep:
                total_h += cond_height + cond_gap_px
            total_h += p.eq_height
            if i < len(parsed) - 1:
                total_h += line_spacing_px
        total_h += MARGIN

        img = Image.new("RGB", (img_w, int(total_h)), (255, 255, 255))
        draw = ImageDraw.Draw(img)

        y = float(MARGIN)
        for i, p in enumerate(parsed):
            if has_cond and p.sep:
                sep_center_x = MARGIN + p.left_w + p.sep_w / 2.0
                cond_x = sep_center_x - cond_width / 2.0
                draw.text((int(cond_x), int(y)), cond_text,
                          font=small_font, fill=(128, 128, 128))
                y += cond_height + cond_gap_px

            if p.sep:
                draw.text((MARGIN, int(y)), p.left, font=font,
                          fill=(0, 0, 0))
                draw.text((MARGIN + int(p.left_w), int(y)),
                          ' ' + p.sep + ' ', font=font, fill=(0, 0, 0))
                draw.text((MARGIN + int(p.left_w + p.sep_w), int(y)),
                          p.right, font=font, fill=(0, 0, 0))
            else:
                draw.text((MARGIN, int(y)), p.left, font=font,
                          fill=(0, 0, 0))

            y += p.eq_height
            if i < len(parsed) - 1:
                y += line_spacing_px

        return img

    def _render_image(lines: list[str], condition: str,
                      filepath: str) -> None:
        """将渲染好的图片保存到文件。"""
        DPI = 300
        img = _render_to_pil_image(lines, condition)
        img.save(filepath, "PNG", dpi=(DPI, DPI))
        messagebox.showinfo("导出成功",
                            f"图片已保存到：\n{filepath}",
                            parent=root)

    def _show_preview(lines: list[str], condition: str) -> bool:
        """预览窗口 — 直接显示已渲染的 PNG 图片（WYSIWYG）。
        返回 True 当用户点击'确认导出'。"""
        confirmed = [False]
        preview = tk.Toplevel(root)
        preview.title("导出预览")
        preview.transient(root)
        preview.resizable(False, False)

        # 渲染原图
        full_img = _render_to_pil_image(lines, condition, silent=True)

        # 缩放至合适预览尺寸
        PREVIEW_W = 720
        if full_img.width > PREVIEW_W:
            ratio = PREVIEW_W / full_img.width
            new_size = (PREVIEW_W, int(full_img.height * ratio))
            preview_img = full_img.resize(new_size, Image.LANCZOS)
        else:
            preview_img = full_img

        photo = ImageTk.PhotoImage(preview_img)
        # 保持引用，防止被 GC
        preview._photo = photo  # type: ignore[attr-defined]

        frame = tk.Frame(preview, bg="white")
        frame.pack(padx=8, pady=(8, 6), fill="both", expand=True)
        tk.Label(frame, image=photo, bg="white").pack()

        btn_frm = tk.Frame(preview)
        btn_frm.pack(pady=(4, 10))

        def _confirm():
            confirmed[0] = True
            preview.destroy()

        tk.Button(btn_frm, text="确认导出", width=12,
                  command=_confirm).pack(side="left", padx=10)
        tk.Button(btn_frm, text="取消", width=12,
                  command=preview.destroy).pack(side="left", padx=10)

        # 居中并按比例设置窗口尺寸
        preview.update_idletasks()
        pad_top = 8 + 6
        btn_h = btn_frm.winfo_reqheight() + 4 + 10
        ph = preview_img.height + pad_top + btn_h + 20
        pw = preview_img.width + 16
        px = root.winfo_x() + (root.winfo_width() - pw) // 2
        py = root.winfo_y() + (root.winfo_height() - ph) // 2
        preview.geometry(f"{pw}x{ph}+{px}+{py}")

        preview.wait_window()
        return confirmed[0]  

    def _do_export() -> None:
        """导出图片：检查依赖 → 读取内容 → 预览 → 保存。"""
        if not _PIL_AVAILABLE:
            messagebox.showerror(
                "缺少依赖",
                "请执行 pip install Pillow 安装图片导出库。",
                parent=root,
            )
            return

        text = output.get("1.0", tk.END).strip()
        if not text:
            messagebox.showwarning("提示", "结果框为空", parent=root)
            return

        # 跳过纯错误消息
        non_error_lines = []
        has_error = False
        for line in text.split('\n'):
            stripped = line.strip()
            if stripped and not any(
                kw in stripped for kw in ('错误', '─')
            ):
                non_error_lines.append(stripped)
            elif stripped.startswith("错误"):
                has_error = True

        if has_error and not non_error_lines:
            messagebox.showwarning("提示",
                                   "结果仅包含错误信息，无法导出。",
                                   parent=root)
            return

        lines = non_error_lines if non_error_lines else \
                [l.strip() for l in text.split('\n') if l.strip()]

        condition = cond_entry.get().strip()

        if not _show_preview(lines, condition):
            return

        filepath = filedialog.asksaveasfilename(
            parent=root,
            defaultextension=".png",
            filetypes=[("PNG 图片", "*.png")],
            initialfile="balanced_equation.png",
        )
        if not filepath:
            return

        _render_image(lines, condition, filepath)

    # === 核心逻辑 ===

    def do_balance():
        """执行配平并显示结果。支持多行批量输入。"""
        raw = entry.get("1.0", tk.END).strip()
        if not raw:
            messagebox.showwarning("提示", "请输入方程式")
            return

        # 输出区域清空旧内容并重置撤销栈
        output.config(state="normal")
        output.delete("1.0", tk.END)
        output.edit_reset()

        lines = raw.split('\n')
        unicode_lines = []
        needs_sep = False

        for line in lines:
            line = line.strip()
            if not line:
                continue

            plain = line.translate(_UNICODE_TO_PLAIN)

            if needs_sep:
                output.insert(tk.END, "\n" + "─" * 55 + "\n")
            needs_sep = True

            success, msg = balance_wrapper(plain)

            if success:
                output.insert(tk.END, format_formula_with_unicode(msg))
            else:
                output.insert(tk.END, msg, "error")

            unicode_lines.append(format_formula_with_unicode(line))

        # 输入框同步为 Unicode 角标版
        entry.delete("1.0", tk.END)
        if unicode_lines:
            entry.insert("1.0", '\n'.join(unicode_lines))

    def do_clear():
        """清空输入框和输出区域。"""
        entry.delete("1.0", tk.END)
        output.config(state="normal")
        output.delete("1.0", tk.END)
        output.edit_reset()
        cond_entry.delete(0, tk.END)

    # 取消/重做安全包装
    def _safe_undo(_event=None):
        try:
            output.edit_undo()
        except tk.TclError:
            pass
        return "break"

    def _safe_redo(_event=None):
        try:
            output.edit_redo()
        except tk.TclError:
            pass
        return "break"

    # --- 创建主窗口 ---
    root = tk.Tk()
    root.title("化学方程式配平工具")
    root.resizable(False, False)

    # --- 菜单栏 ---
    menubar = tk.Menu(root)
    root.config(menu=menubar)

    file_menu = tk.Menu(menubar, tearoff=0)
    file_menu.add_command(label="退出", command=root.destroy)
    menubar.add_cascade(label="文件", menu=file_menu)

    def _show_help() -> None:
        """弹出使用指南对话框。"""
        help_dlg = tk.Toplevel(root)
        help_dlg.title("使用指南")
        help_dlg.resizable(False, False)
        help_dlg.transient(root)

        help_dlg.update_idletasks()
        dhw, dhh = 520, 380
        phx = root.winfo_x() + (root.winfo_width() - dhw) // 2
        phy = root.winfo_y() + (root.winfo_height() - dhh) // 2
        help_dlg.geometry(f"{dhw}x{dhh}+{phx}+{phy}")

        guide = ScrolledText(help_dlg, font=("Microsoft YaHei", 10),
                              wrap="word", width=60, height=18)
        guide.pack(padx=12, pady=12, fill="both", expand=True)
        guide.tag_config("title", font=("Microsoft YaHei", 12, "bold"))
        guide.tag_config("section", font=("Microsoft YaHei", 10, "bold"))

        guide.insert(tk.END, "化学方程式配平工具 — 使用指南\n\n", "title")

        guide.insert(tk.END, "一、基本输入\n", "section")
        guide.insert(tk.END,
            "  格式：Fe2+ + MnO4- + H+ = Fe3+ + Mn2+ + H2O\n"
            "  支持 -> 或 = 作为分隔符，支持 (s)/(g)/(l)/(aq) 状态标记。\n"
            "  输入完成后按 Ctrl+Enter 即可配平。\n\n")

        guide.insert(tk.END, "二、批量配平\n", "section")
        guide.insert(tk.END,
            "  在输入框中每行写一个方程式，按 Ctrl+Enter 一次性配平全部。\n"
            "  各行结果之间以分隔线隔开。\n\n")

        guide.insert(tk.END, "三、工具栏\n", "section")
        guide.insert(tk.END,
            "  符号工具栏支持快速插入：\n"
            "   · ↑ (气体上升)、↓ (沉淀下降)、⇌ (可逆符号)\n"
            "   · [O] (氧化)、[H] (还原)\n"
            "   · 加热/点燃等反应条件请直接在【条件】输入框中书写，"
            "导出图片时自动显示在等号/箭头正上方（符合化学方程式书写规范）。\n"
            "   · ⇌ 按钮在结果框选中文本时替换 =/->/→ 为 ⇌。\n\n")

        guide.insert(tk.END, "四、自定义分子式\n", "section")
        guide.insert(tk.END,
            "  点击主窗口的【自定义分子式】按钮打开管理窗口：\n"
            "   · 输入分子式后点击【保存】，自动按字母排序存入本地。\n"
            "   · 在列表中点击选中，再点【插入】添加到输入框光标位置；\n"
            "     双击列表项也可直接插入。自动在非行首位置补加 '+' 号。\n"
            "   · 【插入等号】在光标处添加 ' = '。\n"
            "   · 【换行】在光标处换行，方便输入多个方程式。\n"
            "   · 【删除选中】从列表中移除不再需要的分子式。\n\n")

        guide.insert(tk.END, "五、导出图片\n", "section")
        guide.insert(tk.END,
            "  点击【导出为图片】将结果框内容渲染为 300 DPI PNG 图片。\n"
            "   · 支持 Times New Roman 字体 + Unicode 上下标。\n"
            "   · 导出前可预览，确认后再保存。\n\n")

        guide.insert(tk.END, "六、注意事项\n", "section")
        guide.insert(tk.END,
            "   · 化学式须严格区分大小写（Co 为钴，CO 为一氧化碳）。\n"
            "   · 输入框中请勿添加前置系数（如 2H2O），程序自动计算。\n"
            "   · 无解或歧义时，结果区会以红色字体提示错误。\n"
            "   · 结果框可编辑：Ctrl+Z 撤销、Ctrl+Y 重做。\n"
            "   · 【复制纯文本】将结果复制到系统剪贴板。\n"
            "   · 图片导出需要 Pillow 库（pip install Pillow）。\n\n")

        guide.insert(tk.END, "─" * 52 + "\n")
        guide.insert(tk.END, "作者：张曾继明\n")
        guide.insert(tk.END, "联系方式：1617528813@qq.com\n")

        guide.config(state="disabled")

        btn_close = tk.Button(help_dlg, text="关闭", width=10,
                               command=help_dlg.destroy)
        btn_close.pack(pady=(0, 10))

    help_menu = tk.Menu(menubar, tearoff=0)
    help_menu.add_command(label="使用指南", command=_show_help)
    menubar.add_cascade(label="帮助", menu=help_menu)

    # 居中显示（窗口加大以容纳工具栏和底部按钮）
    root.update_idletasks()
    win_w = 620
    win_h = 560
    screen_w = root.winfo_screenwidth()
    screen_h = root.winfo_screenheight()
    x = (screen_w - win_w) // 2
    y = (screen_h - win_h) // 2
    root.geometry(f"{win_w}x{win_h}+{x}+{y}")

    # --- 输入区域 ---
    lbl_input = tk.Label(root, text="请在每个新行输入一个化学方程式（Ctrl+Enter 配平）：")
    lbl_input.pack(anchor="w", padx=15, pady=(12, 4))

    entry = ScrolledText(root, width=55, height=5,
                         font=("Consolas", 11), wrap="none")
    entry.pack(padx=15, fill="x")
    entry.bind("<Control-Return>", lambda _event: do_balance())
    entry.focus_set()

    # --- 符号工具栏 + 条件输入 ---
    toolbar = tk.Frame(root)
    toolbar.pack(padx=15, pady=(4, 2), fill="x")

    tk.Label(toolbar, text="符号：", font=("Microsoft YaHei", 9)).pack(
        side="left", padx=(0, 2))

    _TOOLBAR_ITEMS = [
        ("↑",   '↑'),
        ("↓",   '↓'),
        ("⇌",   None),   # 特殊处理：选中替换 =/->/→ 为 ⇌
        ("[O]",  '[O]'),
        ("[H]",  '[H]'),
    ]   # △ 和"点燃"已移至条件输入框
    for label, sym in _TOOLBAR_ITEMS:
        if sym is None:
            btn = tk.Button(toolbar, text=label, width=3,
                            font=("Consolas", 9),
                            command=_handle_reversible)
        else:
            btn = tk.Button(toolbar, text=label, width=3,
                            font=("Consolas", 9),
                            command=lambda s=sym: _insert_symbol(s))
        btn.pack(side="left", padx=1)

    tk.Label(toolbar, text="  条件：",
             font=("Microsoft YaHei", 9)).pack(side="left", padx=(10, 2))
    cond_entry = tk.Entry(toolbar, width=16, font=("Microsoft YaHei", 10))
    cond_entry.pack(side="left")

    # --- 操作按钮 ---
    btn_frame = tk.Frame(root)
    btn_frame.pack(pady=8)

    btn_balance = tk.Button(btn_frame, text="配平", width=10,
                            command=do_balance)
    btn_balance.pack(side="left", padx=5)

    btn_clear = tk.Button(btn_frame, text="清空", width=10,
                          command=do_clear)
    btn_clear.pack(side="left", padx=5)

    btn_formula = tk.Button(btn_frame, text="自定义分子式", width=12,
                            command=lambda: manage_formulas(root, entry))
    btn_formula.pack(side="left", padx=5)

    # --- 输出区域（可编辑 + 支持撤销/重做） ---
    lbl_output = tk.Label(root, text="配平结果（可编辑，Ctrl+Z 撤销 / Ctrl+Y 重做）：")
    lbl_output.pack(anchor="w", padx=15, pady=(6, 4))

    output = ScrolledText(root, width=62, height=8,
                          font=("Consolas", 11), wrap="word",
                          undo=True)
    output.pack(padx=15, pady=(0, 4), fill="both", expand=True)
    output.bind("<Control-z>", _safe_undo)
    output.bind("<Control-y>", _safe_redo)

    # 配置红色错误标签
    output.tag_config("error", foreground="red")

    # --- 底部按钮 ---
    bottom_frame = tk.Frame(root)
    bottom_frame.pack(pady=(4, 10))

    tk.Button(bottom_frame, text="复制纯文本", width=14,
              command=_copy_plain_text).pack(side="left", padx=5)

    tk.Button(bottom_frame, text="导出为图片", width=14,
              command=_do_export).pack(side="left", padx=5)

    root.mainloop()


# ============================================================
# 入口
# ============================================================

if __name__ == '__main__':
    gui_main()
