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
from tkinter import messagebox
from tkinter.scrolledtext import ScrolledText


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
    启发式规则：单元素离子 → 数字作为电荷上标；多元素 → 数字作为下标，
    仅末尾的正负号上标。
    """
    # 分离前置系数
    m = re.match(r'(\d*)(.*)', term)
    coeff = m.group(1)
    formula = m.group(2)

    m = re.search(r'(\d*)([+-]+)$', formula)
    if not m:
        return coeff + _digits_to_sub(formula)

    charge_digits = m.group(1)
    charge_signs = m.group(2)
    body = formula[:m.start()]

    if charge_digits:
        # 数一数主体中的大写字母（即几种元素）
        upper_count = len(set(re.findall(r'[A-Z]', body)))
        if upper_count <= 1:
            # 单元素离子：Fe2+ → Fe²⁺, Mn2+ → Mn²⁺, S2- → S²⁻
            body_unicode = _digits_to_sub(body)
            charge_unicode = _digits_to_super(charge_digits) + _signs_to_super(charge_signs)
        else:
            # 多元素：MnO4- → MnO₄⁻, Cr2O72- → Cr₂O₇₂⁻
            body_unicode = _digits_to_sub(body + charge_digits)
            charge_unicode = _signs_to_super(charge_signs)
    else:
        # 仅有正负号：H+ → H⁺
        body_unicode = _digits_to_sub(body)
        charge_unicode = _signs_to_super(charge_signs)

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

    def do_balance():
        """执行配平并显示结果。支持多行批量输入。"""
        raw = entry.get("1.0", tk.END).strip()
        if not raw:
            messagebox.showwarning("提示", "请输入方程式")
            return

        # 输出区域清空旧内容
        output.delete("1.0", tk.END)

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
        output.delete("1.0", tk.END)

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

        guide.insert(tk.END, "三、自定义分子式\n", "section")
        guide.insert(tk.END,
            "  点击主窗口的【自定义分子式】按钮打开管理窗口：\n"
            "   · 输入分子式后点击【保存】，自动按字母排序存入本地。\n"
            "   · 在列表中点击选中，再点【插入】添加到输入框光标位置；\n"
            "     双击列表项也可直接插入。自动在非行首位置补加 '+' 号。\n"
            "   · 【插入等号】在光标处添加 ' = '。\n"
            "   · 【换行】在光标处换行，方便输入多个方程式。\n"
            "   · 【删除选中】从列表中移除不再需要的分子式。\n\n")

        guide.insert(tk.END, "四、注意事项\n", "section")
        guide.insert(tk.END,
            "   · 化学式须严格区分大小写（Co 为钴，CO 为一氧化碳）。\n"
            "   · 输入框中请勿添加前置系数（如 2H2O），程序自动计算。\n"
            "   · 无解或歧义时，结果区会以红色字体提示错误。\n"
            "   · 分子式列表保存在同目录的 custom_formulas.json 文件中。\n\n")

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

    # 居中显示
    root.update_idletasks()
    win_w = 560
    win_h = 420
    screen_w = root.winfo_screenwidth()
    screen_h = root.winfo_screenheight()
    x = (screen_w - win_w) // 2
    y = (screen_h - win_h) // 2
    root.geometry(f"{win_w}x{win_h}+{x}+{y}")

    # --- 输入区域 ---
    lbl_input = tk.Label(root, text="请在每个新行输入一个化学方程式（Ctrl+Enter 配平）：")
    lbl_input.pack(anchor="w", padx=15, pady=(12, 4))

    entry = ScrolledText(root, width=55, height=6,
                         font=("Consolas", 11), wrap="none")
    entry.pack(padx=15, fill="x")
    entry.bind("<Control-Return>", lambda _event: do_balance())
    entry.focus_set()

    # --- 按钮区域 ---
    btn_frame = tk.Frame(root)
    btn_frame.pack(pady=10)

    btn_balance = tk.Button(btn_frame, text="配平", width=10,
                            command=do_balance)
    btn_balance.pack(side="left", padx=5)

    btn_clear = tk.Button(btn_frame, text="清空", width=10,
                          command=do_clear)
    btn_clear.pack(side="left", padx=5)

    btn_formula = tk.Button(btn_frame, text="自定义分子式", width=12,
                            command=lambda: manage_formulas(root, entry))
    btn_formula.pack(side="left", padx=5)

    # --- 输出区域 ---
    lbl_output = tk.Label(root, text="配平结果：")
    lbl_output.pack(anchor="w", padx=15, pady=(4, 4))

    output = ScrolledText(root, width=60, height=6,
                          font=("Consolas", 11), wrap="word")
    output.pack(padx=15, pady=(0, 12), fill="both", expand=True)

    # 配置红色错误标签
    output.tag_config("error", foreground="red")

    root.mainloop()


# ============================================================
# 入口
# ============================================================

if __name__ == '__main__':
    gui_main()
