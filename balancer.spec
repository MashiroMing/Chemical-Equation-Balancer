# -*- mode: python ; coding: utf-8 -*-

a = Analysis(
    ['balancer.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=[
        'tkinter',
        'tkinter.messagebox',
        'tkinter.scrolledtext',
        'tkinter.filedialog',
        # Pillow 图片导出
        'PIL',
        'PIL.Image',
        'PIL.ImageDraw',
        'PIL.ImageFont',
        'PIL.ImageTk',
        'PIL._imagingtk',
        'PIL._tkinter_finder',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        # 科学计算
        'numpy', 'scipy', 'pandas', 'matplotlib',
        # 图形与图像（保留 Pillow，仅排除不需要的）
        'cv2', 'cairo',
        # GUI 框架（非 tkinter）
        'wx', 'PyQt5', 'PySide2', 'PySide6', 'PyQt6', 'kivy',
        # 网络
        'urllib3', 'requests', 'aiohttp', 'xmlrpc',
        # 数据库
        'sqlite3', 'sqlalchemy', 'pymongo',
        # 邮件
        'smtplib', 'poplib', 'imaplib',
        # 其他大模块
        'curses', 'ctypes.test', 'distutils', 'setuptools',
        'test', 'unittest', 'doctest', 'pdb',
        'multiprocessing', 'concurrent.futures',
        'asyncio', 'argparse',
        'ssl', 'socket', 'ftplib',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=None,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=None)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='化学方程式配平工具',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=r'icon.ico',
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='化学方程式配平工具',
)
