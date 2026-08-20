import pathlib
import sys


if len(sys.argv) != 3:
    raise SystemExit("usage: review_v9_python_runner.py <script.py> <isolated_workdir>")

script = pathlib.Path(sys.argv[1]).resolve()
workdir = pathlib.Path(sys.argv[2]).resolve()
code = script.read_text(encoding="utf-8")
needle = 'MR = "D:/R_ex/MR"'
if code.count(needle) != 1:
    raise SystemExit(f"expected exactly one hard-coded MR assignment in {script.name}")
code = code.replace(needle, f"MR = {str(workdir)!r}")
exec(compile(code, str(script), "exec"), {"__name__": "__main__"})
