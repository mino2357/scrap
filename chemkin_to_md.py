import re
import sys
from pathlib import Path


def format_formula(formula: str) -> str:
    """Convert a CHEMKIN reaction string to LaTeX-friendly form."""
    # determine arrow type
    if '<=>' in formula or ('=' in formula and '=>' not in formula):
        arrow = '\\rightleftharpoons'
        formula = formula.replace('<=>', ' -> ').replace('=', ' -> ')
    elif '=>' in formula:
        arrow = '\\rightarrow'
        formula = formula.replace('=>', ' -> ')
    else:
        arrow = '\\rightarrow'

    # insert spaces around + and parentheses
    formula = formula.replace('+', ' + ')
    formula = formula.replace('(', ' (')
    formula = re.sub(r'\s+', ' ', formula).strip()
    # tidy spaces inside third-body parentheses: '( + M)' -> '(+M)'
    formula = re.sub(r'\( \+ ([^\)]+)\)', r'(+\1)', formula)
    formula = formula.replace('->', arrow)
    # replace digits that follow element symbols or ')' with subscripts,
    # leaving leading stoichiometric coefficients untouched
    formula = re.sub(r'(?<=[A-Za-z)])(\d+)', r'_{\1}', formula)
    return formula


def parse_reactions(text: str):
    lines = text.splitlines()
    reactions = []
    in_section = False
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not in_section:
            if line.upper() == 'REACTIONS':
                in_section = True
            i += 1
            continue
        if line.upper() == 'END':
            break
        if '=' in line and not line.startswith('!'):
            # grab the reaction formula before the rate coefficients.
            # Coefficients are typically separated from the formula by two or
            # more spaces, while single spaces may exist inside the formula.
            line_no_comment = line.split('!')[0].rstrip()
            formula = re.split(r'\s{2,}', line_no_comment)[0].strip()
            dup = False
            j = i + 1
            while j < len(lines) and lines[j].strip() == '':
                j += 1
            if j < len(lines) and lines[j].strip().upper().startswith('DUP'):
                dup = True
                i = j  # skip the DUP line
            formatted = format_formula(formula)
            reactions.append((formatted, dup))
        i += 1
    return reactions


def convert_file(path: Path) -> str:
    content = path.read_text()
    reactions = parse_reactions(content)
    md_lines = []
    md_lines.append('以下は、`{}` の **REACTIONS** セクションに記載された化学反応式を Markdown 形式で整理した一覧です\n'.format(path))
    for idx, (formula, dup) in enumerate(reactions, start=1):
        line = f"{idx}. $\\mathrm{{{formula}}}$"
        if dup:
            line += " *(DUP)*"
        md_lines.append(line)
    md_lines.append('')
    return '\n'.join(md_lines)


if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print('Usage: python chemkin_to_md.py <chem.inp> [output.md]')
        sys.exit(1)
    inp = Path(sys.argv[1])
    md_text = convert_file(inp)
    if len(sys.argv) == 3:
        Path(sys.argv[2]).write_text(md_text)
    else:
        print(md_text)
