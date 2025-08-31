import re
import subprocess
import sys
from pathlib import Path
import shutil


def parse_reactions(text: str):
    lines = text.splitlines()
    edges = set()
    in_section = False
    for line in lines:
        line = line.strip()
        if not in_section:
            if line.upper() == 'REACTIONS':
                in_section = True
            continue
        if line.upper() == 'END':
            break
        if '=' not in line or line.startswith('!'):
            continue
        line_no_comment = line.split('!')[0].rstrip()
        formula = re.split(r"\s{2,}", line_no_comment)[0].strip()
        reversible = False
        if '<=>' in formula or ('=' in formula and '=>' not in formula):
            lhs, rhs = re.split(r'<=>|=', formula)
            reversible = True
        elif '=>' in formula:
            lhs, rhs = formula.split('=>')
        elif '->' in formula:
            lhs, rhs = formula.split('->')
        else:
            continue
        reactants = split_species(lhs)
        products = split_species(rhs)
        for r in reactants:
            for p in products:
                edges.add((r, p))
                if reversible:
                    edges.add((p, r))
    return edges


def split_species(side: str):
    # remove third-body notation like (+M), (+N2)
    side = re.sub(r'\(\+[^)]+\)', '', side)
    species = []
    for token in side.split('+'):
        token = token.strip()
        if not token:
            continue
        token = re.sub(r'^\d+', '', token).strip()
        if not token or token == 'M':
            continue
        species.append(token)
    return species


def build_dot(edges) -> str:
    lines = ['digraph G {']
    for src, dst in sorted(edges):
        lines.append(f'    "{src}" -> "{dst}";')
    lines.append('}')
    return '\n'.join(lines)


def convert_file(inp: Path, out: Path):
    content = inp.read_text()
    edges = parse_reactions(content)
    dot = build_dot(edges)
    suffix = out.suffix.lower()
    if suffix == '.dot':
        out.write_text(dot)
    elif suffix == '.png':
        dot_cmd = shutil.which('dot')
        if not dot_cmd:
            raise RuntimeError('Graphviz "dot" executable not found; install graphviz to output PNG')
        subprocess.run([dot_cmd, '-Tpng', '-o', str(out)], input=dot.encode(), check=True)
    else:
        raise ValueError('Output file must end with .dot or .png')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python chemkin_to_graph.py <chem.inp> <output.{dot,png}>')
        sys.exit(1)
    convert_file(Path(sys.argv[1]), Path(sys.argv[2]))