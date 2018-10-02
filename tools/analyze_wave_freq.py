#!/usr/bin/env python3

from argparse import *
from sympy import *

parser = ArgumentParser(
	description='Barotropic wave frequency analysis codes.',
	formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-c', '--coord-type', dest='coord_type',
	help='Choose coordinate system type.', choices=['plane', 'sphere'], default='plane')
parser.add_argument('--viewer',
	help='Choose viewer to display result.', choices=['skim'])
args = parser.parse_args()

# Define symbols.
f  = symbols('f',  real=True)
w  = symbols('w',  real=True, positive=True)
k  = symbols('k',  real=True, positive=True, integer=True)
l  = symbols('l',  real=True, positive=True, integer=True)
H  = symbols('H',  real=True, positive=True)
a  = symbols('a',  real=True, positive=True)
dx = symbols('dx', real=True, positive=True)
dy = symbols('dy', real=True, positive=True)

# Artifacts
skx    = sin(k * dx) / (k * dx)
ckx    = cos(k * dx) / (k * dx)
sky    = sin(k * dy) / (k * dy)
cky    = cos(k * dy) / (k * dy)
slx    = sin(l * dx) / (l * dx)
clx    = cos(l * dx) / (l * dx)
sly    = sin(l * dy) / (l * dy)
cly    = cos(l * dy) / (l * dy)
skx2   = sin(k * dx / 2) / (k * dx / 2)
ckx2   = cos(k * dx / 2) / (k * dx / 2)
sky2   = sin(k * dy / 2) / (k * dy / 2)
cky2   = cos(k * dy / 2) / (k * dy / 2)
slx2   = sin(l * dx / 2) / (l * dx / 2)
clx2   = cos(l * dx / 2) / (l * dx / 2)
sly2   = sin(l * dy / 2) / (l * dy / 2)
cly2   = cos(l * dy / 2) / (l * dy / 2)
ckxly2 = cos(k * dx / 2) * cos(l * dy / 2)

freq_relation_expr = {}
if args.coord_type == 'plane':
	M = Matrix(
		[
			[            -I * w,                -f,        I * k ],
			[                 f,            -I * w,        I * l ],
			[         I * k * H,         I * l * H,       -I * w ]
		])
	freq_relation_expr['\\section{Continuous case}'] = simplify(M.det())
	M = Matrix(
		[
			[            -I * w,                -f,  skx * I * k ],
			[                 f,            -I * w,  sly * I * l ],
			[   skx * I * k * H,   sly * I * l * H,       -I * w ]
		])
	freq_relation_expr['\\section{A-grid case}'] = simplify(M.det())
	M = Matrix(
		[
			[            -I * w,       -ckxly2 * f, skx2 * I * k ],
			[       -ckxly2 * f,            -I * w, sly2 * I * l ],
			[  skx2 * I * k * H,  sly2 * I * l * H,       -I * w ]
		])
	freq_relation_expr['\\section{C-grid case}'] = simplify(M.det())
elif args.coord_type == 'sphere':
	pass

if args.viewer == 'skim':
	viewer = '/Applications/Skim.app/Contents/MacOS/Skim'
	preview(str.join('\n', [
		title + latex(expr, mode='equation*',
			symbol_names={
				w: '\\omega',
				dx: '\\Delta{x}',
				dy: '\\Delta{y}'
			}) for title, expr in freq_relation_expr.items()]),
		output='pdf', viewer=viewer,
		preamble='''
			\\documentclass{article}
			\\usepackage[a4paper, landscape]{geometry}
			\\usepackage{amsmath,amsfonts}
			\\begin{document}
		''')
else:
	pprint(freq_relation_expr)
