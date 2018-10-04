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
x  = symbols('x',  real=True)
y  = symbols('y',  real=True)
dx = symbols('dx', real=True, positive=True)
dy = symbols('dy', real=True, positive=True)
t  = symbols('t',  real=True, positive=True)

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

def wave(xd = 0, yd = 0):
	return exp(I * (k * (x + xd) + l * (y + yd) - w * t))

freq_relation_expr = {}
if args.coord_type == 'plane':
	# Continuous linearized barotropic equations
	eqn = [
		[ diff(wave(), t), - f * wave(), H * diff(wave(), x) ],
		[ diff(wave(), t),   f * wave(), H * diff(wave(), y) ],
		[ diff(wave(), t),   H * diff(wave(), x), H * diff(wave(), y) ]
	]
	M = trigsimp(simplify(Matrix([
		[ eqn[0][0] / wave(), eqn[0][1] / wave(), eqn[0][2] / wave() ],
		[ eqn[1][1] / wave(), eqn[1][0] / wave(), eqn[1][2] / wave() ],
		[ eqn[2][1] / wave(), eqn[2][2] / wave(), eqn[2][0] / wave() ]
	])).expand(complex=True), method='fu')
	freq_relation_expr['\\section{Continuous case}'] = [ M, simplify(M.det()) ]
	# A-grid linearized barotropic equations
	eqn = [
		[ diff(wave(), t), - f * wave(), H * (wave(xd=dx) - wave(xd=-dx)) / (2 * dx) ],
		[ diff(wave(), t),   f * wave(), H * (wave(yd=dy) - wave(yd=-dy)) / (2 * dy) ],
		[ diff(wave(), t),   H * (wave(xd=dx) - wave(xd=-dx)) / (2 * dx), H * (wave(yd=dy) - wave(yd=-dy)) / (2 * dy) ]
	]
	M = trigsimp(simplify(Matrix([
		[ eqn[0][0] / wave(), eqn[0][1] / wave(), eqn[0][2] / wave() ],
		[ eqn[1][1] / wave(), eqn[1][0] / wave(), eqn[1][2] / wave() ],
		[ eqn[2][1] / wave(), eqn[2][2] / wave(), eqn[2][0] / wave() ]
	])).expand(complex=True), method='fu')
	freq_relation_expr['\\section{A-grid case}'] = [ M, simplify(M.det()) ]
	# C-grid linearized barotropic equations
	eqn = [
		[ diff(wave(xd=-dx/2), t), - f / 4 * (wave(xd=-dx,   yd=-dy/2) + wave(yd=-dy/2) + wave(xd=-dx,  yd=dy/2) + wave(yd=dy/2)), H * (wave() - wave(xd=-dx)) / dx ],
		[ diff(wave(yd=-dy/2), t),   f / 4 * (wave(xd=-dx/2, yd=-dy  ) + wave(xd=-dx/2) + wave(xd=dx/2, yd=-dy ) + wave(xd=dx/2)), H * (wave() - wave(yd=-dy)) / dy ],
		[ diff(wave(),         t),   H * (wave(xd=dx/2) - wave(xd=-dx/2)) / dx, H * (wave(yd=dy/2) - wave(yd=-dy/2)) / dy ]
	]
	M = trigsimp(simplify(Matrix([
		[ eqn[0][0] / wave(xd=-dx/2), eqn[0][1] / wave(xd=-dx/2), eqn[0][2] / wave(xd=-dx/2) ],
		[ eqn[1][1] / wave(yd=-dy/2), eqn[1][0] / wave(yd=-dy/2), eqn[1][2] / wave(yd=-dy/2) ],
		[ eqn[2][1] / wave(),         eqn[2][2] / wave(),         eqn[2][0] / wave() ]
	])).expand(complex=True), method='fu')
	# freq_relation_expr['\\section{C-grid case}'] = collect(trigsimp(simplify(M.det()).expand(complex=True), method='fu'), [f, k, l])
	freq_relation_expr['\\section{C-grid case}'] = [ M, simplify(M.det()) ]
elif args.coord_type == 'sphere':
	pass

if args.viewer == 'skim':
	viewer = '/Applications/Skim.app/Contents/MacOS/Skim'
	preview(str.join('\n', [
		title + str.join('\n', [latex(expression, mode='equation*',
			symbol_names={
				w: '\\omega',
				dx: '\\Delta{x}',
				dy: '\\Delta{y}'
			}) for expression in expressions]) for title, expressions in freq_relation_expr.items()]),
		output='pdf', viewer=viewer,
		preamble='''
			\\documentclass{article}
			\\usepackage[a4paper]{geometry}
			\\usepackage{amsmath,amsfonts}
			\\begin{document}
		''')
else:
	pprint(freq_relation_expr)
