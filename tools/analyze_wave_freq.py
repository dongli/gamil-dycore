#!/usr/bin/env python3

from argparse import *
from sympy import *
import sys

parser = ArgumentParser(
	description='Barotropic wave frequency analysis codes.',
	formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-c', '--coord-type', dest='coord_type',
	help='Choose coordinate system type.', choices=['plane', 'sphere'], default='plane')
parser.add_argument('--viewer',
	help='Choose viewer to display result.', choices=['pdf'])
args = parser.parse_args()

# Define symbols.
f  = symbols('f',  real=True)
w  = symbols('w',  real=True, positive=True)
k  = symbols('k',  real=True, positive=True, integer=True)
l  = symbols('l',  real=True, positive=True, integer=True)
H  = symbols('H',  real=True, positive=True)
x  = symbols('x',  real=True)
y  = symbols('y',  real=True)
dx = symbols('dx', real=True, positive=True)
dy = symbols('dy', real=True, positive=True)
t  = symbols('t',  real=True, positive=True)
a  = symbols('a',  real=True, positive=True)

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
	va = (wave(xd=-dx,   yd=-dy/2) + wave(yd=-dy/2) + wave(xd=-dx,  yd=dy/2) + wave(yd=dy/2)) / 4
	ua = (wave(xd=-dx/2, yd=-dy  ) + wave(xd=-dx/2) + wave(xd=dx/2, yd=-dy ) + wave(xd=dx/2)) / 4
	eqn = [
		[ diff(wave(xd=-dx/2), t), - f * va, H * (wave() - wave(xd=-dx)) / dx ],
		[ diff(wave(yd=-dy/2), t),   f * ua, H * (wave() - wave(yd=-dy)) / dy ],
		[ diff(wave(),         t),   H * (wave(xd=dx/2) - wave(xd=-dx/2)) / dx, H * (wave(yd=dy/2) - wave(yd=-dy/2)) / dy ]
	]
	M = trigsimp(simplify(Matrix([
		[ eqn[0][0] / wave(xd=-dx/2), eqn[0][1] / wave(xd=-dx/2), eqn[0][2] / wave(xd=-dx/2) ],
		[ eqn[1][1] / wave(yd=-dy/2), eqn[1][0] / wave(yd=-dy/2), eqn[1][2] / wave(yd=-dy/2) ],
		[ eqn[2][1] / wave(),         eqn[2][2] / wave(),         eqn[2][0] / wave() ]
	])).expand(complex=True), method='fu')
	freq_relation_expr['\\section{C-grid case}'] = [ M, simplify(M.det()) ]

	# C-grid linearized barotropic equations with modified Coriolis terms
	va = (
		(1 - a) * (a * (wave(xd=-dx,  yd=-3*dy/2) + wave(yd=-3*dy/2)) + (1 - a) * (wave(xd=-2*dx, yd=-3*dy/2) + wave(xd=dx, yd=-3*dy/2))) +
		     a  * (a * (wave(xd=-dx,  yd=  -dy/2) + wave(yd=  -dy/2)) + (1 - a) * (wave(xd=-2*dx, yd=  -dy/2) + wave(xd=dx, yd=  -dy/2))) +
		     a  * (a * (wave(xd=-dx,  yd=   dy/2) + wave(yd=   dy/2)) + (1 - a) * (wave(xd=-2*dx, yd=   dy/2) + wave(xd=dx, yd=   dy/2))) +
		(1 - a) * (a * (wave(xd=-dx,  yd= 3*dy/2) + wave(yd= 3*dy/2)) + (1 - a) * (wave(xd=-2*dx, yd= 3*dy/2) + wave(xd=dx, yd= 3*dy/2)))
	)
	ua = (
		(1 - a) * (a * (wave(xd=-dx/2, yd=-2*dy) + wave(xd=dx/2, yd=-2*dy)) + (1 - a) * (wave(xd=-3*dx/2, yd=-2*dy) + wave(xd=3*dx/2, yd=-2*dy))) +
		     a  * (a * (wave(xd=-dx/2, yd=  -dy) + wave(xd=dx/2, yd=  -dy)) + (1 - a) * (wave(xd=-3*dx/2, yd=  -dy) + wave(xd=3*dx/2, yd=  -dy))) +
		     a  * (a * (wave(xd=-dx/2          ) + wave(xd=dx/2          )) + (1 - a) * (wave(xd=-3*dx/2          ) + wave(xd=3*dx/2          ))) +
		(1 - a) * (a * (wave(xd=-dx/2, yd=   dy) + wave(xd=dx/2, yd=   dy)) + (1 - a) * (wave(xd=-3*dx/2, yd=   dy) + wave(xd=3*dx/2, yd=   dy)))
	)
	eqn = [
		[ diff(wave(xd=-dx/2), t), - f * va, H * (wave() - wave(xd=-dx)) / dx ],
		[ diff(wave(yd=-dy/2), t),   f * ua, H * (wave() - wave(yd=-dy)) / dy ],
		[ diff(wave(),         t),   H * (wave(xd=dx/2) - wave(xd=-dx/2)) / dx, H * (wave(yd=dy/2) - wave(yd=-dy/2)) / dy ]
	]
	M = trigsimp(simplify(Matrix([
		[ eqn[0][0] / wave(xd=-dx/2), eqn[0][1] / wave(xd=-dx/2), eqn[0][2] / wave(xd=-dx/2) ],
		[ eqn[1][1] / wave(yd=-dy/2), eqn[1][0] / wave(yd=-dy/2), eqn[1][2] / wave(yd=-dy/2) ],
		[ eqn[2][1] / wave(),         eqn[2][2] / wave(),         eqn[2][0] / wave() ]
	])).expand(complex=True), method='fu')
	freq_relation_expr['\\section{Modified C-grid case}'] = [ M, simplify(M.det()) ]
elif args.coord_type == 'sphere':
	pass

if args.viewer == 'pdf':
	if sys.platform == 'darwin':
		viewer = '/Applications/Skim.app/Contents/MacOS/Skim'
	elif sys.platform == 'linux':
		viewer = 'evince'
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
			\\usepackage[a4paper, paperwidth=60cm, margin=2cm]{geometry}
			\\usepackage{amsmath,amsfonts}
			\\begin{document}
		''')
else:
	pprint(freq_relation_expr)
