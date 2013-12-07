from sympy import symbols, eye, collect
from sympy.matrices import Matrix
from sympy.utilities.codegen import codegen

b1, b2 = symbols('b1, b2', real=True)
a11, a12, a21, a22 = symbols('a11, a12, a21, a22', real=True)
L = symbols('L', real=True)

b = Matrix([b1, b2])
A = Matrix([[a11, a12],
            [a21, a22]])

C = L * eye(2) + (1 - L) * A

p = collect((L * (1 - L) * b.T * C.adjugate() * b)[0].expand(), L)
p_prime = p.diff(L)

q = collect(C.det(), L)
q_prime = q.diff(L)

h = collect((p_prime * q - p * q_prime).expand(), L)

f = p / q

h0 = h.coeff(L, 0)
h1 = h.coeff(L, 1)
h2 = h.coeff(L, 2)
h3 = h.coeff(L, 3)
h4 = h.coeff(L, 4)

codegen([("f", f), ("h0", h0), ("h1", h1), ("h2", h2), ("h3", h3), ("h4", h4)], "C", "parametric_function_roots", to_files=True, header=False, empty=True)