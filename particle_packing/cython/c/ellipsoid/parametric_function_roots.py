from sympy import symbols, eye, collect
from sympy.matrices import Matrix
from sympy.utilities.codegen import codegen

b1, b2, b3 = symbols('b1, b2, b3', real=True)
a11, a12, a13, a21, a22, a23, a31, a32, a33 = symbols('a11, a12, a13, a21, a22, a23, a31, a32, a33', real=True)
L = symbols('L', real=True)

b = Matrix([b1, b2, b3])
A = Matrix([[a11, a12, a13],
            [a21, a22, a23],
            [a31, a32, a33]])

C = L * eye(3) + (1 - L) * A

p = collect((L * (1 - L) * b.T * C.adjugate() * b)[0].expand(), L)
p_prime = p.diff(L)

q = collect(C.det(), L)
q_prime = q.diff(L)

f = p / q

h = collect((p_prime * q - p * q_prime).expand(), L)

h0 = h.coeff(L, 0)
h1 = h.coeff(L, 1)
h2 = h.coeff(L, 2)
h3 = h.coeff(L, 3)
h4 = h.coeff(L, 4)
h5 = h.coeff(L, 5)
h6 = h.coeff(L, 6)

codegen([("f", f), ("h0", h0), ("h1", h1), ("h2", h2), ("h3", h3), ("h4", h4), ("h5", h5), ("h6", h6)], "C", "parametric_function_roots", to_files=True, header=False, empty=True)