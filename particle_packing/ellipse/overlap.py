import numpy as np

def overlap_potential_py(rA, radiiA, phiA, rB, radiiB, phiB):
    
    """
    
    X = Q^T * O^-2 * Q, where Q is rotation matrix given by angle phi and
    O is diagonal matrix with radii values along main diagonal
    
    What's needed is XA^-1 = Q^T * O^2 * Q
    
    """
    
    
    O = np.diag(radiiA ** 2)
    
    Q = np.array([[np.cos(phiA), np.sin(phiA)],
              [-np.sin(phiA), np.cos(phiA)]])
    
    XA = Q.T.dot(O).dot(Q)
    

    
    
    """
    
    X = Q^T * O^-2 * Q, where Q is rotation matrix given by angle phi and
    O is diagonal matrix with radii values along main diagonal
    
    What's needed is XB^1/2 = Q^T * O^-1 * Q
    
    """
    
    
    O = np.diag(radiiB ** -1)
    
    Q = np.array([[np.cos(phiB), np.sin(phiB)],
              [-np.sin(phiB), np.cos(phiB)]])
    
    XB = Q.T.dot(O).dot(Q)
    
    
    # Find rAB
    rAB = rB - rA
    
    # Find A_AB
    A_AB = XB.dot(XA).dot(XB)
    
    # Find a_AB
    a_AB = XB.dot(rAB)
    
    
    coeffs = h_coeffs(A_AB, a_AB)
    roots = np.roots(coeffs)
    
    F = None
    for r in roots:
        if np.isreal(r) and (0. <= r <= 1.):
            #print 'max in (0, 1)'
            F =  f_AB(A_AB, a_AB, np.real(r))
    
    if F is None:
        if f_AB(A_AB, a_AB, 0) > f_AB(A_AB, a_AB, 1):
            #print 'max at lambda = 0'
            F = f_AB(A_AB, a_AB, 0)
        else:
            #print 'max at lambda = 1'
            F = f_AB(A_AB, a_AB, 1)

    
    return F



def p_AB(A, b, L):
    a11 = A[0, 0]
    a12 = A[0, 1]
    a21 = A[1, 0]
    a22 = A[1,1]
    b1 = b[0, 0]
    b2 = b[1, 0]
    
    return L**3*(a11*b2**2 - a12*b1*b2 - a21*b1*b2 + a22*b1**2 - b1**2 - b2**2) + L**2*(-2*a11*b2**2 + 2*a12*b1*b2 + 2*a21*b1*b2 - 2*a22*b1**2 + b1**2 + b2**2) + L*(a11*b2**2 - a12*b1*b2 - a21*b1*b2 + a22*b1**2)

def q_AB(A, b, L):
    a11 = A[0, 0]
    a12 = A[0, 1]
    a21 = A[1, 0]
    a22 = A[1,1]
    b1 = b[0, 0]
    b2 = b[1, 0]
    
    return L**2*(a11*a22 - a11 - a12*a21 - a22 + 1) + L*(-2*a11*a22 + a11 + 2*a12*a21 + a22) + a11*a22 - a12*a21

def f_AB(A, b, L):
    return p_AB(A, b, L) / q_AB(A, b, L)

def h_coeffs(A, b):
    a11 = A[0, 0]
    a12 = A[0, 1]
    a21 = A[1, 0]
    a22 = A[1,1]
    b1 = b[0, 0]
    b2 = b[1, 0]
    
    coeff = [a11**2*a22*b2**2 - a11**2*b2**2 - a11*a12*a21*b2**2 - a11*a12*a22*b1*b2 + a11*a12*b1*b2 - a11*a21*a22*b1*b2 + a11*a21*b1*b2 + a11*a22**2*b1**2 - 2*a11*a22*b1**2 - 2*a11*a22*b2**2 + a11*b1**2 + 2*a11*b2**2 + a12**2*a21*b1*b2 + a12*a21**2*b1*b2 - a12*a21*a22*b1**2 + a12*a21*b1**2 + a12*a21*b2**2 + a12*a22*b1*b2 - a12*b1*b2 + a21*a22*b1*b2 - a21*b1*b2 - a22**2*b1**2 + 2*a22*b1**2 + a22*b2**2 - b1**2 - b2**2,
             -4*a11**2*a22*b2**2 + 2*a11**2*b2**2 + 4*a11*a12*a21*b2**2 + 4*a11*a12*a22*b1*b2 - 2*a11*a12*b1*b2 + 4*a11*a21*a22*b1*b2 - 2*a11*a21*b1*b2 - 4*a11*a22**2*b1**2 + 6*a11*a22*b1**2 + 6*a11*a22*b2**2 - 2*a11*b1**2 - 2*a11*b2**2 - 4*a12**2*a21*b1*b2 - 4*a12*a21**2*b1*b2 + 4*a12*a21*a22*b1**2 - 4*a12*a21*b1**2 - 4*a12*a21*b2**2 - 2*a12*a22*b1*b2 - 2*a21*a22*b1*b2 + 2*a22**2*b1**2 - 2*a22*b1**2 - 2*a22*b2**2, 
             6*a11**2*a22*b2**2 - a11**2*b2**2 - 6*a11*a12*a21*b2**2 - 6*a11*a12*a22*b1*b2 + a11*a12*b1*b2 - 6*a11*a21*a22*b1*b2 + a11*a21*b1*b2 + 6*a11*a22**2*b1**2 - 6*a11*a22*b1**2 - 6*a11*a22*b2**2 + a11*b1**2 + 6*a12**2*a21*b1*b2 + 6*a12*a21**2*b1*b2 - 6*a12*a21*a22*b1**2 + 5*a12*a21*b1**2 + 5*a12*a21*b2**2 + a12*a22*b1*b2 + a12*b1*b2 + a21*a22*b1*b2 + a21*b1*b2 - a22**2*b1**2 + a22*b2**2, 
             -4*a11**2*a22*b2**2 + 4*a11*a12*a21*b2**2 + 4*a11*a12*a22*b1*b2 + 4*a11*a21*a22*b1*b2 - 4*a11*a22**2*b1**2 + 2*a11*a22*b1**2 + 2*a11*a22*b2**2 - 4*a12**2*a21*b1*b2 - 4*a12*a21**2*b1*b2 + 4*a12*a21*a22*b1**2 - 2*a12*a21*b1**2 - 2*a12*a21*b2**2, 
             a11**2*a22*b2**2 - a11*a12*a21*b2**2 - a11*a12*a22*b1*b2 - a11*a21*a22*b1*b2 + a11*a22**2*b1**2 + a12**2*a21*b1*b2 + a12*a21**2*b1*b2 - a12*a21*a22*b1**2]
    return coeff
