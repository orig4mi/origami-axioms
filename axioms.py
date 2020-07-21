"""
Created on Sat Mar 21 09:39:08 2020

@author: Jorge C. Lucero

Simulation of origami axioms in the plane
"""

from math import sqrt, cos, acos, isclose, pi, sinh, asinh, cosh, acosh,\
                 copysign
from sys import float_info

TOL = float_info.epsilon


class Point(object):

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return '(' + str(self.x) + ',' + str(self.y) + ')'

    def __repr__(self):
        return 'Point(' + str(self.x) + ',' + str(self.y) + ')'

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def __truediv__(self, k):
        return Point(self.x/k, self.y/k)

    def __mul__(self, k):
        return Point(self.x*k, self.y*k)

    def __matmul__(self, other):
        return self.x*other.x + self.y*other.y

    def __eq__(self, other):
        return isclose(self.x, other.x) and isclose(self.y, other.y)

    def __neg__(self):
        return Point(-self.x, -self.y)

    def getCoords(self):
        return self.x, self.y

    def getNormal(self):
        return Point(self.y, -self.x)

    def norm(self):
        return sqrt(self @ self)

    def dist(self, other):
        return sqrt((self - other) @ (self - other))

    def sqNorm(self):
        return self @ self

    def getProp(self, other):
        return (self @ other)/self.sqNorm()

    def vectProd(self, other):
        return self.x*other.y - self.y*other.x

    def isParallel(self, other):
        return isSmall(self.vectProd(other))

    def isNormal(self, other):
        return isSmall(self @ other)


class Line(object):

    def __init__(self, a, b, c):

        if isSmall(a) and isSmall(b):
            raise ValueError('Line has a null normal vector')

        self.a = a
        self.b = b
        self.c = c

    def __str__(self):
        if self.b >= 0:
            str_b = '+ ' + str(self.b)
        else:
            str_b = '- ' + str(-self.b)

        if self.c >= 0:
            str_c = '+ ' + str(self.c)
        else:
            str_c = '- ' + str(-self.c)

        return str(self.a) + 'x ' + str_b + 'y ' + str_c + ' = 0'

    def __repr__(self):
        return ('Line(' + str(self.a) + ', ' + str(self.b) + ', '
                + str(self.c) + ')')

    def __eq__(self, other):
        k = self.getNormalProp(other)
        return (self.isParallel(other) and isclose(self.c, other.c*k))

    def __contains__(self, point):
        return isSmall(self.a*point.x + self.b*point.y + self.c)

    def vectProd(self, other):
        return self.getNormal().vectProd(other.getNormal())

    def getNormalProp(self, other):
        return self.getNormal().getProp(other.getNormal())

    def getCoefs(self):
        return self.a, self.b, self.c

    def getNormal(self):
        return Point(self.a, self.b)

    def getParallel(self):
        return self.getNormal().getNormal()

    def isParallel(self, other):
        return isSmall(self.vectProd(other))

    def intersection(self, other):
        if self == other:
            print('Lines are coincident')
            return None
        if self.isParallel(other):
            print('Lines do not intersect')
            return None
        else:
            d = self.vectProd(other)
            x, y = ((self.b*other.c - self.c*other.b)/d,
                    (self.c*other.a - self.a*other.c)/d)
            return Point(x, y)

    def dist(self, other):
        if not self.isParallel(other):
            return 0
        else:
            k = self.getNormalProp(other)
            return abs(self.c - other.c*k)/self.getNormal().norm()


def cubicRoot(x):
    '''
    Cubic root with sign
    '''
    if x >= 0:
        return x**(1/3)
    elif x < 0:
        return -(abs(x)**(1/3))


def makeLine(normal, x0):
    '''
    Create a line from a normal vector and a point
    '''
    return Line(normal.x, normal.y, -normal @ x0)


def solveCubic(a0, b0, c0, d0):
    '''
    Get the real roots of a cubic equation by trigonometric/hyperbolic methods.
    Source: Wikipedia
    '''
    b = b0/a0
    c = c0/a0
    d = d0/a0

    tau = 2*pi
    b13 = b/3

    p = (3*c - b*b)/3
    q = (2*b*b*b - 9*b*c + 27*d)/27
    d = 4*p*p*p + 27*q*q

    if isSmall(d):
        if isSmall(p):
            return [-b13]
        else:
            t = 3*q/p
            return [t - b13, -t/2 - b13]
    elif d < 0:
        u = sqrt(-4*p/3)
        k = acos(-4*q/(u*u*u))
        return[u*cos(k/3) - b13, u*cos((k - tau)/3) - b13,
               u*cos((k + tau)/3) - b13]
    else:
        if isSmall(p):
            return [-cubicRoot(q) - b13]
        elif p < 0:
            u = copysign(sqrt(-4*p/3), -q)
            return [u*cosh(acosh(-4*q/(u*u*u))/3) - b13]
        else:
            u = -sqrt(4*p/3)
            return [u*sinh(asinh(-4*q/(u*u*u))/3) - b13]


def isSmall(x):
    return isclose(x, 0, abs_tol=TOL)


def axiom_1(p, q):
    '''
    Given two distinct points p and q, fold to place p onto q.
    '''
    if p == q:
        raise ValueError('axiom_1: points are coincident')

    return [makeLine(p - q, (p + q)/2)]


def axiom_2(m, n):
    '''
    Given two distinct lines m and n, fold to align m and n.
    '''
    if m == n:
        raise ValueError('axiom_2: lines are coincident')

    if m.isParallel(n):
        k = m.getNormalProp(n)
        return [Line(m.a, m.b, (m.c + n.c*k)/2)]

    x0 = m.intersection(n)
    normal_m = m.getNormal()
    normal_n = n.getNormal()
    normal = normal_n*normal_m.norm() + normal_m*normal_n.norm()
    return [makeLine(normal, x0), makeLine(normal.getNormal(), x0)]


def axiom_3(m):
    '''
    Fold along a given line m
    '''
    return m


def axiom_4(p, q):
    '''
    Given two distinct points p and q, fold along a line passing through p
    and q.
    '''
    if p == q:
        raise ValueError('axiom_4: points are coincident')

    return [makeLine((p - q).getNormal(), q)]


def axiom_5(p, m):
    '''
    Given a line m and a point p, fold along a line passing through p to
    reflect half of m onto its other half.
    '''
    return [makeLine(m.getParallel(), p)]


def axiom_6(p, q, m):
    '''
    Given a line m, a point p not on m and a point q, fold along a line
    passing through q to place p onto m.
    '''
    if p in m:
        raise ValueError('axiom_6: point p is in line m')

    A = p.dist(q)
    B = m.a*q.x + m.b*q.y + m.c
    M = m.getNormal().sqNorm()
    disc = M*A*A - B*B

    if isSmall(disc):
        r = q - m.getNormal()*B/M
        return [makeLine(p - r, q)]
    elif disc < 0:
        return []
    else:
        r1 = q - (m.getNormal()*B + m.getParallel()*sqrt(disc))/M
        r2 = q - (m.getNormal()*B - m.getParallel()*sqrt(disc))/M
        return [makeLine(p - r1, q), makeLine(p - r2, q)]


def axiom_7(p, q, m, n):
    '''
    Given two lines m and n, a point p not on m, and a point q not on n, where
    m and n are distinct or p and q are distinct, fold to place p onto m, and
    q onto n.
    '''

    if p in m:
        raise ValueError('axiom_7: point p is in line m')
    if q in n:
        raise ValueError('axiom_7: point q is in line n')
    if (m == n and p == q):
        raise ValueError('axiom_7: both points and lines are coincident')

    solution = []

    if m.isParallel(n):
        A = p.dist(q)
        k = m.getNormalProp(n)
        B = n.c*k - m.c
        N = m.getNormal()
        M = N.sqNorm()
        disc = M*A*A - B*B
        if isSmall(disc):
            vec = Point(q.x - p.x - m.a*B/M, q.y - p.y - m.b*B/M)
            o = makeLine(vec, p)
            solution += axiom_1(p, o.intersection(m))
#        elif disc < 0:
#            return solution
        elif disc > 0:
            solution = []
            vec1 = p - q + (N*B + m.getParallel()*sqrt(disc))/M
            vec2 = p - q + (N*B - m.getParallel()*sqrt(disc))/M
            if not vec1.isParallel(N):
                solution += axiom_1(p, makeLine(vec1, p).intersection(m))
            if not vec2.isParallel(N):
                solution += axiom_1(p, makeLine(vec2, p).intersection(m))
#            return solution
    else:
        if isSmall(m.b):
            p, q, m, n = q, p, n, m

        t = (p - q)*2
        B1 = m.a*p.x + m.b*p.y + m.c
        B2 = n.a*q.x + n.b*q.y + n.c

        cA = (m.a*m.a + m.b*m.b)*(m.a*n.b - n.a*m.b)
        cB = (-(n.a*t.x + B2)*m.b*m.b*m.b
              + (B1*n.b + m.a*(n.a*t.y + n.b*t.x))*m.b*m.b
              - (2*B1*n.a + m.a*(t.y*n.b + B2))*m.a*m.b + 3*B1*m.a*m.a*n.b)
        cC = ((n.a*t.y + n.b*t.x)*m.b*m.b - (2*m.a*(n.b*t.y + B2) + B1*n.a)*m.b
              + 3*B1*m.a*n.b)*B1
        cD = -((n.b*t.y + B2)*m.b - B1*n.b)*B1*B1

        raizes = solveCubic(cA, cB, cC, cD)

        for raiz in raizes:
            xa = raiz + p.x
            ya = -(m.a*raiz + B1)/m.b + p.y
            solution += axiom_1(p, Point(xa, ya))

    return solution


def axiom_8(p, m, n):
    '''
    Given two lines m and n, and a point p not on m, fold to place P onto m,
    and to reflect half of n onto its other half.
    '''

    if p in m:
        raise ValueError('axiom_8: point p is in line m')

    if m.isParallel(n):
        return []

    d = n.vectProd(m)
    t = n.getNormal() @ p
    r = (m.getParallel()*t + n.getParallel()*m.c)/d
    return axiom_1(p, r)


if __name__ == '__main__':
    '''
    Compute cubic root of 2
    '''

    # Set points p, q, and lines m, n

    p = Point(0, -1)
    q = Point(-2, 0)
    m = Line(0, 1, -1)
    n = Line(1, 0, -2)

    # Place p onto m and q onto n

    fold = axiom_7(p, q, m, n)[0]  # Get the fold line
    a, b, c = fold.getCoefs()  # Get the coeficients

    # Compute the x intercept of the fold line

    x = -c/a
    print('The cubic root of 2 is {:.5f}\n'.format(x))
