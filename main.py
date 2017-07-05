import math


def degree(a):
    if a == 0:
        return 0
    return int(math.log(a, 2)) + 1


def galois_add(a, b):
    return a ^ b


def galois_multiply(a, b, p):
    bits = degree(p)
    temp_copy = a
    msb = 1 << (bits - 1)
    shifting_bit = 1
    r = 0

    for i in range(bits):
        if b & shifting_bit:
            r = galois_add(r, temp_copy)
        temp_copy = temp_copy << 1
        shifting_bit = shifting_bit << 1
        if temp_copy & msb:
            temp_copy = temp_copy ^ p
    return r


def galois_multiply_no_reduce(a, b, p):
    bits = degree(p)
    temp_copy = a
    shifting_bit = 1
    r = 0

    for i in range(bits):
        if b & shifting_bit:
            r = galois_add(r, temp_copy)
        temp_copy = temp_copy << 1
        shifting_bit = shifting_bit << 1
    return r


def galois_quotient(a, b):
    if b == 0:
        raise ZeroDivisionError

    if a == 0:
        return 0

    ret = 0
    degree_diff = degree(a) - degree(b)
    while degree_diff >= 0:
        ret = ret ^ 1 << degree_diff
        a = galois_add(a, b << degree_diff)
        degree_diff = degree(a) - degree(b)
    return ret


def galois_div(a, b, p):
    return galois_multiply(a, galois_multiplicative_inverse(b, p), p)


def galois_multiplicative_inverse(a, p):
    t = 0
    newt = 1
    r = p
    newr = a

    while newr:
        quotient = galois_quotient(r, newr)
        r, newr = newr, galois_add(r, galois_multiply_no_reduce(quotient, newr, p))
        t, newt = newt, galois_add(t, galois_multiply_no_reduce(quotient, newt, p))

    return t


class GaliosFieldElement(object):
    '''Encapsulating the values and functions for galois arithmetic so that I can use override operations.
       For the purpose of this class, I'm only allowing exponents of 2'''
    def __init__(self, value, p):
        self.value = value
        self.p = p

    def __add__(self, other):
        if type(other) == int:
            return GaliosFieldElement(galois_add(self.value, other), self.p)
        return GaliosFieldElement(galois_add(self.value, other.value), self.p)

    def __mul__(self, other):
        if type(other) == int:
            return GaliosFieldElement(galois_multiply(self.value, other, self.p), self.p)
        return GaliosFieldElement(galois_multiply(self.value, other.value, self.p), self.p)

    def __pow__(self, power):
        assert power == 2
        return GaliosFieldElement(galois_multiply(self.value, self.value, self.p), self.p)

    def __eq__(self, other):
        if type(other) == int:
            return self.value == other
        else:
            return self.value == other.value

    def __div__(self, other):
        if type(other) == int:
            return GaliosFieldElement(galois_div(self.value, other, self.p), self.p)
        return GaliosFieldElement(galois_div(self.value, other.value, self.p), self.p)

    def __rdiv__(self, other):
        assert other == 1
        return GaliosFieldElement(galois_multiplicative_inverse(self.value, self.p), self.p)

    def __str__(self):
        return str(self.value)


# Elliptic Curve operations E: y^2 + xy = x^3 + ax^2 + b
class ECCPoint(object):
    def __init__(self, curve, x, y):
        '''x: GaloisFieldElement
           y: GaloisFieldElement,
           curve: ECC'''
        self.x = x
        self.y = y
        self.curve = curve

    def __add__(self, other):
        assert type(other) == ECCPoint
        if self == 0:
            return other

        if other == 0:
            return self

        if self.x == other.x:
            if self.y + other.y == self.x:
                return self.curve.generate_zero_point()
            else:
                lamb = other.x + other.y/other.x
                x2 = lamb**2 + lamb + self.curve.a
                y2 = other.x**2 + (lamb+1) * x2

        else:
            lamb = (self.y + other.y)/(self.x + other.x)
            x2 = lamb**2 + lamb + self.x + other.x + self.curve.a
            y2 = (other.x + x2)*lamb + x2 + other.y

        return self.curve.generate_point_from_galois_elements(x2, y2)

    def __sub__(self, other):
        p = self.curve.generate_point_from_galois_elements(self.x, self.x+self.y)
        return other + p

    def copy(self):
        return ECCPoint(self.curve, self.x, self.y)

    def __mul__(self, n):
        c = n
        Q = self.curve.generate_zero_point()
        p0 = self.copy()

        while c > 0:
            if c % 2:
                u = 2 - c % 4
                c = c - u
                if u == 1:
                    Q = Q + p0
                if u == -1:
                    Q = Q - p0
            c = c / 2
            p0 = p0 + p0
        return Q

    def __rmul__(self, n):
        return self * n

    def __radd__(self, other):
        return self + other

    def __eq__(self, other):
        if other == 0:
            return self.x == 0 and self.y == 0

        return self.x == other.x and self.y == other.y


class ECC(object):
    def __init__(self, a, b, p):
        self.a = GaliosFieldElement(a, p)
        self.b = GaliosFieldElement(b, p)
        self.p = p

    def generate_point_from_ints(self, x, y):
        return ECCPoint(self, GaliosFieldElement(x, self.p), GaliosFieldElement(y, self.p))

    def generate_point_from_galois_elements(self, x, y):
        return ECCPoint(self, x, y)

    def generate_zero_point(self):
        return ECCPoint(self, GaliosFieldElement(0, self.p), GaliosFieldElement(0, self.p))

    def test_point(self, point):
        return point.y**2 + point.x * point.y == \
               point.x * point.x * point.x + self.a * point.x * point.x + self.b


if __name__ == '__main__':
    poly = 2**163 + 2**7 + 2**6 + 2**3 + 1
    a = 1
    b = 1
    Gx = 0x2fe13c0537bbc11acaa07d793de4e6d5e5c94eee8
    Gy = 0x289070fb05d38ff58321f2e800536d538ccdaa3d9
    ecc = ECC(a, b, poly)
    point1 = ecc.generate_point_from_ints(Gx, Gy)
    point = ecc.generate_point_from_ints(Gx, Gy)

    print ecc.test_point(point * 734892174932789417298478913274987348903721890479018237498032790841723908740192)
