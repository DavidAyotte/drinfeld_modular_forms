from drinfeld_modules import DrinfeldModule

from sage.rings.all import FiniteField
from sage.rings.real_double import RDF
from sage.rings.power_series_ring import PowerSeriesRing

from sage.functions.log import log


class CarlitzModule(DrinfeldModule):
    def __init__(self, polynomial_base_ring, name='𝜏'):
        one = polynomial_base_ring.one()
        DrinfeldModule.__init__(self, one, name=name)

    def _repr_(self):
        return "Carlitz module over %s defined by:\n%s |--> %s" % (self.ring_of_constants(), self.base_ring().gen(), self.operator_polynomial())

    def inverse_cyclotomic_polynomial(self, a, name='X'):
        act = self.action_polynomial(a, name=name)
        if a.is_zero():
            return act.parent().zero()
        X = act.parent().gen()
        N = self.q() ** a.degree()
        return act.subs(X ** (-1)) * X ** N

    def Br(self, n):
        if n <= 0:
            raise ValueError("the intger `n` (=%s) must be postive." % (n))
        T = self.base_field().gen()
        return T ** (self.q() ** n) - T

    def D(self, n):
        if n < 0:
            raise ValueError("the integer `n` (=%s) must be non-negative." & (n))
        ans = self.base_field().one()
        if n == 0:
            return ans
        for j in range(1, n+1):
            ans = ans * self.Br(j) ** (self.q() ** (n - j))
        return ans

    def L(self, n):
        if n < 0:
            raise ValueError("the integer `n` (=%s) must be non-negative." & (n))
        ans = self.base_field().one()
        if n == 0:
            return ans
        for j in range(1, n+1):
            ans = ans * self.Br(j)
        return ans

    def goss_polynomial(self, n, name='X'):
        R = self.base_polynomial_ring()['X']
        X = R.gen()
        q = self.q()
        pol = R.zero()
        if n == 0:
            return pol
        if n <= q - 1:
            return X ** n
        if (n % q) == 0:
            return self.goss_polynomial(n/q) ** q
        for j in range(0, RDF(log(n, q)).floor() + 1):
            pol = pol + X * self.goss_polynomial(n - q ** j) * (1/self.D(j))
        return pol

    def ta(self, polynomial, prec, name='u'):
        # TODO: add polynomial type verification
        if not polynomial:
            return ValueError("the polynomial must be non-zero")
        R = PowerSeriesRing(self.base_polynomial_ring(), name=name, default_prec=prec)
        d = polynomial.degree()
        fpol = R(self.inverse_cyclotomic_polynomial(polynomial, name=name))
        u = R.gen()
        return u ** (self.q() ** d) * (1/fpol)
