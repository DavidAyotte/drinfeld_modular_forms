"""
Base class for Carlitz modules (rank 1 Drinfeld module)

AUTHORS:

- David Ayotte (2021): initial version
"""

# ****************************************************************************
#       Copyright (C) 2021 DAVID AYOTTE <davidayotte94@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

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
        r"""
        Return the string representation of ``self``.

        TESTS::

            sage: from drinfeld_modules.all import *
            sage: A.<T> = GF(3)['T']
            sage: C = CarlitzModule(A)
            sage: C
            Carlitz module over Finite Field of size 3 defined by:
            T |--> 𝜏 + T
        """
        return "Carlitz module over %s defined by:\n%s |--> %s" % (self.ring_of_constants(), self.base_ring().gen(), self.operator_polynomial())

    def inverse_cyclotomic_polynomial(self, a, name='X'):
        r"""
        Return the inverse cyclotomic polynomial in ``name`` of a polynomial `a`.

        EXAMPLES::

            sage: from drinfeld_modules.all import *
            sage: A.<T> = GF(3)['T']
            sage: C = CarlitzModule(A)
            sage: C.inverse_cyclotomic_polynomial(A.one())
            1
            sage: C.inverse_cyclotomic_polynomial(T)
            T*X^2 + 1
            sage: C.inverse_cyclotomic_polynomial(T^2)
            T^2*X^8 + (T^3 + T)*X^6 + 1
            sage: C.inverse_cyclotomic_polynomial(T^2 + T + 1)
            (T^2 + T + 1)*X^8 + (T^3 + T + 1)*X^6 + 1
        """
        act = self.action_polynomial(a, name=name)
        if a.is_zero():
            return act.parent().zero()
        X = act.parent().gen()
        N = self.q() ** a.degree()
        return act.subs(X ** (-1)) * X ** N

    def bracket(self, n):
        r"""
        Return the element `[n] = T^(q^n) - T` where `T` is the generator of the
        polynomial base ring of self.

        EXAMPLES::

            sage: from drinfeld_modules.all import *
            sage: A.<T> = GF(3)['T']
            sage: C = CarlitzModule(A)
            sage: C.bracket(1)
            T^3 + 2*T
            sage: C.bracket(2)
            T^9 + 2*T
            sage: C.bracket(0)
            Traceback (most recent call last):
            ...
            ValueError: the intger n (=0) must be postive.
        """
        if n <= 0:
            raise ValueError("the intger n (=%s) must be postive." % (n))
        T = self.base_field().gen()
        return T ** (self.q() ** n) - T

    def product_of_monic_polynomials(self, n):
        r"""
        Return the product of all monic polynomials in `\mathbb{F}_q[T]` of
        degree `n`.

        An alias of this method is ``D``.

        EXAMPLES::

            sage: from drinfeld_modules.all import *
            sage: A.<T> = GF(3)['T']
            sage: C = CarlitzModule(A)
            sage: C.product_of_monic_polynomials(0)
            1
            sage: C.product_of_monic_polynomials(1)
            T^3 + 2*T
            sage: f = C.product_of_monic_polynomials(2); f
            T^18 + 2*T^12 + 2*T^10 + T^4
            sage: f.factor()
            T^4 * (T + 1)^4 * (T + 2)^4 * (T^2 + 1) * (T^2 + T + 2) * (T^2 + 2*T + 2)
        """
        if n < 0:
            raise ValueError("the integer n (=%s) must be non-negative." & (n))
        ans = self.base_field().one()
        if n == 0:
            return ans
        for j in range(1, n+1):
            ans = ans * self.bracket(j) ** (self.q() ** (n - j))
        return ans

    D = product_of_monic_polynomials # alias

    def lcm_of_monic_polynomials(self, n):
        r"""
        Return the least common multiple of all monic polynomials in
        `\mathbb{F}_q[T]` of degree `n`.

        An alias of this method is ``L``.

        EXAMPLES::

            sage: from drinfeld_modules.all import *
            sage: A.<T> = GF(3)['T']
            sage: C = CarlitzModule(A)
            sage: C.lcm_of_monic_polynomials(1)
            T^3 + 2*T
            sage: C.lcm_of_monic_polynomials(2)
            T^12 + 2*T^10 + 2*T^4 + T^2
            sage: C.lcm_of_monic_polynomials(3)
            T^39 + 2*T^37 + 2*T^31 + T^29 + 2*T^13 + T^11 + T^5 + 2*T^3
        """
        if n < 0:
            raise ValueError("the integer n (=%s) must be non-negative." & (n))
        ans = self.base_field().one()
        if n == 0:
            return ans
        for j in range(1, n+1):
            ans = ans * self.bracket(j)
        return ans

    L = lcm_of_monic_polynomials # alias

    def goss_polynomial(self, n, name='X'):
        r"""
        Return the `n`-th Goss polynomial for the Carlitz module.

        EXAMPLES::

            sage: from drinfeld_modules.all import *
            sage: A.<T> = GF(3)['T']
            sage: C = CarlitzModule(A)
            sage: C.goss_polynomial(1)
            X
            sage: C.goss_polynomial(2)
            X^2
            sage: C.goss_polynomial(3)
            X^3
            sage: C.goss_polynomial(4)
            X^4 + (1/(T^3 + 2*T))*X^2
            sage: C.goss_polynomial(5)
            X^5 + (2/(T^3 + 2*T))*X^3
            sage: C.goss_polynomial(6)
            X^6
        """
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

    def ta(self, a, prec=10, name='t'):
        r"""
        Return the function `t(az)` as a power series in `t` truncated up to
        `prec`.

        EXAMPLES::

            sage: from drinfeld_modules.all import *
            sage: A.<T> = GF(3)['T']
            sage: C = CarlitzModule(A)
            sage: C.ta(A.one())
            t
            sage: C.ta(T)
            t^3 + 2*T*t^5 + T^2*t^7 + 2*T^3*t^9 + T^4*t^11 + O(t^13)
            sage: C.ta(T^2)
            t^9 + (2*T^3 + 2*T)*t^15 + 2*T^2*t^17 + O(t^19)
        """
        # TODO: add polynomial type verification
        # TODO: fix default precision issue
        if not a:
            return ValueError("the polynomial must be non-zero")
        R = PowerSeriesRing(self.base_polynomial_ring(), name=name, default_prec=prec) # understand default_prec
        d = a.degree()
        fpol = R(self.inverse_cyclotomic_polynomial(a, name=name))
        u = R.gen()
        return u ** (self.q() ** d) * (1/fpol)
