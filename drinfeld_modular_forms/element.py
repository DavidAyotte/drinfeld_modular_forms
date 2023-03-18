r"""
Elements of Drinfeld modular forms rings

The module defines the element of the class
:class:`~drinfeld_modular_forms.ring.DrinfeldModularFormsRing` from the
module :mod:`~drinfeld_modular_forms.ring`. See the aforementioned
module for definitions.

EXAMPLES::

    sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
    sage: q = 3
    sage: A = GF(q)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularFormsRing(K, 2)  # rank 2
    sage: g0, g1 = M.gens()
    sage: g0.parent()
    Ring of Drinfeld modular forms of rank 2 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3

AUTHORS:

- David Ayotte (2022): initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 DAVID AYOTTE <davidayotte94@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.lazy_import import lazy_import

from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp, op_NE, op_EQ

from .expansions import compute_delta_rank_2, compute_eisentein_serie_rank_2

lazy_import('sage.rings.lazy_series_ring', 'LazyPowerSeriesRing')

class DrinfeldModularFormsRingElement(ModuleElement):
    r"""
    Element class of rings of Drinfeld modular forms.

    EXAMPLES::

        sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
        sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
        sage: M = DrinfeldModularFormsRing(K, 2)
        sage: M.0
        g0
        sage: M.1
        g1
        sage: (T^2 + 1)*(M.0 + M.0 * M.1)
        (T^2 + 1)*g0*g1 + (T^2 + 1)*g0
        sage: (M.0).parent()
        Ring of Drinfeld modular forms of rank 2 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3
        sage: M.1 in M
        True
    """
    def __init__(self, parent, polynomial):
        # TODO: add checks
        self.polynomial = polynomial

        ModuleElement.__init__(self, parent)

    def _repr_(self):
        r"""
        Return the string representation of self.

        TESTS::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: (M.0)._repr_()
            'g0'
            sage: M.0 + M.1
            g1 + g0
        """
        return str(self.polynomial)

    def _add_(self, other):
        r"""
        Return the addition of self with other.

        TESTS::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.0 + M.1  # indirect doctest
            g1 + g0
        """
        return self.__class__(self.parent(), self.polynomial + other.polynomial)

    def _mul_(self, other):
        r"""
        Return the multiplication of self with other.

        TESTS::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.0*M.1  # indirect doctest
            g0*g1
            sage: M.0*(M.0 + M.1)
            g0*g1 + g0^2
            sage: (M.0 + M.1)*M.0
            g0*g1 + g0^2
        """
        return self.__class__(self.parent(), self.polynomial * other.polynomial)

    def _lmul_(self, c):
        r"""
        Return the scalar multiplication of self by `c`.

        TESTS::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: (T^2 + T + 2) * M.0  # indirect doctest
            (T^2 + T - 1)*g0
            sage: M.1 * (T^5 + T^2)
            (T^5 + T^2)*g1
            sage: 0 * M.1
            0
            sage: M.0 * 0
            0
        """
        return self.__class__(self.parent(), c * self.polynomial)

    def __neg__(self):
        r"""
        Return the negation of self.

        TESTS::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: -M.0  # indirect doctest
            -g0
        """
        return self.__class__(self.parent(), -self.polynomial)

    def __bool__(self):
        r"""
        Return True whether self is nonzero.

        TESTS::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: bool(M.0)
            True
        """
        return bool(self.polynomial)

    def _richcmp_(self, other, op):
        r"""
        Return the comparison of self with other.

        TESTS::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.0 == M.1
            False
            sage: M.0 != M.1
            True
            sage: M.0 == M.0
            True
        """
        if op != op_EQ and op != op_NE:
            raise TypeError('invalid comparison between modular forms ring elements')
        return richcmp(self.polynomial, other.polynomial, op)

    def __getitem__(self, n):
        r"""
        Return the `n`-coefficient of the Drinfeld modular form.

        This method is only implemented when the rank is 2.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: g0, g1 = M.gens()
            sage: g0[2]
            2*T^3 + T
            sage: g0[0:16]
            [1, 0, 2*T^3 + T, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*T^3 + T, 0]
            sage: g1[24]
            T^9 + 2*T^3
        """
        if self.parent()._rank != 2:
            raise NotImplementedError("expansion only implemented in rank 2")
        return self.expansion()[n]

    coefficient = __getitem__  # alias

    def rank(self):
        r"""
        Return the rank of graded Drinfeld modular form.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M2 = DrinfeldModularFormsRing(K, 2)
            sage: (M2.0).rank()
            2
            sage: M5 = DrinfeldModularFormsRing(K, 5)
            sage: (M5.0 + M5.3).rank()
            5
        """
        return self.parent()._rank

    def is_one(self):
        r"""
        Return ``True`` whether the given graded Drinfeld form is the
        multiplicative identity.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: u = M.one()
            sage: u.is_one()
            True
            sage: (M.0).is_one()
            False
        """
        return self.polynomial.is_one()

    def is_zero(self):
        r"""
        Return ``True`` whether the given graded Drinfeld form is the additive
        identity.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: z = M.zero()
            sage: z.is_zero()
            True
            sage: f = M.0
            sage: f.is_zero()
            False
            sage: (f - f).is_zero()
            True
            sage: (0 * M.0).is_zero()
            True
        """
        return not bool(self)

    def is_drinfeld_modular_form(self):
        r"""
        Return whether ``self`` is a Drinfeld modular form.

        We recall that elements of Drinfeld modular forms ring are not
        necessarily modular forms as they may have mixed weight components.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: g0, g1 = M.gens()
            sage: f = g0^5*g1^2  # homogeneous polynomial
            sage: f.is_drinfeld_modular_form()
            True
            sage: g = g0 + g1  # mixed weight components
            sage: g.is_drinfeld_modular_form()
            False
        """
        return self.polynomial.is_homogeneous()

    def expansion(self, name='t'):
        r"""
        Return the expansion at infinity of the graded Drinfeld form.

        OUTPUT: a lazy power series over the base ring.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: g0, g1 = M.gens()
            sage: g0.expansion()
            1 + ((2*T^3+T)*t^2) + O(t^7)
            sage: g1.expansion()
            t^2 + 2*t^6 + O(t^8)
            sage: F = (g0 + g1)*g0
            sage: F.expansion()
            1 + ((T^3+2*T+1)*t^2) + ((T^6+T^4+2*T^3+T^2+T)*t^4) + 2*t^6 + O(t^7)
        """
        A = self.base_ring().base()
        degs = self.polynomial.degrees()
        L = LazyPowerSeriesRing(self.base_ring(), name)
        poly_ring = self.parent()._poly_ring
        g0, g1 = poly_ring.gens()
        sub_dict = {}
        if not degs[0] and degs[1]:
            D = compute_delta_rank_2(A, name)
            E = D.parent().one()
        elif degs[0] and not degs[1]:
            E = compute_eisentein_serie_rank_2(A, name)
            D = E.parent().one()
        elif degs[0] and degs[1]:
            E = compute_eisentein_serie_rank_2(A, name)
            D = compute_delta_rank_2(A, name)
        else:
            return L(self.polynomial)
        t_exp = L.zero()
        for c, (n, m) in zip(self.polynomial.coefficients(), self.polynomial.exponents()):
            t_exp += c*(E**n)*(D**m)
        return t_exp

    def weight(self):
        r"""
        Return the weight of self.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: g0, g1 = M.gens()
            sage: g0.weight()
            2
            sage: g1.weight()
            8
            sage: f = g0^5*g1^2
            sage: f.weight()
            26

        If the form is not modular, then the method returns an error::

            sage: f = g0 + g1
            sage: f.weight()
            Traceback (most recent call last):
            ...
            ValueError: the given ring element is not a Drinfeld modular form
        """
        if not self.is_drinfeld_modular_form():
            raise ValueError("the given ring element is not a Drinfeld modular form")
        return self.polynomial.degree()
