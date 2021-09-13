"""
The ring of additive polynomial over a finite field

AUTHORS:

- DAVID AYOTTE (2021): initial version

TODO: add documentation and doctests!
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

from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.all import PolynomialRing
from sage.rings.function_field.function_field import FunctionField

from sage.structure.parent import Parent
from sage.structure.element import RingElement
from sage.structure.richcmp import op_EQ, op_NE

from .defaults import DEFAULT_VARIABLE, DEFAULT_FROBENIUS

class AdditivePolynomialElement(RingElement):
    r"""
    Elements of additive polynomials ring.

    An additive polynomial is an endomorphism over a rational function field `K`
    given by a polynomial in `\tau : x \mapsto x^q`, the Frobenius endomorphism,

    .. MATH::

        a_0 \tau^0 + a_1 \tau + \cdots + a_r \tau^r.

    Note that multiplication is given by composition:

    .. MATH::

        \tau \alpha = \alpha^q \tau

    for any `\alpha \in K`.
    """
    def __init__(self, parent, polynomial):
        if not isinstance(polynomial, Polynomial):
            raise ValueError("the given argument should be a univariate polynomial")
        if polynomial.base_ring() != parent.base_ring():
            raise ValueError("the base ring of the given polynomial is not consistent with the base ring of the parent")
        self._polynomial = polynomial
        RingElement.__init__(self, parent)

    def _repr_(self):
        return self._polynomial._repr_()

    def _add_(self, other):
        return self.__class__(self.parent(), self._polynomial + other._polynomial)

    def __neg__(self):
        return self.__class__(self.parent(), -self._polynomial)

    def _mul_(self, other):
        P = self._polynomial.parent()
        t = P.gen()
        q = self.base_ring().base_ring().cardinality()
        pol_list = []
        for i in range(self._polynomial.degree() + 1):
            c = self._polynomial[i]
            pol = P.zero()
            for j in range(other._polynomial.degree() + 1):
                k = other._polynomial[j]
                pol += c * k ** (q ** i) * t ** (j + i)
            pol_list.append(pol)
        mul_other = sum(p for p in pol_list)
        return self.__class__(self.parent(), mul_other)

    def _lmul_(self, c):
        return self.__class__(self.parent(), c * self._polynomial)

    def _rmul_(self, c):
        q = self.base_ring().base_ring().cardinality()
        rmul_pol = self._polynomial.zero()
        t = self._polynomial.parent().gen()
        for i in range(self._polynomial.degree() + 1):
            coef = self._polynomial[i]
            rmul_pol += c ** (q ** i) * coef * t ** i
        return self.__class__(self.parent(), rmul_pol)

    def _richcmp_(self, other, op):
        if op != op_EQ and op != op_NE:
            raise ValueError("invalid comparison between additive polynomial ring element")
        return self._polynomial == other._polynomial

    def __call__(a):
        q = self.base_ring().cardinality()
        return self._polynomial.subs(a ** q)




class AdditivePolynomials(Parent):
    r"""
    Base class for the endomorphism ring of additive polynomials
    """

    Element = AdditivePolynomialElement

    def __init__(self, ring_of_constants, frobenius_name=DEFAULT_FROBENIUS, name=DEFAULT_VARIABLE):

        if not isinstance(frobenius_name, str) or not isinstance(name, str):
            raise TypeError("the names must be a string")
        if len(frobenius_name.split()) != 1 or len(name.split()) != 1:
            raise ValueError("the names must be a single word without spaces")
        self._frobenius_name = frobenius_name
        self._name = name
        #ring_of_constants = polynomial_ring.base_ring()
        if not ring_of_constants.is_finite() or not ring_of_constants.is_field():
            raise ValueError("the ring of constants must be a finite field")
        self._ring_of_constants = ring_of_constants
        base_ring = PolynomialRing(self._ring_of_constants, self._name)

        self._polynomial_ring = PolynomialRing(base_ring, self._frobenius_name)
        Parent.__init__(self, base=base_ring)

    def ring_of_constants(self):
        r"""
        Return the ring of constants of the the given additive polynomial ring.

        EXAMPLES::

            sage: from drinfeld_modules import *
            sage: AdditivePolynomials(GF(3^6)).ring_of_constants()
            Finite Field in z6 of size 3^6
            sage: AdditivePolynomials(GF(5)).ring_of_constants()
            Finite Field of size 5
            sage: AdditivePolynomials(GF(7)).ring_of_constants()
            Finite Field of size 7
        """
        return self._ring_of_constants

    def _repr_(self):
        r"""
        Return the string representation of self.

        TESTS::

            sage: from drinfeld_modules import *
            sage: AdditivePolynomials(GF(5))
            Endomorphism ring of Additive Polynomial over Finite Field of size 5 generated by the Frobenius tau
        """
        return "Endomorphism ring of Additive Polynomial over %s generated by the Frobenius %s" % (self._ring_of_constants, self._frobenius_name)

    def _element_constructor_(self, x):
        r"""
        Construct an element living in self.

        TESTS::

            sage: from drinfeld_modules import *
            sage: Ktau = AdditivePolynomials(GF(5))
            sage: Ktau(0)
            0
            sage: Ktau(1)
            1
            sage: A.<T> = PolynomialRing(GF(5))
            sage: B.<tau> = PolynomialRing(A)
            sage: Ktau(T + tau)
            tau + T
        """
        if isinstance(x, Polynomial):
            elt = self.element_class(self, x)
        else:
            elt = self.element_class(self, self._polynomial_ring.coerce(x))
        return elt

    def _coerce_map_from_(self, other_ring):
        r"""
        Code to make the coercion framework work.

        TESTS::

            sage: from drinfeld_modules import *
            sage: Ktau = AdditivePolynomials(GF(5))
            sage: A.<T> = PolynomialRing(GF(5))
            sage: Ktau.has_coerce_map_from(A)
            True
            sage: Ktau.has_coerce_map_from(GF(5))
            True
            sage: Ktau.has_coerce_map_from(ZZ)
            True
            sage: Ktau.has_coerce_map_from(QQ)
            False
        """
        return self._polynomial_ring.has_coerce_map_from(other_ring)

    def gen(self):
        r"""
        Return the generator of the ring, that is the frobenius endomorphism
        `\tau:x\mapsto x^q`.

        EXAMPLES::

            sage: from drinfeld_modules import *
            sage: Ktau = AdditivePolynomials(GF(7))
            sage: tau = Ktau.gen(); tau
            tau
        """
        return self(self._polynomial_ring.gen())

    frobenius = gen # alias
    frobenius_endomorphism = gen # alias

    def polygen(self):
        r"""
        Return the generator of the underlying polynomial ring as an element of
        the additive polynomials ring.

        EXAMPLES::

            sage: from drinfeld_modules import *
            sage: Ktau = AdditivePolynomials(GF(13^3))
            sage: T = Ktau.polygen(); T
            T
            sage: T.parent()
            Endomorphism ring of Additive Polynomial over Finite Field in z3 of size 13^3 generated by the Frobenius tau
        """
        T = self.base_ring().gen()
        return self(self._polynomial_ring(T))

    def one(self):
        r"""
        Return the multiplicative identity of this ring.

        EXAMPLES::

            sage: from drinfeld_modules import *
            sage: Ktau = AdditivePolynomials(GF(19^2))
            sage: u = Ktau.one(); u
            1
            sage: u * u
            1
            sage: tau = Ktau.gen(); T = Ktau.polygen()
            sage: u * tau
            tau
            sage: (T + tau) * u
            tau + T
        """
        return self(self._polynomial_ring.one())

    def zero(self):
        r"""
        Return the zero element of this ring.

        EXAMPLES::

            sage: from drinfeld_modules import *
            sage: Ktau = AdditivePolynomials(GF(11))
            sage: z = Ktau.zero(); z
            0
            sage: z + z
            0
            sage: tau = Ktau.gen(); T = Ktau.polygen()
            sage: z + tau
            tau
            sage: T + z + tau
            tau + T
        """
        return self(self._polynomial_ring.zero())
