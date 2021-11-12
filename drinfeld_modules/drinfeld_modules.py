"""
Base class for Drinfeld modules

AUTHORS:

- DAVID AYOTTE (2021): initial version
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

from sage.rings.all import Polynomial, ZZ
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.rings.polynomial.ore_polynomial_element import OrePolynomial
from sage.rings.fraction_field_element import is_FractionFieldElement
from sage.rings.fraction_field import is_FractionField


from sage.structure.parent import Parent

class DrinfeldModule(Parent):
    def __init__(self, *args, name='𝜏'):
        if not isinstance(name, str):
            raise TypeError('the name of the Frobenius must be string')
        if len(args) == 1:
            args = args[0]
            if isinstance(args, OrePolynomial):
                K = args.base_ring()
                if not is_FractionField(K):
                    raise ValueError("the base ring of the given Ore polynomial should be the fraction field of a univariate polynomial ring")
                if not K.characteristic():
                    raise ValueError("the characteristic should be non-zero")
                if not args.constant_coefficient() == K.gen():
                    raise ValueError("the constant coefficient should be %s" % (A.gen()))
                rank = args.degree()
                if rank == 0:
                    raise ValueError("the degree of the defining Ore polynomial should be > 0")
                operator_polynomial = args
                ore_polynomial_ring = operator_polynomial.parent()
            elif isinstance(args, list):
                if len(args) == 0:
                    raise ValueError("the defining list must be nonempty")
                rank = ZZ(len(args))
                for idx, c in enumerate(args):
                    if not is_FractionFieldElement(c) and not isinstance(c, Polynomial):
                        raise ValueError("the elements of the list must be a fraction field element or a univariate polynomial ring element")
                    if isinstance(c, Polynomial):
                        R = c.parent().fraction_field()
                        c = R.coerce(c)
                        args[idx] = c
                    if not args[0].parent() == c.parent():
                        raise ValueError("inconsistent parent type in the defining list")
                K = args[0].parent()
                if not K.characteristic():
                    raise ValueError("the characteristic should be non-zero")
                Frob = K.frobenius_endomorphism()
                ore_polynomial_ring = OrePolynomialRing(K, Frob, name)
                args.insert(0, K.gen())
                operator_polynomial = ore_polynomial_ring(args)
            elif is_FractionFieldElement(args) or isinstance(args, Polynomial):
                if isinstance(args, Polynomial):
                    R = args.parent().fraction_field()
                    args = R.coerce(args)
                rank = ZZ(1)
                K = args.parent()
                if not K.characteristic():
                    raise ValueError("the characteristic should be non-zero")
                Frob = K.frobenius_endomorphism()
                ore_polynomial_ring = OrePolynomialRing(K, Frob, name)
                operator_polynomial = ore_polynomial_ring([K.gen(), args])
            else:
                raise TypeError("the given argument should be an Ore polynomial, a single fraction field element or a list of fraction field elements")
        elif len(args) > 1:
            rank = len(args)
            coeff_list = []
            for c in args:
                if not is_FractionFieldElement(c) and not isinstance(c, Polynomial):
                    raise ValueError("every argument must elements of a fraction field or a polynomial")
                if isinstance(c, Polynomial):
                    R = c.parent().fraction_field()
                    c = R.coerce(c)
                coeff_list.append(c)
                if not coeff_list[0].parent() == c.parent():
                    raise ValueError("inconsistent parent type in the arguments")
            K = coeff_list[0].parent()
            if not K.characteristic():
                raise ValueError("the characteristic should be non-zero")
            Frob = K.frobenius_endomorphism()
            ore_polynomial_ring = OrePolynomialRing(K, Frob, name)
            coeff_list.insert(0, K.gen())
            operator_polynomial = ore_polynomial_ring(coeff_list)

        Parent.__init__(self, base=ore_polynomial_ring.base_ring())

        self._name = name
        self._rank = rank
        self._base_field = K
        self._operator_polynomial = operator_polynomial
        self._ore_polynomial_ring = ore_polynomial_ring

    def rank(self):
        r"""
        Return the rank of the given Drinfeld module.

        EXAMPLES::

            sage: K.<T> = FunctionField(GF(3^2))
            sage: DrinfeldModule(T).rank()
            1
            sage: DrinfeldModule(T, T^2, T).rank()
            3
        """
        return self._rank

    def operator_polynomial(self):
        r"""
        Return the Ore polynomial corresponding to the action of `T`

        EXAMPLES::

            sage: K.<T> = FunctionField(GF(3^2))
            sage: phi_T = DrinfeldModule(T).operator_polynomial(); phi_T
            T*Frob + T
            sage: phi_T.parent()
            Ore Polynomial Ring in Frob over Rational function field in T over Finite Field in z2 of size 3^2 twisted by Frob
        """
        return self._operator_polynomial

    def base_field(self):
        return self._base_field

    regular_functions_at_infinity = base_field # alias

    def constant_base_field(self):
        r"""
        Return the field of constant of the base function field.

        EXAMPLES::

            sage: K.<T> = FunctionField(GF(7))
            sage: DrinfeldModule(T).constant_base_field()
            Finite Field of size 7
        """
        return self._base_field.base_ring()

    ring_of_constants = constant_base_field # alias

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: K.<T> = FunctionField(GF(13**2))
            sage: DrinfeldModule(T, name='tau')
            Drinfeld Module over Rational function field in T over Finite Field in z2 of size 13^2 defined by:
            T |--> T*tau + T
        """
        return "Drinfeld Module of rank %s over %s defined by:\n%s |--> %s" % (self.rank(), self.ring_of_constants(), self.base_ring().gen(), self.operator_polynomial())

    def action_endomorphism(self, a):
        r"""
        Return the Ore polynomial corresponding to the action of `a`.

        INPUT:

        - ``a`` (Polynomial) -- A polynomial living in the ring of regular function of the base function field.
        """
        if not isinstance(a, Polynomial):
            raise TypeError("the input must be a univariate polynomial")
        if not a.base_ring() == self.constant_base_field():
            raise ValueError("the the constants field (%s) of the input is inconsistent with that of the Drinfeld module (%s)" % (a.base_ring(), self.constant_base_field()))
        f = self.operator_polynomial()
        action_endomorphism = self._ore_polynomial_ring.zero()
        for idx, c in enumerate(a.coefficients(sparse=False)):
            action_endomorphism += c * (f ** idx)
        return action_endomorphism

    def _element_constructor_(self, a):
        return self.action_endomorphism(a)
