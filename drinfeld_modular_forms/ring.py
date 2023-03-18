r"""
This module defines a class named :class:`DrinfeldModularFormsRing`.

Currently, the implementation only supports the full group modular group
`\mathrm{GL}_r(A)` where `A = \mathbb{F}_q[T]`.

The implementation is based on the following identification:

.. MATH::

    M^r(\mathrm{GL}_r(A))
    = \mathbb{C}_{\infty}[g_1, \ldots, g_{r-1}, g_{r}].

where `g_i` the `i`-th coefficient form of weight `q^{i} - 1`.

.. RUBRIC:: Examples of computations

In this package, the generators of the ring are the coefficients forms:

EXAMPLES::

    sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
    sage: q = 3
    sage: A = GF(q)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularFormsRing(K, 2)  # rank 2
    sage: M.gens()  # generators
    [g0, g1]
    sage: g0, g1 = M.gens()
    sage: g0.weight()
    2
    sage: g1.weight()
    8

::

    sage: M = DrinfeldModularFormsRing(K, 3)  # rank 3
    sage: M.gens()
    [g0, g1, g2]
    sage: [g.weight() for g in M.gens()]  # list of the weights
    [2, 8, 26]

One can compute basis for any subspace of given weight::

    sage: M = DrinfeldModularFormsRing(K, 4)
    sage: M.basis_of_weight(q^4 - 1)
    [g3,
     g1^10,
     g0*g2^3,
     g0^2*g1^3*g2^2,
     g0^3*g1^6*g2,
     g0^4*g1^9,
     g0^6*g1^2*g2^2,
     g0^7*g1^5*g2,
     g0^8*g1^8,
     g0^10*g1*g2^2,
     g0^11*g1^4*g2,
     g0^12*g1^7,
     g0^14*g2^2,
     g0^15*g1^3*g2,
     g0^16*g1^6,
     g0^19*g1^2*g2,
     g0^20*g1^5,
     g0^23*g1*g2,
     g0^24*g1^4,
     g0^27*g2,
     g0^28*g1^3,
     g0^32*g1^2,
     g0^36*g1,
     g0^40]

We note that the elements of this ring may not be *modular forms* as
as they may have mixed weight components::

    sage: M = DrinfeldModularFormsRing(K, 4)
    sage: g0, g1, g2, g3 = M.gens()
    sage: F = g0 + g1 + g2 + g3
    sage: F.is_drinfeld_modular_form()
    False

This is why we call these elements *graded Drinfeld modular forms*.

.. RUBRIC:: The rank 2 case

In rank 2 case, one can also compute the expansion at infinity of Drinfeld
modular forms.

EXAMPLES::

    sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
    sage: q = 3
    sage: A = GF(q)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularFormsRing(K, 2)  # rank 2
    sage: g0, g1 = M.gens()
    sage: g0.expansion()
    1 + ((2*T^3+T)*t^2) + O(t^7)
    sage: g1.expansion()
    t^2 + 2*t^6 + O(t^8)

The returned series is a lazy power series, meaning that it can compute
any coefficient at any precision on demands::

    sage: g1[6]
    2
    sage: g1[24]
    T^9 + 2*T^3
    sage: g1[36]
    2*T^9 + T^3
    sage: g1[702]  # long time
    2*T^252 + T^246 + 2*T^90 + T^84 + 2*T^36 + T^30 + T^18 + T^12 + T^6

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

from sage.categories.graded_algebras import GradedAlgebras

from sage.structure.parent import Parent

from sage.rings.function_field.function_field import FunctionField
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.structure.sequence import Sequence
from sage.misc.lazy_import import lazy_import

from .element import DrinfeldModularFormsRingElement
from .expansions import compute_petrov_expansion

lazy_import('sage.rings.lazy_series', 'LazyPowerSeries')
lazy_import('sage.rings.power_series_ring_element', 'PowerSeries')

class DrinfeldModularFormsRing(Parent):
    r"""
    Base class for the graded Drinfeld modular forms ring.
    """

    Element = DrinfeldModularFormsRingElement

    def __init__(self, base_ring, rank=2, group=1, names='g'):

        if group != 1:
            raise NotImplementedError("the ring of Drinfeld modular forms is only implemented for the full group")

        # if not isinstance(base_ring, PolynomialRing):
        #     raise ValueError

        self._rank = rank
        self._base_ring = base_ring
        q = base_ring.base_ring().cardinality()
        degs = [q ** i - 1 for i in range(1, rank + 1, 1)]
        self._poly_ring = PolynomialRing(base_ring, rank, names=names,
                                         order=TermOrder('wdeglex', degs))

        Parent.__init__(self, base=base_ring, category=GradedAlgebras(base_ring))

    def _repr_(self):
        r"""
        Return the string representation of self.
        """
        return ("Ring of Drinfeld modular forms of rank %s over %s"
                % (self._rank, self._base_ring))

    def from_expansion(self, expansion, k):
        r"""
        Return the Drinfeld modular form which corresponds to the given
        expansion.

        INPUT:

        - ``expansion`` -- a lazy power series, a power series, a list or a
          tuple. The precision or the length must be at least the Sturm bound
          of the weight `k` subspace.
        - ``k`` -- an integer representing the weight of the expected Drinfeld
          modular form.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: q = 3
            sage: A = GF(q)['T']
            sage: K.<T> = Frac(A)
            sage: M = DrinfeldModularFormsRing(K)
            sage: M.sturm_bound(q - 1)
            1
            sage: M.from_expansion([K.one()], q - 1)
            g0
            sage: f = (M.1).expansion()
            sage: M.from_expansion(f, (M.1).weight())
            g1
        """
        if self._rank != 2:
            raise NotImplementedError
        basis = self.basis_of_weight(k)
        bound = self.sturm_bound(k)
        if isinstance(expansion, (list, tuple)):
            coefficients = Sequence(expansion)[0:bound]
            if coefficients.universe() != self.base_ring():
                raise ValueError("incorrect parent for given list")
        elif isinstance(expansion, PowerSeries):
            if expansion.parent().base_ring() != self.base_ring():
                raise ValueError("incorrect parent for given power series")
            coefficients = expansion.coefficients()[0:bound]
        elif isinstance(expansion, LazyPowerSeries):
            if expansion.parent().base_ring() != self.base_ring():
                raise ValueError("incorrect parent for given lazy power series")
            coefficients = expansion[0:bound]
        else:
            raise TypeError("expansion must be a lazy power series, a power "
                            "series, a list or a tuple")
        if len(coefficients) < bound:
            raise ValueError("not enough coefficients, please provide"
                                f"at least {bound} coefficients")
        v = Matrix([coefficients])
        coeff_basis = [b.expansion()[0:bound] for b in basis]
        try:
            c = Matrix(coeff_basis).solve_left(v)
        except ValueError:
            raise ValueError("the given expansion does not correspond to a"
                             "form of the given weight")
        return sum(c[0][i]*b for i, b in enumerate(basis))

    def gen(self, n):
        r"""
        Return the `n`-th generator of the ring.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.0
            g0
            sage: M.1
            g1
        """
        return self(self._poly_ring.gen(n))

    def gens(self):
        r"""
        Return a list of generators for this ring.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 5)
            sage: M.gens()
            [g0, g1, g2, g3, g4]
        """
        return [self(g) for g in self._poly_ring.gens()]

    def ngens(self):
        r"""
        Return the number of generators of the ring.

        Note that the number of generators is equal to the rank.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 5)
            sage: M.ngens()
            5
        """
        return self._rank

    def _element_constructor_(self, polynomial):
        r"""
        Return the element corresponding to the given polynomial.
        """
        pol_ring = polynomial.parent()
        if pol_ring != self._poly_ring:
            raise ValueError("cannot convert the given polynomial")
        return self.element_class(self, polynomial)

    def rank(self):
        r"""
        Return the rank of the ring of Drinfeld modular forms.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A);
            sage: DrinfeldModularFormsRing(K, 2).rank()
            2
            sage: DrinfeldModularFormsRing(K, 3).rank()
            3
            sage: DrinfeldModularFormsRing(K, 4).rank()
            4
        """
        return self._rank

    def one(self):
        r"""
        Return the multiplicative identity of the ring.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.one()
            1
            sage: M.one() * M.0
            g0
            sage: M.one().is_one()
            True
        """
        return self(self._poly_ring.one())

    def zero(self):
        r"""
        Return the additive identity of the ring.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.zero()
            0
            sage: M.zero() + M.1
            g1
            sage: M.zero() * M.1
            0
            sage: M.zero().is_zero()
            True
        """
        return self(self._poly_ring.zero())

    def basis_of_weight(self, k):
        r"""
        Return a list of Drinfeld modular forms which forms a basis for the
        subspace of weight `k`.

        Note that if `k\not\equiv 0` modulo `q-1`, then the subspace is 0.

        An alias of this method is ``basis``.

        INPUT:

        - ``k`` -- an integer.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: q = 3; A = GF(q)['T']; K = Frac(A);
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.basis_of_weight(q - 1)
            [g0]
            sage: M.basis_of_weight(q^2 - 1)
            [g1, g0^4]
            sage: M.basis_of_weight(q^3 - 1)
            [g0*g1^3, g0^5*g1^2, g0^9*g1, g0^13]
            sage: M.basis_of_weight(19*(q-1))
            [g0^3*g1^4, g0^7*g1^3, g0^11*g1^2, g0^15*g1, g0^19]
        """
        return [self(mon) for mon in self._poly_ring.monomials_of_degree(k)]

    basis = basis_of_weight  # alias

    def petrov_expansion(self, k, n):
        r"""
        Return a Drinfeld modular form which admits the `A`-expansion
        for `k` and `n` as defined by Petrov.

        Recall that it is defined by:

        .. MATH::

            f_{k, i}(z) := \sum_{\substack{a\in \mathbb{F}_q[T] \\ a\text{ monic}}} a^{k - i}G_i(t(az))

        where `k` and `n` are 2 integers which are divisible by `q - 1` and
        `n \leq p^{v_{p}(k - n)}` (`p` is the characteristic).

        INPUT:

        - ``k`` -- an integer divisible by `q-1`.
        - ``n`` -- an integer divisible by `q-1` which is stricly less than
          `p^{v_{p}(k - n)}`.
        - ``name`` -- string, the name of the parameter at infinity.

        OUTPUT: a Drinfeld modular form of weight `k`.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: q = 3
            sage: A = GF(q)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.petrov_expansion((q + 1)*(q - 1), q - 1)
            g1
            sage: M.petrov_expansion((q^2 + 1)*(q - 1), q - 1)
            g0^6*g1
            sage: M.petrov_expansion((q^3 + 1)*(q - 1), q - 1)
            g0^24*g1 + (T^27 - T^9)*g0^12*g1^4 + (T^54 + T^36 + T^18)*g1^7
        """
        if self._rank != 2:
            raise NotImplementedError("A-expansions are only known in rank 2 for the moment")
        if k not in ZZ or n not in ZZ:
            raise TypeError("k and n must be integers")
        if k == n:
            raise ValueError("k must be different from n")
        q = self.base_ring().base_ring().cardinality()
        if k%(q-1):
            raise ValueError("k must be divisible by q - 1")
        if n%(q-1):
            raise ValueError("n must be divisible by q - 1")
        k = ZZ(k)
        n = ZZ(n)
        p = self.base_ring().characteristic()
        if n > p**(k-n).valuation(p):
            raise ValueError("n must be less or equal to p^v_p(k-n)")
        expansion = compute_petrov_expansion(k, n, self.base_ring().base())
        return self.from_expansion(expansion, k)

    def sturm_bound(self, k):
        r"""
        Return the Sturm bound of the subspace of weight `k`.

        INPUT:

        - ``k`` -- an integer

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: q = 3; A = GF(q)['T']; K = Frac(A);
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.sturm_bound(q - 1)
            1
            sage: M.sturm_bound(q^2 - 1)
            3
            sage: M.sturm_bound(q^9 - 1)
            4921
        """
        if k not in ZZ:
            raise TypeError("input must be an integer")
        k = ZZ(k)
        q = self._base_ring.base_ring().cardinality()
        if k%(q-1):
            return ZZ.zero()
        return ZZ((k/(q + 1)).floor() + 1)

    def eisenstein_series(self, k):
        r"""
        Return the Drinfeld Eisenstein series of weight `k`.

        The method is currently only implemented for rank 2 and when `k` is of
        of the form `q^v - 1`. If `k = 0`, the method returns 0.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: q = 3
            sage: A = GF(q)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.eisenstein_series(0)
            0
            sage: M.eisenstein_series(q - 1)
            g0
            sage: M.eisenstein_series(q^2 - 1)
            g0^4
            sage: M.eisenstein_series(q^3 - 1)
            g0^13 + (-T^9 + T)*g0*g1^3
        """
        if self._rank != 2:
            raise NotImplementedError
        if k not in ZZ:
            raise TypeError("k must be an integer")
        k = ZZ(k)
        if k < 0:
            raise ValueError("the integer k (=%s) should be nonnegative" % (k))
        q = self._base_ring.base_ring().cardinality()
        if k == 0:
            return self.zero()
        if not (k + 1).is_power_of(q):
            raise NotImplementedError("only implemented when k is of the form q^v - 1")
        k = (k + 1).log(q)
        if k == 1:
            return self.gen(0)
        T = self._base_ring.gen()
        sqb = T**(q**(k - 1)) - T
        return -sqb*self.eisenstein_series(q**(k - 2) - 1)*self.gen(1)**(q**(k - 2)) + self.eisenstein_series(q**(k - 1) - 1)*self.gen(0)**(q**(k - 1))
