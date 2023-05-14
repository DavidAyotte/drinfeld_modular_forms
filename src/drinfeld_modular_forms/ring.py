r"""
This module defines a class named :class:`DrinfeldModularFormsRing`.

Currently, the implementation only supports the full group modular group
`\mathrm{GL}_r(A)` where `A = \mathbb{F}_q[T]`.

The implementation is based on the following identification:

.. MATH::

    M^r(\mathrm{GL}_r(A))
    = \mathbb{C}_{\infty}[g_1, \ldots, g_{r-1}, g_{r}].

where `g_i` the `i`-th coefficient form of weight `q^{i} - 1`.

EXAMPLES::

    sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
    sage: q = 3
    sage: A = GF(q)['T']
    sage: K.<T> = Frac(A)
    sage: M = DrinfeldModularFormsRing(K, 2)  # rank 2
    sage: M.gens()  # generators
    [g1, g2]
    sage: M.inject_variables()  # assign variables
    Defining g1, g2
    sage: g1.weight()
    2
    sage: g2.weight()
    8

::

    sage: M = DrinfeldModularFormsRing(K, 3)  # rank 3
    sage: M.gens()
    [g1, g2, g3]
    sage: [g.weight() for g in M.gens()]  # list of the weights
    [2, 8, 26]
    sage: M.inject_variables()
    Defining g1, g2, g3
    sage: g1.weight() == 3 - 1
    True
    sage: g2.weight() == 3^2 - 1
    True
    sage: g3.weight() == 3^3 - 1
    True

One can compute basis for any subspace of given weight::

    sage: M = DrinfeldModularFormsRing(K, 4)
    sage: M.basis_of_weight(q^4 - 1)
    [g4,
     g2^10,
     g1*g3^3,
     g1^2*g2^3*g3^2,
     g1^3*g2^6*g3,
     g1^4*g2^9,
     g1^6*g2^2*g3^2,
     g1^7*g2^5*g3,
     g1^8*g2^8,
     g1^10*g2*g3^2,
     g1^11*g2^4*g3,
     g1^12*g2^7,
     g1^14*g3^2,
     g1^15*g2^3*g3,
     g1^16*g2^6,
     g1^19*g2^2*g3,
     g1^20*g2^5,
     g1^23*g2*g3,
     g1^24*g2^4,
     g1^27*g3,
     g1^28*g2^3,
     g1^32*g2^2,
     g1^36*g2,
     g1^40]

We note that the elements of this ring may not be *modular forms* as
as they may have mixed weight components::

    sage: M = DrinfeldModularFormsRing(K, 4)
    sage: M.inject_variables()
    Defining g1, g2, g3, g4
    sage: F = g1 + g2 + g3 + g4
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
    sage: M.inject_variables()
    Defining g1, g2
    sage: g1.expansion()
    1 + ((2*T^3+T)*t^2) + O(t^7)
    sage: g2.expansion()
    t^2 + 2*t^6 + O(t^8)

The returned series is a lazy power series, meaning that it can compute
any coefficient at any precision on demands::

    sage: g2[6]
    2
    sage: g2[24]
    T^9 + 2*T^3
    sage: g2[36]
    2*T^9 + T^3
    sage: g2[702]  # long time
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
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.lazy_import import lazy_import

from .element import DrinfeldModularFormsRingElement
from .expansions import compute_petrov_expansion

lazy_import('sage.rings.lazy_series', 'LazyPowerSeries')
lazy_import('sage.rings.power_series_ring_element', 'PowerSeries')

class DrinfeldModularFormsRing(Parent, UniqueRepresentation):
    r"""
    Base class for the graded Drinfeld modular forms ring.

    INPUT:

    - ``base_ring`` -- The fraction field of a univariate polynomial
      ring over `\mathbb{F}_q`.
    - ``rank`` (integer, default: 2) -- the rank of the ring
    - ``group`` (NoneType) -- the group of self. The current
      implementation only supports the full group
      `\mathrm{GL}_r(A)`.
    - ``has_type`` (bool, default: ``False``) -- if set to True, returns
      the graded ring of arbitrary type. Currently only implemented in
      rank two.
    - ``names`` (string, default: ``'g'``) -- a single character or a
      comma seperated string of character representing the names of the
      generators.

    TESTS::

        sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
        sage: K = Frac(GF(3)['T'])
        sage: TestSuite(DrinfeldModularFormsRing(K)).run()
        sage: TestSuite(DrinfeldModularFormsRing(K, 3)).run()
        sage: TestSuite(DrinfeldModularFormsRing(K, 4)).run()
        sage: K = Frac(GF(7)['T'])
        sage: TestSuite(DrinfeldModularFormsRing(K)).run()
        sage: TestSuite(DrinfeldModularFormsRing(K, 3)).run()
        sage: TestSuite(DrinfeldModularFormsRing(K, 4)).run()
    """

    Element = DrinfeldModularFormsRingElement

    def __init__(self, base_ring, rank=2, group=None, has_type=False, names='g'):
        if group is not None:
            raise NotImplementedError("the ring of Drinfeld modular "
                                      "forms is only implemented for "
                                      "the full group")
        if not isinstance(names, str):
            raise TypeError("names must be string type")
        q = base_ring.base_ring().cardinality()
        if has_type:
            if rank != 2:
                raise NotImplementedError("ring with type are not "
                                          "implemented in rank =/= 2")
            if len(names) == 1:
                names += "1, h"
            degs = [q - 1, q + 1]
        else:
            if len(names) == 1:
                n = names
                names += "1, "
                for i in range(2, rank, 1):
                    names += n + str(i) + ", "
                names += n + str(rank)
            else:
                if len(names.split()) != rank:
                    raise ValueError("the rank does not corresponds to"
                                    " the number of generators")
            degs = [q ** i - 1 for i in range(1, rank + 1, 1)]
        self._has_type = has_type
        self._rank = rank
        self._base_ring = base_ring
        self._poly_ring = PolynomialRing(base_ring, rank, names=names,
                                         order=TermOrder('wdeglex', degs))
        self._assign_names(names)

        Parent.__init__(self, base=base_ring, category=GradedAlgebras(base_ring))

    def _an_element_(self):
        r"""
        Return an element of self.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: q = 3
            sage: A = GF(q)['T']
            sage: K.<T> = Frac(A)
            sage: M = DrinfeldModularFormsRing(K)
            sage: M.an_element()
            g1
        """
        return self.element_class(self, self._poly_ring.an_element())

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
            g1
            sage: f = (M.1).expansion()
            sage: M.from_expansion(f, (M.1).weight())
            g2
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

    def coefficient_form(self, i):
        r"""
        Return the `i`-th coefficient form of the universal Drinfeld
        module over `\Omega^r(\mathbb{C}_{\infty})`.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: q = 3
            sage: A = GF(q)['T']
            sage: K.<T> = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 3)
            sage: M.coefficient_form(1)
            g1
            sage: M.coefficient_form(2)
            g2
            sage: M.coefficient_form(3)
            g3

        ::

            sage: M = DrinfeldModularFormsRing(K, 2, has_type=True)
            sage: M.coefficient_form(1)
            g1
            sage: M.coefficient_form(2)
            h^2
        """
        if i not in ZZ:
            raise TypeError("i must be an integer")
        i = ZZ(i)
        if i < 1:
            raise ValueError("i must be >= 1")
        if i == 1:
            return self.gen(0)
        if self._has_type and i == 2:
            q = self._base_ring.base_ring().cardinality()
            return self.gen(1)**(q - 1)
        return self.gen(i - 1)

    def gen(self, n):
        r"""
        Return the `n`-th generator of the ring.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.0
            g1
            sage: M.1
            g2
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
            [g1, g2, g3, g4, g5]
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
            g1
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
            g2
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
            [g1]
            sage: M.basis_of_weight(q^2 - 1)
            [g2, g1^4]
            sage: M.basis_of_weight(q^3 - 1)
            [g1*g2^3, g1^5*g2^2, g1^9*g2, g1^13]
            sage: M.basis_of_weight(19*(q-1))
            [g1^3*g2^4, g1^7*g2^3, g1^11*g2^2, g1^15*g2, g1^19]
        """
        return [self(mon) for mon in self._poly_ring.monomials_of_degree(k)]

    basis = basis_of_weight  # alias

    def polynomial_ring(self):
        r"""
        Return the multivariate polynomial ring over the base ring where
        each variable corresponds to a generator of a ring.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: q = 3; A = GF(q)['T']; K = Frac(A);
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: P = M.polynomial_ring()
            sage: P
            Multivariate Polynomial Ring in g1, g2 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3

        The degree of the variables corresponds to the weight of the
        associated generator::

            sage: P.inject_variables()
            Defining g1, g2
            sage: g1.degree()
            2
            sage: g2.degree()
            8
        """
        return self._poly_ring

    def petrov_expansion(self, k, n, name='t'):
        r"""
        Return a Drinfeld modular form which admits the `A`-expansion
        for `k` and `n` as defined by Petrov.

        Recall that it is defined by:

        .. MATH::

            f_{k, i}(z) := \sum_{\substack{a\in \mathbb{F}_q[T] \\ a\text{ monic}}} a^{k - i}G_i(t(az))

        where `k` and `n` are 2 integers which are divisible by `q - 1` and
        `n \leq p^{v_{p}(k - n)}` (`p` is the characteristic).

        INPUT:

        - ``k`` -- an integer equal to the weight of the resulting form.
        - ``n`` -- an integer congruent to the type modulo `q-1` which is
          stricly less than `p^{v_{p}(k - n)}`.
        - ``name`` -- string (default: ``'T'``), the name of the
          parameter at infinity.

        OUTPUT: a Drinfeld modular form of weight `k`.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: q = 3
            sage: A = GF(q)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.petrov_expansion((q + 1)*(q - 1), q - 1)
            g2
            sage: M.petrov_expansion((q^2 + 1)*(q - 1), q - 1)
            g1^6*g2
            sage: M.petrov_expansion((q^3 + 1)*(q - 1), q - 1)
            g1^24*g2 + (T^27 - T^9)*g1^12*g2^4 + (T^54 + T^36 + T^18)*g2^7

        ::

            sage: M = DrinfeldModularFormsRing(K, 2, has_type=True)
            sage: M.petrov_expansion(q + 1, 1)
            h
            sage: M.petrov_expansion(q^2 - 1, 1)
            g1^2*h
            sage: M.petrov_expansion(q^3 - 1, 1)
            g1^11*h + (T^21 - T^19 + T^15 - T^13 + T^9 - T^7)*g1^7*h^3 + (-T^24 - T^22 - T^20 - T^18 - T^16 - T^14 - T^12 - T^10 - T^8)*g1^3*h^5
        """
        if self._rank != 2:
            raise NotImplementedError("A-expansions are only known in rank 2 for the moment")
        if k not in ZZ or n not in ZZ:
            raise TypeError("k and n must be integers")
        if k == n:
            raise ValueError("k must be different from n")
        q = self.base_ring().base_ring().cardinality()
        k = ZZ(k)
        n = ZZ(n)
        p = self.base_ring().characteristic()
        if k - 2*n <= 0:
            raise ValueError(f"k - 2n must be > 0 (current value: {k - 2*n})")
        if self._has_type:
            if (k - 2*n)%(q - 1):
                raise ValueError("k - 2n must be a multiple of q - 1")
        else:
            if k%(q-1):
                raise ValueError("k must be divisible by q - 1")
            if n%(q-1):
                raise ValueError("n must be divisible by q - 1")
        if n > p**((k-n).valuation(p)):
            raise ValueError(f"n (={n}) must be less than p^val(p, k - n)"
                             f"(={p**((k-n).valuation(p))})")
        expansion = compute_petrov_expansion(k, n, self.base_ring().base(), name=name, check=False)
        return self.from_expansion(expansion, k)

    def sturm_bound(self, k):
        r"""
        Return the Sturm bound of the subspace of weight `k`.

        Currently only implemented in rank 2.

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
        if self._rank != 2:
            raise NotImplementedError
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
        of the form `q^i - 1`. If `k = 0`, the method returns 0.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: q = 3
            sage: A = GF(q)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.eisenstein_series(0)
            0
            sage: M.eisenstein_series(q - 1)
            g1
            sage: M.eisenstein_series(q^2 - 1)
            g1^4
            sage: M.eisenstein_series(q^3 - 1)
            g1^13 + (-T^9 + T)*g1*g2^3
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
            raise NotImplementedError("only implemented when k is of the form q^i - 1")
        k = (k + 1).log(q)
        if k == 1:
            return self.coefficient_form(1)
        T = self._base_ring.gen()
        sqb = T**(q**(k - 1)) - T
        return -sqb*self.eisenstein_series(q**(k - 2) - 1)*self.coefficient_form(2)**(q**(k - 2)) + self.eisenstein_series(q**(k - 1) - 1)*self.coefficient_form(1)**(q**(k - 1))
