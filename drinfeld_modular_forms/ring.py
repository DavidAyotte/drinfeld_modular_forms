from sage.categories.graded_algebras import GradedAlgebras

from sage.structure.parent import Parent

from sage.rings.function_field.function_field import FunctionField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.integer import Integer

from .element import DrinfeldModularFormsRingElement

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
            sage: B = M.basis_of_weight(19*(q-1))
            [g0^3*g1^4, g0^7*g1^3, g0^11*g1^2, g0^15*g1, g0^19]
        """
        return [self(mon) for mon in self._poly_ring.monomials_of_degree(k)]

    def weighted_eisenstein_serie(self, k):
        r"""
        Return the Drinfeld Eisenstein serie of weight `q^k - 1`.

        The method is currently only implemented for rank 2.

        EXAMPLES::

            sage: from drinfeld_modular_forms import DrinfeldModularFormsRing
            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.weighted_eisenstein_serie(0)
            0
            sage: M.weighted_eisenstein_serie(1)
            g0
            sage: M.weighted_eisenstein_serie(2)
            g0^4
            sage: M.weighted_eisenstein_serie(3)
            g0^13 + (-T^9 + T)*g0*g1^3
        """
        if not isinstance(k, (Integer, int)):
            raise TypeError("k should be an integer")
        if k < 0:
            raise ValueError("the integer k (=%s) should be nonnegative" % (k))
        if self._rank != 2:
            raise NotImplementedError
        q = self._base_ring.base_ring().cardinality()
        if k == 0:
            return self.zero()
        if k == 1:
            return self.gen(0)
        T = self._base_ring.gen()
        sqb = T ** (q ** (k - 1)) - T
        return -sqb * self.weighted_eisenstein_serie(k - 2) * self.gen(1) ** (q ** (k - 2)) + self.weighted_eisenstein_serie(k - 1) * self.gen(0) ** (q ** (k - 1))
