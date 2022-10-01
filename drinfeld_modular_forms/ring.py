from sage.categories.graded_algebras import GradedAlgebras

from sage.structure.parent import Parent

from sage.rings.function_field.function_field import FunctionField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder

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
        self._poly_ring = PolynomialRing(base_ring, rank, names=names, order=TermOrder('wdeglex', degs))

        Parent.__init__(self, base=base_ring, category=GradedAlgebras(base_ring))

    def _repr_(self):
        r"""
        Return the string representation of self.
        """
        return "Ring of Drinfeld modular forms of rank %s over %s" % (self._rank, self._base_ring)

    def gen(self, n):
        r"""
        Return the `n`-th generator of the ring.

        EXAMPLES::

            sage: from drinfeld_modular_forms.all import *
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

            sage: from drinfeld_modular_forms.all import *
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

            sage: from drinfeld_modular_forms.all import *
            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
            sage: M = DrinfeldModularFormsRing(K, 5)
            sage: M.ngens()
            5
        """
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

            sage: from drinfeld_modular_forms.all import *
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

            sage: from drinfeld_modular_forms.all import *
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
