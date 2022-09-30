from sage.categories.graded_algebras import GradedAlgebras

from sage.structure.parent import Parent

from sage.rings.function_field.function_field import FunctionField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from .element import DrinfeldModularFormsRingElement

class DrinfeldModularFormsRing(Parent):
    r"""
    Base class for the graded Drinfeld modular forms ring.
    """

    Element = DrinfeldModularFormsRingElement

    def __init__(self, base_ring, rank=2, group=1, names='g'):

        if rank != 2:
            raise NotImplementedError("only rank 2 is implemented")

        if group != 1:
            raise NotImplementedError("the ring of Drinfeld modular forms is only implemented for the full group")

        # if not isinstance(base_ring, PolynomialRing):
        #     raise ValueError

        self._rank = rank
        self._base_ring = base_ring
        self._poly_ring = PolynomialRing(base_ring, rank, names=names)

        Parent.__init__(self, base=base_ring, category=GradedAlgebras(base_ring))

    def _repr_(self):
        r"""
        Return the string representation of self.
        """
        return "Ring of Drinfeld modular forms of rank %s over %s" % (self._rank, self._base_ring)

    def gen(self, n):
        return self(self._poly_ring.gen(n))

    def gens(self):
        return [self(g) for g in self._poly_ring.gens()]

    def _element_constructor_(self, polynomial):
        r"""
        Return the element coresponding to the given polynomial.
        """
        pol_ring = polynomial.parent()
        if pol_ring != self._poly_ring:
            raise ValueError("cannot convert the given polynomial")
        return self.element_class(self, polynomial)

    def one(self):
        return self(self._poly_ring.one())

    def zero(self):
        return self(self._poly_ring.zero())
