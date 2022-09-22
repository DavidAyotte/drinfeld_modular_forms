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

        if not isinstance(base_ring, FunctionField):
            raise ValueError

        self._base_ring = base_ring
        self._poly_ring = PolynomialRing(base_ring, r, names=names)

        Parent.__init__(base_ring)
