from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp, op_NE, op_EQ

class DrinfeldModularFormsRingElement(ModuleElement):
    r"""
    Element class of rings of Drinfeld modular forms.
    """
    def __init__(self, parent, polynomial):
        self.polynomial = polynomial

        ModuleElement.__init__(self, parent)

    def _repr_(self):
        r"""
        Return the string representation of self.
        """
        return str(self.polynomial)

    def _add_(self, other):
        return self.__class__(self.parent(), self.polynomial + other.polynomial)

    def _mul_(self, other):
        return self.__class__(self.parent(), self.polynomial * other.polynomial)

    def _lmul_(self, c):
        return self.__class__(self.parent(), c * self.polynomial)

    def __neg__(self):
        return self.__class__(self.parent(), -self.polynomial)

    def __bool__(self):
        return bool(self.polynomial)

    def _richcmp_(self, other, op):
        if op != op_EQ and op != op_NE:
            raise TypeError('invalid comparison between modular forms ring elements')
        return richcmp(self.polynomial, other.polynomial, op)

    def is_one(self):
        return self.polynomial.is_one()

    def is_zero(self):
        return self.polynomial.is_zero()
