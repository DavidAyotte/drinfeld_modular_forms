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

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: (M.0)._repr_()
            'g0'
            sage: M.0 + M.1
            g0 + g1
        """
        return str(self.polynomial)

    def _add_(self, other):
        r"""
        Return the addition of self with other.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.0 + M.1  # indirect doctest
            g0 + g1
        """
        return self.__class__(self.parent(), self.polynomial + other.polynomial)

    def _mul_(self, other):
        r"""
        Return the multiplication of self with other.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.0*M.1  # indirect doctest
            g0*g1
            sage: M.0*(M.0 + M.1)
            g0^2 + g0*g1
            sage: (M.0 + M.1)*M.0
            g0^2 + g0*g1
        """
        return self.__class__(self.parent(), self.polynomial * other.polynomial)

    def _lmul_(self, c):
        r"""
        Return the scalar multiplication of self by `c`.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
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

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: -M.0  # indirect doctest
            -g0
        """
        return self.__class__(self.parent(), -self.polynomial)

    def __bool__(self):
        r"""
        Return True wether self is zero.

        TESTS::

            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: bool(M.0)
            False
        """
        return bool(self.polynomial)

    def _richcmp_(self, other, op):
        r"""
        Return the comparison of self with other.

        TESTS::

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

    def is_one(self):
        r"""
        Return ``True`` wether the given graded Drinfeld form is the
        multiplicative identity.

        EXAMPLES::

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
        Return ``True`` wether the given graded Drinfeld form is the additive
        identity.

        EXAMPLES::

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
        return self.polynomial.is_zero()
