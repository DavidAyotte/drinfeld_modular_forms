from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp, op_NE, op_EQ

from .expansions import compute_delta, compute_first_coefficient
from .carlitz_module import CarlitzModule

class DrinfeldModularFormsRingElement(ModuleElement):
    r"""
    Element class of rings of Drinfeld modular forms.

    EXAMPLES::

        sage: from drinfeld_modular_forms.all import *
        sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
        sage: M = DrinfeldModularFormsRing(K, 2)
        sage: M.0
        g0
        sage: M.1
        g1
        sage: (T^2 + 1)*(M.0 + M.0 * M.1)
        (T^2 + 1)*g0*g1 + (T^2 + 1)*g0
        sage: (M.0).parent()
        Ring of Drinfeld modular forms of rank 2 over Fraction Field of Univariate Polynomial Ring in T over Finite Field of size 3
        sage: M.1 in M
        True
    """
    def __init__(self, parent, polynomial):
        self.polynomial = polynomial

        ModuleElement.__init__(self, parent)

    def _repr_(self):
        r"""
        Return the string representation of self.

        TESTS::

            sage: from drinfeld_modular_forms.all import *
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: (M.0)._repr_()
            'g0'
            sage: M.0 + M.1
            g1 + g0
        """
        return str(self.polynomial)

    def _add_(self, other):
        r"""
        Return the addition of self with other.

        TESTS::

            sage: from drinfeld_modular_forms.all import *
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.0 + M.1  # indirect doctest
            g1 + g0
        """
        return self.__class__(self.parent(), self.polynomial + other.polynomial)

    def _mul_(self, other):
        r"""
        Return the multiplication of self with other.

        TESTS::

            sage: from drinfeld_modular_forms.all import *
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: M.0*M.1  # indirect doctest
            g0*g1
            sage: M.0*(M.0 + M.1)
            g0*g1 + g0^2
            sage: (M.0 + M.1)*M.0
            g0*g1 + g0^2
        """
        return self.__class__(self.parent(), self.polynomial * other.polynomial)

    def _lmul_(self, c):
        r"""
        Return the scalar multiplication of self by `c`.

        TESTS::

            sage: from drinfeld_modular_forms.all import *
            sage: A = GF(3)['T']; K = Frac(A); T = K.gen()
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

            sage: from drinfeld_modular_forms.all import *
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: -M.0  # indirect doctest
            -g0
        """
        return self.__class__(self.parent(), -self.polynomial)

    def __bool__(self):
        r"""
        Return True wether self is nonzero.

        TESTS::

            sage: from drinfeld_modular_forms.all import *
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: bool(M.0)
            True
        """
        return bool(self.polynomial)

    def _richcmp_(self, other, op):
        r"""
        Return the comparison of self with other.

        TESTS::

            sage: from drinfeld_modular_forms.all import *
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

    def rank(self):
        r"""
        Return the rank of graded Drinfeld modular form.

        EXAMPLES::

            sage: from drinfeld_modular_forms.all import *
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M2 = DrinfeldModularFormsRing(K, 2)
            sage: (M2.0).rank()
            2
            sage: M5 = DrinfeldModularFormsRing(K, 5)
            sage: (M5.0 + M5.3).rank()
            5
        """
        return self.parent()._rank

    def is_one(self):
        r"""
        Return ``True`` wether the given graded Drinfeld form is the
        multiplicative identity.

        EXAMPLES::

            sage: from drinfeld_modular_forms.all import *
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

            sage: from drinfeld_modular_forms.all import *
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
        return not bool(self)

    def t_expansion(self, max_deg=6, prec=6, name='t'):
        r"""
        Return the `t`-expansion of the graded Drinfeld form.

        EXAMPLES::

            sage: from drinfeld_modular_forms.all import *
            sage: A = GF(3)['T']; K = Frac(A)
            sage: M = DrinfeldModularFormsRing(K, 2)
            sage: g0, g1 = M.gens()
            sage: g0.t_expansion()  # long time
            1/(T^3 + 2*T) + t^2 + O(t^12)
            sage: g1.t_expansion()  # long time
            t^2 + 2*t^6 + (T^3 + 2*T)*t^8 + O(t^12)
            sage: F = (g0 + g1)*g0
            sage: F.t_expansion()  # long time
            1/(T^6 + T^4 + T^2) + 2*t^4 + (2/(T^3 + 2*T))*t^6 + (T^3 + 2*T)*t^10 + O(t^12)
        """
        C = CarlitzModule(self.base_ring().base())
        degs = self.polynomial.degrees()
        poly_ring = self.parent()._poly_ring
        g0, g1 = poly_ring.gens()
        sub_dict = {}
        if degs[0]:
            sub_dict[g0] = compute_first_coefficient(C, max_deg, prec, name)
        if degs[1]:
            sub_dict[g1] = compute_delta(C, max_deg, prec, name)
        return self.polynomial.subs(sub_dict)
