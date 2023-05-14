r"""
In analogy with the classical theory, any Drinfeld modular form
`f:\Omega^2(\mathbb{C}_{\infty}) \rightarrow \mathbb{C}_{\infty}` admits
an *expansion at infinity* of the form:

.. MATH::

    f(w) = \sum_{i = 0}^{\infty} a_i(f) t(w)^i

where `t(w) := 1/e(w)` and `e(w)` is the exponential of the Carlitz
module `T\mapsto T + \tau`.

We say that a Drinfeld modular forms of weight `k` admits a
*Petrov expansion* or an `A`-*expansion* if there exists an integer `n`
and coefficients `c_{a}(f)` such that

.. MATH::

    f =
    \sum_{\substack{a\in \mathbb{F}_q[T] \\ a\text{ monic}}}
    c_a(f)G_n(t(az))

where `G_n(X)` is the `n`-th Goss polynomial of the Carlitz module.

Petrov showed that there exists an infinite family of Drinfeld modular
forms with `A`-expansion:

.. MATH::

    f_{k, n} =
    \sum_{\substack{a\in \mathbb{F}_q[T] \\ a\text{ monic}}}
    a^{k - n}G_n(t(az))

provided that `k` and `n` are integers such that `k - 2n \equiv 0`
modulo `q - 1` and `n \leq p^{v_p(k - n)}`.

This module defines functions that compute the expansion at infinity of
the forms `f_{k, n}`.

EXAMPLES::

    sage: from drinfeld_modular_forms import compute_petrov_expansion
    sage: q = 3
    sage: D = compute_petrov_expansion(q + 1, 1, q); D
    t + t^5 + ((2*T^3+T)*t^7) + O(t^8)

It is possible to compute on demands any coefficients of the above
series::

    sage: D[59]  # 59-th coefficient
    2*T^27 + T^9 + T^3 + 2*T
    sage: D[0:10]  # first 10 coefficients
    [0, 1, 0, 0, 0, 1, 0, 2*T^3 + T, 0, 1]

AUTHORS:

- David Ayotte (2021): initial version
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

from drinfeld_modular_forms.drinfeld_modules import DrinfeldModule
from drinfeld_modular_forms.goss_polynomials import bracket, _carlitz_module_goss_polynomial

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.finite_field_constructor import GF

from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ

lazy_import('sage.rings.lazy_series_ring', 'LazyPowerSeriesRing')

def _format_param(param):
    r"""
    Return a tuple (q, A) where q is the cardinality of the field of
    constants and A is the functions ring.

    This utility function also checks the validity of the given
    parameter.

    INPUT:

    - ``param`` -- a prime power, a univariate polynomial ring over
      `\mathbb{F}_q` or the fraction field of such polynomial ring.

    .. doctest::
       :hide:

        sage: from drinfeld_modular_forms.expansions import _format_param
        sage: _format_param(5)
        (5, Univariate Polynomial Ring in T over Finite Field of size 5)
        sage: A = GF(11^2)['T']
        sage: _format_param(A)
        (121, Univariate Polynomial Ring in T over Finite Field in z2 of size 11^2)
        sage: _format_param(Frac(A))
        (121, Univariate Polynomial Ring in T over Finite Field in z2 of size 11^2)

    ::

        sage: _format_param(6)
        Traceback (most recent call last):
        ...
        ValueError: param must be a prime power
        sage: _format_param(Frac(QQ['T']))
        Traceback (most recent call last):
        ...
        ValueError: field of constants must be finite
        sage: _format_param(Frac(GF(3)['X, Y']))
        Traceback (most recent call last):
        ...
        TypeError: base of the fraction field must be a univariate polynomial ring
        sage: _format_param('x')
        Traceback (most recent call last):
        ...
        TypeError: param must be a prime power, a univariate polynomial ring over Fq or a fraction field of such ring
    """
    from sage.rings.fraction_field import FractionField_generic
    from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
    if param in ZZ:
        q = ZZ(param)
        if not q.is_prime_power():
            raise ValueError("param must be a prime power")
        polynomial_ring = PolynomialRing(GF(q), name='T')
    elif isinstance(param, FractionField_generic):
        polynomial_ring = param.base()
        if not isinstance(polynomial_ring, PolynomialRing_general):
            raise TypeError("base of the fraction field must be a univariate polynomial ring")
        if not polynomial_ring.base_ring().is_field() or not polynomial_ring.base_ring().is_finite():
            raise ValueError("field of constants must be finite")
        q = polynomial_ring.base_ring().cardinality()
    elif isinstance(param, PolynomialRing_general):
        polynomial_ring = param
        if not polynomial_ring.base_ring().is_field() or not polynomial_ring.base_ring().is_finite():
            raise ValueError("field of constants must be finite")
        q = polynomial_ring.base_ring().cardinality()
    else:
        raise TypeError("param must be a prime power, a univariate polynomial "
                        "ring over Fq or a fraction field of such ring")
    return (q, polynomial_ring)


def inverse_cyclotomic_polynomial(a, name='X'):
    r"""
    Return the polynomial `f_a(X) = \rho_a(X^{-1}) X^{q^{\mathrm{deg}(a)}}`
    where `\rho_a` is the Carlitz module.

    This function is used in :func:`.parameter_at_infinity`

    INPUT:

    - ``a`` (polynomial) -- univariate polynomial over a finite field.
    - ``name`` (string, default: ``'X'``) -- the name of the variable.

    EXAMPLES::

        sage: from drinfeld_modular_forms.expansions import inverse_cyclotomic_polynomial
        sage: A.<T> = GF(3)['T']
        sage: inverse_cyclotomic_polynomial(A.one())
        1
        sage: inverse_cyclotomic_polynomial(T)
        T*X^2 + 1
        sage: inverse_cyclotomic_polynomial(T^2)
        T^2*X^8 + (T^3 + T)*X^6 + 1
        sage: inverse_cyclotomic_polynomial(T^2 + T + 1)
        (T^2 + T + 1)*X^8 + (T^3 + T + 1)*X^6 + 1
    """
    A = a.parent()
    q = A.base_ring().cardinality()
    D = DrinfeldModule(A.one())
    act = D.action_polynomial(a, name=name)
    if a.is_zero():
        return act.parent().zero()
    X = act.parent().gen()
    N = q ** a.degree()
    return act.subs(X ** (-1)) * X ** N

def parameter_at_infinity(a, name='t'):
    r"""
    Return the function `t(aw)` as a power series in `t`.

    INPUT:

    - ``a`` (polynomial) -- univariate polynomial over a finite field.
    - ``name`` (string, default: ``'t'``) -- the name of the lazy power
      series ring generator.

    EXAMPLES::

        sage: from drinfeld_modular_forms import parameter_at_infinity
        sage: A.<T> = GF(3)['T']
        sage: parameter_at_infinity(A.one())
        t
        sage: parameter_at_infinity(T)
        t^3 + 2*T*t^5 + T^2*t^7 + 2*T^3*t^9 + O(t^10)
        sage: parameter_at_infinity(T^2)
        t^9 + ((2*T^3+2*T)*t^15) + O(t^16)
    """
    if a.is_zero():
        return ValueError("the polynomial must be non-zero")
    A = a.parent()
    q = A.base_ring().cardinality()
    R = LazyPowerSeriesRing(A.fraction_field(), name)
    d = a.degree()
    fpol = R(inverse_cyclotomic_polynomial(a, name=name))
    u = R.gen()
    return u**(q**d)*fpol.inverse_of_unit()

def coefficient_petrov_expansion(k, n, i, param):
    r"""
    Return the `i`-th coefficient of the expansion for the form
    `f_{k, n}`.

    INPUT:

    - ``k`` -- an integer representing the weight of the resulting form.
    - ``n`` -- an integer. The type of the resulting form will be congruent to
      this integer modulo `q-1`.
    - ``i`` -- an integer representing the index of the coefficient to compute.
    - ``param`` -- a prime power, a univariate polynomial ring over
      `\mathbb{F}_q` or the fraction field of such polynomial ring.

    EXAMPLES::

        sage: from drinfeld_modular_forms import coefficient_petrov_expansion
        sage: q = 3
        sage: coefficient_petrov_expansion(q + 1, 1, 1, q)
        1
        sage: coefficient_petrov_expansion(q + 1, 1, 0, q)
        0
        sage: coefficient_petrov_expansion(q + 1, 1, 1, q)
        1
        sage: coefficient_petrov_expansion(q + 1, 1, 7, q)
        2*T^3 + T

    .. doctest::
       :hide:

        sage: from drinfeld_modular_forms import coefficient_petrov_expansion
        sage: q = 5
        sage: A = GF(q)['T']
        sage: coefficient_petrov_expansion(q + 1, 1, 1, A)
        1
        sage: coefficient_petrov_expansion(q + 1, 1, 1, Frac(A))
        1
    """
    if k not in ZZ:
        raise TypeError("k must be an integer")
    if n not in ZZ:
        raise TypeError("n must be an integer")
    if i not in ZZ:
        raise TypeError("i must be an integer")
    q, polynomial_ring = _format_param(param)
    i = ZZ(i)
    n_th_goss_pol = _carlitz_module_goss_polynomial(n, polynomial_ring)
    m, G_m = next(((d, c) for d, c in enumerate(n_th_goss_pol.list()) if c))

    # compute part 1
    part1 = polynomial_ring.zero()
    if not i%m:
        if ZZ(i/m).is_power_of(q):
            d = ZZ(i/m).log(q)
            part1 = sum(a**(k - n)*G_m for a in polynomial_ring.monics(of_degree=d))

    # compute part 2
    part2 = polynomial_ring.zero()
    for d in range(1, i+1, 1):
        s = polynomial_ring.zero()
        if i >= (q**d + 1)*m:
            for a in polynomial_ring.monics(of_degree=d):
                G_ta = n_th_goss_pol.subs({n_th_goss_pol.parent().gen(): parameter_at_infinity(a)})
                dn = G_ta[i]
                s += a**(k - n)*dn
        part2 += s
    return part1 + part2

def compute_petrov_expansion(k, n, param, name='t', check=True):
    r"""
    Return the `A`-expansion of the form `f_{k, n}` as a lazy power
    series.

    Recall that it is defined by:

    .. MATH::

        f_{k, n} =
        \sum_{\substack{a\in \mathbb{F}_q[T] \\ a\text{ monic}}}
        a^{k - n}G_n(t(az))

    INPUT:

    - ``k`` -- an integer equal to the weight of the resulting form.
    - ``n`` -- an integer congruent to the type modulo `q-1`.
    - ``param`` -- a prime power, a univariate polynomial ring over
      `\mathbb{F}_q` or the fraction field of such polynomial ring.
    - ``name`` (string, default: ``'t'``) -- represent the parameter at
      infinity.
    - ``check`` (Boolean, default: ``True``) -- If this parameter is set
      to ``True``, the code will check if  `k` and `n` verify the three
      conditions:

      - `k - 2n > 0`;

      - `q - 1` divides `k - 2n`;

      - `n > p^{v_p(k - n)}`.

    OUTPUT: A lazy power series.

    EXAMPLES::

        sage: from drinfeld_modular_forms import compute_petrov_expansion
        sage: q = 3
        sage: D = compute_petrov_expansion(q + 1, 1, q); D
        t + t^5 + ((2*T^3+T)*t^7) + O(t^8)

    To obtain more coefficient, one just need to take slices the series::

        sage: D[0:10]  # output the coefficients from i=0...9
        [0, 1, 0, 0, 0, 1, 0, 2*T^3 + T, 0, 1]
        sage: D[39]  # print the 39-th coefficient
        T^9 + 2*T^3
        sage: D[59]  # 59-th coefficient
        2*T^27 + T^9 + T^3 + 2*T

    .. doctest::
       :hide:

        sage: from drinfeld_modular_forms import compute_petrov_expansion
        sage: q = 7
        sage: A = GF(q)['T']; K = Frac(A)
        sage: compute_petrov_expansion(q + 1, 1, A)
        t + O(t^8)
        sage: compute_petrov_expansion(q + 1, 1, K)
        t + O(t^8)
    """
    if k not in ZZ:
        raise TypeError("k must be an integer")
    if n not in ZZ:
        raise TypeError("n must be a integer")
    q, polynomial_ring = _format_param(param)
    Fq = polynomial_ring.base_ring()
    if check:
        if k - 2*n <= 0:
            raise ValueError(f"k - 2n must be > 0 (current value: {k - 2*n})")
        if (k - 2*n)%(q - 1):
            raise ValueError("k - 2n must be a multiple of q - 1")
        p = Fq.characteristic()
        if n > p**((k-n).valuation(p)):
            raise ValueError(f"n (={n}) must be less than p^val(p, k - n) (={p**((k-n).valuation(p))})")
    f = lambda i: coefficient_petrov_expansion(k, n, i, polynomial_ring)
    R = LazyPowerSeriesRing(polynomial_ring.fraction_field(), name)
    return R(f, valuation=1)

def compute_delta_rank_2(param, name='t'):
    r"""
    Return the expansion of the normalized Drinfeld modular discriminant of
    rank 2.

    Recall that for `w\in \Omega^2(\mathbb{C}_{\infty})` the *Drinfeld modular
    discriminant* is the leading coefficient of the Drinfeld module

    .. MATH::

        \phi_{w} : T \mapsto T + g(z)\tau + \Delta(w)\tau^2

    corresponding to the lattice `\Lambda_w := wA + A`. The *normalized*
    discriminant is the form `\Delta_0(w)` such that its first nonzero
    coefficient is 1.

    INPUT:

    - ``param`` -- a prime power, a univariate polynomial ring over
      `\mathbb{F}_q` or the fraction field of such polynomial ring.
    - ``name`` (string, default: ``'t'``) -- represent the parameter at
      infinity.

    OUTPUT: a lazy power series over the base ring.

    .. NOTE::

        We have `\Delta_0(w) = f_{q^2 - 1, q - 1}`.

    EXAMPLES::

        sage: from drinfeld_modular_forms import compute_delta_rank_2
        sage: q = 3
        sage: D = compute_delta_rank_2(q); D
        t^2 + 2*t^6 + O(t^8)
        sage: D[38]
        T^18 + 2*T^12 + 2*T^10 + T^4 + 1
        sage: D[56]
        2*T^27 + 2*T^3 + 2*T

    .. doctest::
       :hide:

        sage: from drinfeld_modular_forms import compute_delta_rank_2
        sage: q = 9
        sage: A = GF(q)['T']; K = Frac(A)
        sage: compute_delta_rank_2(A)
        O(t^8)
        sage: compute_delta_rank_2(K)
        O(t^8)
    """
    q = _format_param(param)[0]
    return compute_petrov_expansion(q**2 - 1, q - 1, param, name, check=False)

def compute_eisentein_series_rank_2(param, name='t'):
    r"""
    Return the expansion fo the normalized Drinfeld Eisenstein series of
    weight `q-1`.

    Recall that this form is defined by:

    .. MATH::

        E_{q-1}(w) := \frac{T^q - T}{\tilde{\pi}^{q-1}} \sum_{(c, d)\in A^2\setminus 0} \frac{1}{(cw + d)^{q-1}}.

    where `\tilde{\pi}` is the Carlitz period.

    INPUT:

    - ``param`` -- a prime power, a univariate polynomial ring over
      `\mathbb{F}_q` or the fraction field of such polynomial ring.
    - ``name`` (string, default: ``'t'``) -- represent the parameter at
      infinity.

    OUTPUT: a lazy power series over the base ring.

    .. NOTE::

        We have `E_{q-1} = 1 - (T^q - T)f_{q - 1, q - 1}`.

    EXAMPLES::

        sage: from drinfeld_modular_forms import compute_eisentein_series_rank_2
        sage: q = 3
        sage: E = compute_eisentein_series_rank_2(q); E
        1 + ((2*T^3+T)*t^2) + O(t^7)
        sage: E[0:15]
        [1, 0, 2*T^3 + T, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*T^3 + T]

    .. doctest::
       :hide:

        sage: from drinfeld_modular_forms import compute_eisentein_series_rank_2
        sage: q = 4
        sage: A = GF(q)['T']; K = Frac(A)
        sage: compute_eisentein_series_rank_2(A)
        1 + ((T^4+T)*t^3) + O(t^7)
        sage: compute_eisentein_series_rank_2(K)
        1 + ((T^4+T)*t^3) + O(t^7)
    """
    q, polynomial_ring = _format_param(param)
    T = polynomial_ring.gen()
    K = polynomial_ring.fraction_field()
    b = K(T**q - T)
    return K.one() - b*compute_petrov_expansion(q - 1, q - 1, param, name, check=False)

def compute_h_rank_2(param, name='t'):
    r"""
    Return the expansion of the function `h` defined by
    `\Delta = h^{q-1}` and having an expansion of the form `t + O(t^2)`.

    INPUT:

    - ``param`` -- a prime power, a univariate polynomial ring over
      `\mathbb{F}_q` or the fraction field of such polynomial ring.
    - ``name`` (string, default: ``'t'``) -- represent the parameter at
      infinity.

    EXAMPLES::

        sage: from drinfeld_modular_forms import compute_h_rank_2
        sage: q = 3
        sage: h = compute_h_rank_2(q); h
        t + t^5 + ((2*T^3+T)*t^7) + O(t^8)
        sage: h[0:15]
        [0, 1, 0, 0, 0, 1, 0, 2*T^3 + T, 0, 1, 0, T^3 + 2*T, 0, T^6 + T^4 + T^2 + 1, 0]

    .. doctest::
       :hide:

        sage: from drinfeld_modular_forms import compute_h_rank_2
        sage: q = 2
        sage: A = GF(q)['T']; K = Frac(A)
        sage: compute_h_rank_2(A)
        t + t^2 + ((T^2+T+1)*t^3) + t^4 + ((T^4+T+1)*t^5) + ((T^4+T^2)*t^6) + ((T^6+T^5+T^4+T^3+1)*t^7) + O(t^8)
        sage: compute_h_rank_2(K)
        t + t^2 + ((T^2+T+1)*t^3) + t^4 + ((T^4+T+1)*t^5) + ((T^4+T^2)*t^6) + ((T^6+T^5+T^4+T^3+1)*t^7) + O(t^8)
    """
    q = _format_param(param)[0]
    return compute_petrov_expansion(q + 1, 1, param, name)

def j_invariant(param, name='t'):
    r"""
    Return the `j`-invariant defined by `j := g_1^{q+1}/g_2`.

    It is normalized so that `j = 1/t^{q - 1} + O(1)`.

    INPUT:

    - ``param`` -- a prime power, a univariate polynomial ring over
      `\mathbb{F}_q` or the fraction field of such polynomial ring.
    - ``name`` (string, default: ``'t'``) -- represent the parameter at
      infinity.

    EXAMPLES::

        sage: from drinfeld_modular_forms import j_invariant
        sage: q = 3
        sage: j_invariant(q)
        1/t^2 + (2*T^3+T) + t^2 + ((2*T^9+2*T^3+2*T)*t^4) + O(t^5)
        sage: q = 4
        sage: j_invariant(q)
        1/t^3 + (T^4+T) + O(t^4)

    .. doctest::
       :hide:

        sage: from drinfeld_modular_forms import j_invariant
        sage: q = 8
        sage: A = GF(q)['T']; K = Frac(A)
        sage: j_invariant(A)
        1/t^7 + O(1)
        sage: j_invariant(K)
        1/t^7 + O(1)
    """
    q = _format_param(param)[0]
    E = compute_eisentein_series_rank_2(q, name)
    D = compute_delta_rank_2(q, name)
    return E**(q+1)/D
