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
    sage: A = GF(3)['T']
    sage: D = compute_petrov_expansion(3+1, 1, A); D
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
from drinfeld_modular_forms.goss_polynomials import bracket, goss_polynomial

from functools import cache

from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ

lazy_import('sage.rings.lazy_series_ring', 'LazyPowerSeriesRing')

def inverse_cyclotomic_polynomial(a, name='X'):
    r"""
    Return the polynomial `f_a(X) = \rho_a(X^{-1}) X^{q^{\mathrm{deg}(a)}}`
    where `\rho_a` is the Carlitz module.

    This function is used in :func:`.parameter_at_infinity`

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
    - ``name`` (Str, default: 't') -- the name of the lazy power series
      ring generator

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

def coefficient_petrov_expansion(k, n, i, polynomial_ring):
    r"""
    Return the `i`-th coefficient of the expansion for the form
    `f_{k, n}`.

    INPUT:

    - ``k`` -- an integer representing the weight of the resulting form.
    - ``n`` -- an integer. The type of the resulting form will be congruent to
      this integer modulo `q-1`.
    - ``i`` -- an integer representing the index of the coefficient to compute.
    - ``polynomial_ring`` -- a univariate polynomial ring over a finite field.

    EXAMPLES::

        sage: from drinfeld_modular_forms import coefficient_petrov_expansion
        sage: A = GF(3)['T']
        sage: coefficient_petrov_expansion(3+1, 1, 1, A)
        1
        sage: coefficient_petrov_expansion(3+1, 1, 0, A)
        0
        sage: coefficient_petrov_expansion(3+1, 1, 1, A)
        1
        sage: coefficient_petrov_expansion(3+1, 1, 7, A)
        2*T^3 + T
    """
    if i not in ZZ:
        raise TypeError("n must be an integer")
    i = ZZ(i)
    Fq = polynomial_ring.base_ring()
    if not Fq.is_field() or not Fq.is_finite():
        raise ValueError("polynomials base ring must be a finite field")
    q = polynomial_ring.base_ring().cardinality()

    n_th_goss_pol = goss_polynomial(n, polynomial_ring)
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

@cache
def compute_petrov_expansion(k, n, polynomial_ring, name='t'):
    r"""
    Return the `A`-expansion of the form `f_{k, n}` as a lazy power series.

    INPUT:

    - ``k`` -- an integer representing the weight of the resulting form.
    - ``n`` -- an integer which will be congruent to the type modulo `q-1`.
    - ``polynomial_ring`` -- a univariate polynomial ring over a finite field.
    - ``name`` -- string (default: 't').

    OUTPUT: A lazy power series over the fraction field of ``polynomial_ring``.

    EXAMPLES::

        sage: from drinfeld_modular_forms import compute_petrov_expansion
        sage: A = GF(3)['T']
        sage: D = compute_petrov_expansion(3+1, 1, A); D
        t + t^5 + ((2*T^3+T)*t^7) + O(t^8)

    To obtain more coefficient, one just need to take slices the series::

        sage: D[0:10]  # output the coefficient from i=0...9
        [0, 1, 0, 0, 0, 1, 0, 2*T^3 + T, 0, 1]
        sage: D[39]  # print the 39-th coefficient
        T^9 + 2*T^3
        sage: D[59]  # 59-th coefficient
        2*T^27 + T^9 + T^3 + 2*T
    """
    f = lambda i: coefficient_petrov_expansion(k, n, i, polynomial_ring)
    R = LazyPowerSeriesRing(polynomial_ring.fraction_field(), name)
    return R(f, valuation=1)

def compute_delta_rank_2(polynomial_ring, name='t'):
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

    - ``polynomial_ring`` -- a univariate polynomial ring over a finite field.
    - ``name`` -- string (default: 't').

    OUTPUT: a lazy power series over the base ring.

    .. NOTE::

        We have `\Delta_0(w) = f_{q^2 - 1, q - 1}`.

    EXAMPLES::

        sage: from drinfeld_modular_forms import compute_delta_rank_2
        sage: A = GF(3)['T']
        sage: D = compute_delta_rank_2(A); D
        t^2 + 2*t^6 + O(t^8)
        sage: D[38]
        T^18 + 2*T^12 + 2*T^10 + T^4 + 1
        sage: D[56]
        2*T^27 + 2*T^3 + 2*T
    """
    q = polynomial_ring.base_ring().cardinality()
    return compute_petrov_expansion(q**2 - 1, q - 1, polynomial_ring, name)

def compute_eisentein_series_rank_2(polynomial_ring, name='t'):
    r"""
    Return the expansion fo the normalized Drinfeld Eisenstein series of
    weight `q-1`.

    Recall that this form is defined by:

    .. MATH::

        E_{q-1}(w) := \frac{T^q - T}{\tilde{\pi}^{q-1}} \sum_{(c, d)\in A^2\setminus 0} \frac{1}{(cw + d)^{q-1}}.

    where `\tilde{\pi}` is the Carlitz period.

    INPUT:

    - ``polynomial_ring`` -- a univariate polynomial ring over a finite field.
    - ``name`` -- string (default: 't').

    OUTPUT: a lazy power series over the base ring.

    .. NOTE::

        We have `E_{q-1} = 1 - (T^q - T)f_{q - 1, q - 1}`.

    EXAMPLES::

        sage: from drinfeld_modular_forms import compute_eisentein_series_rank_2
        sage: A = GF(3)['T']
        sage: E = compute_eisentein_series_rank_2(A); E
        1 + ((2*T^3+T)*t^2) + O(t^7)
        sage: E[0:15]
        [1, 0, 2*T^3 + T, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*T^3 + T]
    """
    T = polynomial_ring.gen()
    q = polynomial_ring.base_ring().cardinality()
    K = polynomial_ring.fraction_field()
    b = K(T**q - T)
    return K.one() - b*compute_petrov_expansion(q - 1, q - 1, polynomial_ring, name)

def compute_h_rank_2(polynomial_ring, name='t'):
    q = polynomial_ring.base_ring().cardinality()
    return compute_petrov_expansion(q + 1, 1, polynomial_ring, name)
