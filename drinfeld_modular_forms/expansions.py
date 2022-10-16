"""
Function for computing the A-expansion as defined by Petrov

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

from functools import cache
from drinfeld_modular_forms.drinfeld_modules import DrinfeldModule
from drinfeld_modular_forms.goss_polynomials import bracket, goss_polynomial
from sage.rings.power_series_ring import PowerSeriesRing

def inverse_cyclotomic_polynomial(a, name='X'):
    r"""
    Return the polynomial `f_a(X) = \rho_a(X^{-1}) X^{q^{\mathrm{deg}(a)}}`
    where `\rho_a` is the Carlitz module.

    This function is used in :func:`.ta`

    EXAMPLES::

        sage: from drinfeld_modular_forms.all import *
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

def ta(a, prec=10, name='t'):
    r"""
    Return the function `t(az)` as a power series in `t`.

    INPUT:

    - `a` -- univariate polynomial over a finite field
    - `prec` (Integer, default: 10) -- the precision for the power series ring
    - `name` (Str, default: 't') -- the name of the power series ring generator

    EXAMPLES::

        sage: from drinfeld_modular_forms.all import *
        sage: A.<T> = GF(3)['T']
        sage: ta(A.one())
        t
        sage: ta(T)
        t^3 + 2*T*t^5 + T^2*t^7 + 2*T^3*t^9 + T^4*t^11 + O(t^13)
        sage: ta(T^2)
        t^9 + (2*T^3 + 2*T)*t^15 + 2*T^2*t^17 + O(t^19)
    """
    # TODO: add polynomial type verification
    # TODO: fix default precision issue
    if not a:
        return ValueError("the polynomial must be non-zero")
    A = a.parent()
    q = A.base_ring().cardinality()
    R = PowerSeriesRing(A, name=name, default_prec=prec)
    d = a.degree()
    fpol = R(inverse_cyclotomic_polynomial(a, name=name))
    u = R.gen()
    return u ** (q ** d) * (1/fpol)

@cache
def compute_A_expansion(k, n, polynomial_ring, max_deg, prec, name='t'):
    r"""
    Return the A-expansion `\sum a^{k-n} G_n (t_a)` by using polynomials of
    degree less than or equal to ``max_deg`` and up to precision ``prec``.

    EXAMPLES::

        sage: from drinfeld_modular_forms.all import *
        sage: A.<T> = GF(3)['T']
        sage: compute_A_expansion(4, 1, A, 3, 10)
        t + t^5 + (2*T^3 + T)*t^7 + t^9 + (T^3 + 2*T)*t^11 + O(t^13)
    """
    F = polynomial_ring.base_ring()
    q = F.cardinality()
    if prec >= q ** (max_deg + 1):
        print("Warning: Cannot ensure the required precision. Choose ``max_deg``, so that ``prec <= q^(max_deg + 1)``.")
    # characteristic = F.characteristic()
    # if n > characteristic ** ((k-n).valuation(characteristic)):
    #     print("Warning: The computed A-expansion is not a modular form.")
    gn = goss_polynomial(n, polynomial_ring)
    normed_goss_pol = gn.coefficients()[0].inverse_of_unit() * gn
    ans = 0
    for j in range(0, max_deg + 1):
        for a in polynomial_ring.monics(of_degree=j):
            ans = ans + a ** (k-n)*normed_goss_pol(ta(a, prec, name))
    return ans

@cache
def compute_delta_rank_2(A, max_deg, prec, name='t'):
    r"""
    Return the `t`-expansion of of the Drinfeld modular discriminant of weight
    `q^2 - 1`.
    """
    q = A.base_ring().cardinality()
    return compute_A_expansion(q ** 2 - 1, q - 1, A, max_deg, prec)

@cache
def compute_eisentein_serie_rank_2(A, max_deg, prec, name='t'):
    r"""
    Return the `t`-expansion of the weight `q - 1` Eisenstein serie.
    """
    q = A.base_ring().cardinality()
    return bracket(1, A) ** (-1) - compute_A_expansion(q - 1, q - 1, A, max_deg, prec)
