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

def ta(a, name='t'):
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
    R = LazyPowerSeriesRing(A.fraction_field(), name)
    d = a.degree()
    fpol = R(inverse_cyclotomic_polynomial(a, name=name))
    u = R.gen()
    return u**(q**d)*fpol.inverse_of_unit()

def coefficient_A_expansion(k, nu, n, polynomial_ring):
    if n not in ZZ:
        raise TypeError("n must be an integer")
    n = ZZ(n)
    Fq = polynomial_ring.base_ring()
    if not Fq.is_field() or not Fq.is_finite():
        raise ValueError("polynomials base ring must be a finite field")
    q = polynomial_ring.base_ring().cardinality()

    Gnu = goss_polynomial(nu, polynomial_ring)
    coeff_G = Gnu.coefficients()[0]

    m, G_m = next(((d, c) for d, c in enumerate(Gnu.list()) if c))

    # compute part 1
    part1 = polynomial_ring.zero()
    if not n%m:
        if ZZ(n/m).is_power_of(q):
            d = ZZ(n/m).log(q)
            part1 = sum(a**(k - nu)*G_m for a in polynomial_ring.monics(of_degree=d))

    # compute part 2
    part2 = polynomial_ring.zero()
    for d in range(1, n+1, 1):
        s = polynomial_ring.zero()
        if n >= (q**d + 1)*m:
            for a in polynomial_ring.monics(of_degree=d):
                G_ta = Gnu.subs({Gnu.parent().gen(): ta(a)})
                dn = G_ta[n]
                s += a**(k - nu)*dn
        part2 += s
    return part1 + part2

def compute_A_expansion(k, nu, polynomial_ring, name='t'):
    f = lambda n: coefficient_A_expansion(k, nu, n, polynomial_ring)
    R = LazyPowerSeriesRing(polynomial_ring.fraction_field(), name)
    return R(f)

def compute_delta_rank_2(polynomial_ring, name='t'):
    q = polynomial_ring.base_ring().cardinality()
    return compute_A_expansion(q**2 - 1, q - 1, polynomial_ring, name)

def compute_eisentein_serie_rank_2(polynomial_ring, name='t'):
    T = polynomial_ring.gen()
    q = polynomial_ring.base_ring().cardinality()
    K = polynomial_ring.fraction_field()
    b = K(T**q - T)
    return b.inverse_of_unit() - compute_A_expansion(q - 1, q - 1, polynomial_ring, name)
