r"""
Goss polynomials for the Carlitz module

AUTHORS:

- David Ayotte (2022): initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 DAVID AYOTTE <davidayotte94@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.real_double import RDF

from sage.functions.log import log

from drinfeld_modular_forms.drinfeld_modules import DrinfeldModule

def bracket(n, polynomial_ring):
    r"""
    Return the element `[n] = T^{q^n} - T` where `T` is the generator of the
    polynomial ring.

    EXAMPLES::

        sage: from drinfeld_modular_forms import bracket
        sage: A.<T> = GF(3)['T']
        sage: bracket(1, A)
        T^3 + 2*T
        sage: bracket(2, A)
        T^9 + 2*T
        sage: bracket(0, A)
        Traceback (most recent call last):
        ...
        ValueError: the integer n (=0) must be postive.
    """
    if n <= 0:
        raise ValueError(f"the integer n (={n}) must be postive.")
    q = polynomial_ring.base_ring().cardinality()
    T = polynomial_ring.gen()
    return T ** (q ** n) - T

def product_of_monic_polynomials(n, polynomial_ring):
    r"""
    Return the product of all monic polynomials of degree `n`.

    EXAMPLES::

        sage: from drinfeld_modular_forms import product_of_monic_polynomials
        sage: A.<T> = GF(3)['T']
        sage: product_of_monic_polynomials(0, A)
        1
        sage: product_of_monic_polynomials(1, A)
        T^3 + 2*T
        sage: f = product_of_monic_polynomials(2, A); f
        T^18 + 2*T^12 + 2*T^10 + T^4
        sage: f.factor()
        T^4 * (T + 1)^4 * (T + 2)^4 * (T^2 + 1) * (T^2 + T + 2) * (T^2 + 2*T + 2)
    """
    if n < 0:
        raise ValueError(f"the integer n (={n}) must be non-negative.")
    ans = polynomial_ring.base_ring().one()
    q = polynomial_ring.base_ring().cardinality()
    if n == 0:
        return ans
    for j in range(1, n+1):
        ans = ans * bracket(j, polynomial_ring) ** (q ** (n - j))
    return ans

def lcm_of_monic_polynomials(n, polynomial_ring):
    r"""
    Return the least common multiple of all monic polynomials of degree `n`.

    EXAMPLES::

        sage: from drinfeld_modular_forms import lcm_of_monic_polynomials
        sage: A.<T> = GF(3)['T']
        sage: lcm_of_monic_polynomials(1, A)
        T^3 + 2*T
        sage: lcm_of_monic_polynomials(2, A)
        T^12 + 2*T^10 + 2*T^4 + T^2
        sage: lcm_of_monic_polynomials(3, A)
        T^39 + 2*T^37 + 2*T^31 + T^29 + 2*T^13 + T^11 + T^5 + 2*T^3
    """
    if n < 0:
        raise ValueError(f"the integer n (={n}) must be non-negative.")
    ans = polynomial_ring.base_ring().one()
    if n == 0:
        return ans
    for j in range(1, n+1):
        ans = ans * bracket(j, polynomial_ring)
    return ans

def goss_polynomial(n, polynomial_ring):
    r"""
    Return the `n`-th Goss polynomial for the Carlitz module.

    EXAMPLES::

        sage: from drinfeld_modular_forms import goss_polynomial
        sage: A.<T> = GF(3)['T']
        sage: goss_polynomial(1, A)
        X
        sage: goss_polynomial(2, A)
        X^2
        sage: goss_polynomial(3, A)
        X^3
        sage: goss_polynomial(4, A)
        X^4 + (1/(T^3 + 2*T))*X^2
        sage: goss_polynomial(5, A)
        X^5 + (2/(T^3 + 2*T))*X^3
        sage: goss_polynomial(6, A)
        X^6
    """
    q = polynomial_ring.base_ring().cardinality()
    R = polynomial_ring['X']
    X = R.gen()
    pol = R.zero()
    if n == 0:
        return pol
    if n <= q - 1:
        return X ** n
    if (n % q) == 0:
        return goss_polynomial(n/q, polynomial_ring) ** q
    for j in range(0, RDF(log(n, q)).floor() + 1):
        pol = pol + X * goss_polynomial(n - q ** j, polynomial_ring) * (1/product_of_monic_polynomials(j, polynomial_ring))
    return pol
