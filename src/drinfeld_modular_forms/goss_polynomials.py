r"""
Functions for computing Goss polynomials and exponentials of general
Drinfeld `\mathbb{F}_q[T]`-modules.

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
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.real_double import RDF
from sage.rings.integer_ring import ZZ

from sage.functions.log import log

from drinfeld_modular_forms.drinfeld_modules import DrinfeldModule
from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.lazy_series_ring', 'LazyPowerSeriesRing')

def bracket(n, K):
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
    q = K.base_ring().cardinality()
    T = K.gen()
    return T ** (q ** n) - T

def product_of_monic_polynomials(n, K):
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
    ans = K.one()
    q = K.base_ring().cardinality()
    if n == 0:
        return ans
    for j in range(1, n+1):
        ans = ans * bracket(j, K) ** (q ** (n - j))
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

def _compute_coefficient_log(coeffs, n):
    r"""
    Return the `n`-th coefficient of the logarithm of the Drinfeld
    module determined by the list ``coeffs``.

    INPUT:

    - ``coeffs`` (list) -- a list of coefficients in a field of
      finite characteristic representing a Drinfeld
      `\mathbb{F}_q[T]`-module.
    - ``n`` (integer) -- the coefficient index

    TESTS::

        sage: from drinfeld_modular_forms.goss_polynomials import _compute_coefficient_log
        sage: A = GF(3)['T']
        sage: K.<T> = Frac(A)
        sage: coeffs = [T, K.one()]  # Carlitz module
        sage: _compute_coefficient_log(coeffs, 1)
        1
        sage: _compute_coefficient_log(coeffs, 3)
        2/(T^3 + 2*T)
        sage: _compute_coefficient_log(coeffs, 3^2)
        1/(T^12 + 2*T^10 + 2*T^4 + T^2)
        sage: _compute_coefficient_log(coeffs, 3^3)
        2/(T^39 + 2*T^37 + 2*T^31 + T^29 + 2*T^13 + T^11 + T^5 + 2*T^3)
    """
    if n not in ZZ:
        raise TypeError("input must be an integer")
    n = ZZ(n)
    T = coeffs[0]
    K = T.parent()
    if n.is_zero():
        return K.zero()
    if n.is_one():
        return K.one()
    r = len(coeffs)
    q = K.base_ring().cardinality()
    if not n.is_power_of(q):
        return K.zero()
    c = K.zero()
    for i in range(n.log(q)):
        j = n.log(q) - i
        if j < r:
            c += _compute_coefficient_log(coeffs, q**i)*coeffs[j]**(q**i)
    return c/(T - T**n)

def drinfeld_logarithm(coeffs, name='z'):
    r"""
    Return the logarithm of the Drinfeld module defined by the list
    ``coeffs``.

    Warning: the method does not check if the list ``coeffs`` is valid.

    INPUT:

    - ``coeffs`` (list) -- a list of coefficients in a field of
      finite characteristic representing a Drinfeld
      `\mathbb{F}_q[T]`-module.
    - ``name`` (str, default: ``'z'``) -- the name of the lazy power
      series ring.

    OUTPUT: a lazy power series in ``name``.

    EXAMPLES::

        sage: from drinfeld_modular_forms import drinfeld_logarithm
        sage: A = GF(3)['T']
        sage: K.<T> = Frac(A)
        sage: coeffs = [T, K.one()]  # Carlitz module
        sage: drinfeld_logarithm(coeffs)
        z + ((2/(T^3+2*T))*z^3) + O(z^8)

    ::

        sage: coeffs = [T, T^5 + T^2 + 2, 1/T, T^3]
        sage: drinfeld_logarithm(coeffs)
        z + (((2*T^5+2*T^2+1)/(T^3+2*T))*z^3) + O(z^8)
    """
    # Add validity checks to the list ``coeffs``.
    if isinstance(coeffs, list):
        if len(coeffs) < 1:
            raise ValueError("input must be of length >= 1")
    else:
        raise TypeError("input must be a list representing a Drinfeld "
                        "module")
    L = LazyPowerSeriesRing(coeffs[0].parent(), name)
    log = lambda k: _compute_coefficient_log(coeffs, k)
    return L(log, valuation=1)

def _compute_coefficient_exp(coeffs, k):
    r"""
    Return the `n`-th coefficient of the exponential of the Drinfeld
    module determined by the list ``coeffs``.

    INPUT:

    - ``coeffs`` (list) -- a list of coefficients in a field of
      finite characteristic representing a Drinfeld
      `\mathbb{F}_q[T]`-module.
    - ``n`` (integer) -- the coefficient index

    TESTS::

        sage: from drinfeld_modular_forms.goss_polynomials import _compute_coefficient_exp
        sage: A = GF(3)['T']
        sage: K.<T> = Frac(A)
        sage: coeffs = [T, K.one()]  # Carlitz module
        sage: _compute_coefficient_exp(coeffs, 1)
        1
        sage: _compute_coefficient_exp(coeffs, 3)
        1/(T^3 + 2*T)
        sage: _compute_coefficient_exp(coeffs, 3^2)
        1/(T^18 + 2*T^12 + 2*T^10 + T^4)
        sage: _compute_coefficient_exp(coeffs, 3^3)
        1/(T^81 + 2*T^63 + 2*T^57 + 2*T^55 + T^39 + T^37 + T^31 + 2*T^13)
    """
    if k not in ZZ:
        raise TypeError("input must be an integer")
    k = ZZ(k)
    K = coeffs[0].parent()
    if k.is_zero():
        return K.zero()
    if k.is_one():
        return K.one()
    q = K.base_ring().cardinality()
    if not k.is_power_of(q):
        return K.zero()
    c = K.zero()
    for i in range(k.log(q)):
        j = k.log(q) - i
        c += _compute_coefficient_exp(coeffs, q**i)*_compute_coefficient_log(coeffs, q**j)**(q**i)
    return -c

def drinfeld_exponential(coeffs, name='z'):
    r"""
    Return the exponential of the Drinfeld module defined by the list
    ``coeffs``.

    Warning: the method does not check if the list ``coeffs`` is valid.

    INPUT:

    - ``coeffs`` (list) -- a list of coefficients in a field of
      finite characteristic representing a Drinfeld
      `\mathbb{F}_q[T]`-module.
    - ``name`` (str, default: ``'z'``) -- the name of the lazy power
      series ring.

    OUTPUT: a lazy power series in ``name``.

    EXAMPLES::

        sage: from drinfeld_modular_forms import drinfeld_exponential
        sage: A = GF(3)['T']
        sage: K.<T> = Frac(A)
        sage: coeffs = [T, K.one()]  # Carlitz module
        sage: drinfeld_exponential(coeffs)
        z + ((1/(T^3+2*T))*z^3) + O(z^8)

    ::

        sage: coeffs = [T, T^5 + T^2 + 2, 1/T, T^3]
        sage: drinfeld_exponential(coeffs)
        z + (((T^5+T^2+2)/(T^3+2*T))*z^3) + O(z^8)
    """
    # Add validity checks to the list ``coeffs``.
    L = LazyPowerSeriesRing(coeffs[0].parent(), name)
    exp = lambda k: _compute_coefficient_exp(coeffs, k)
    return L(exp, valuation=1)

def goss_polynomial(coeffs, n, name='X'):
    r"""
    Return the `n`-th Goss polynomial of the Drinfeld module defined by
    the list ``coeffs``.

    Warning: the method does not check if the list ``coeffs`` is valid.

    INPUT:

    - ``coeffs`` (list) -- a list of coefficients in a field of
      finite characteristic representing a Drinfeld
      `\mathbb{F}_q[T]`-module.
    - ``n`` (integer) -- the index of the Goss polynomial.
    - ``name`` (str, default: ``'X'``) -- the name of polynomial
      variable.

    OUTPUT: a univariate polynomial in ``name``.

    EXAMPLES::

        sage: from drinfeld_modular_forms import goss_polynomial
        sage: A = GF(3)['T']
        sage: K.<T> = Frac(A)
        sage: coeffs = [T, K.one()]  # Carlitz module
        sage: goss_polynomial(coeffs, 1)
        X
        sage: goss_polynomial(coeffs, 2)
        X^2
        sage: goss_polynomial(coeffs, 4)
        X^4 + (1/(T^3 + 2*T))*X^2
        sage: goss_polynomial(coeffs, 5)
        X^5 + (2/(T^3 + 2*T))*X^3
        sage: goss_polynomial(coeffs, 10)
        X^10 + (1/(T^3 + 2*T))*X^8 + (1/(T^6 + T^4 + T^2))*X^6 + (1/(T^9 + 2*T^3))*X^4 + (1/(T^18 + 2*T^12 + 2*T^10 + T^4))*X^2

    ::

        sage: coeffs = [T, 1/(T^2 + 1), T^10 + T^5 + 1, T^2 + 2]
        sage: goss_polynomial(coeffs, 10)
        X^10 + (1/(T^5 + 2*T))*X^8 + (1/(T^10 + T^6 + T^2))*X^6 + (1/(T^15 + 2*T^3))*X^4 + ((T^25 + 2*T^24 + 2*T^22 + T^21 + 2*T^20 + 2*T^19 + T^18 + T^17 + 2*T^15 + T^14 + T^13 + 2*T^8 + T^7 + T^5 + 2*T^4 + T^3 + T^2 + 2*T + 2)/(T^24 + 2*T^23 + 2*T^21 + T^20 + T^19 + T^17 + T^16 + 2*T^12 + T^11 + T^9 + 2*T^8 + 2*T^7 + 2*T^5 + 2*T^4))*X^2
    """
    if n not in ZZ:
        raise TypeError("n must be an integer")
    if not isinstance(coeffs, list):
        raise TypeError("coeffs must be a list")
    if len(coeffs) <= 1:
        raise ValueError("coeffs must be of length > 1")
    n = ZZ(n)
    K = coeffs[0].parent()
    if len(coeffs) == 2 and coeffs[1].is_one():
        return _carlitz_module_goss_polynomial(n, K, name)
    R = K[name]
    X = R.gen()
    if n.is_zero():
        return R.zero()
    q = K.base_ring().cardinality()
    if n <= q - 1:
        return X**n
    if n%q == 0:
        return goss_polynomial(coeffs, ZZ(n/q))**q
    exp = drinfeld_exponential(coeffs)
    pol = sum(exp[q**(i+1)]*goss_polynomial(coeffs, n - q**(i+1))
              for i in range(0, (n.log(q).n()).floor()))
    return X*(goss_polynomial(coeffs, n - 1) + pol)

def _carlitz_module_goss_polynomial(n, K, name='X'):
    r"""
    Return the `n`-th Goss polynomial for the Carlitz module.

    EXAMPLES::

        sage: from drinfeld_modular_forms.goss_polynomials import _carlitz_module_goss_polynomial
        sage: A.<T> = GF(3)['T']
        sage: _carlitz_module_goss_polynomial(1, A)
        X
        sage: _carlitz_module_goss_polynomial(2, A)
        X^2
        sage: _carlitz_module_goss_polynomial(3, A)
        X^3
        sage: _carlitz_module_goss_polynomial(4, A)
        X^4 + (1/(T^3 + 2*T))*X^2
        sage: _carlitz_module_goss_polynomial(5, A)
        X^5 + (2/(T^3 + 2*T))*X^3
        sage: _carlitz_module_goss_polynomial(6, A)
        X^6
    """
    if n not in ZZ:
        raise TypeError("n must be an integer")
    n = ZZ(n)
    q = K.base_ring().cardinality()
    R = K[name]
    X = R.gen()
    pol = R.zero()
    if n.is_zero():
        return pol
    if n <= q - 1:
        return X ** n
    if (n % q) == 0:
        return _carlitz_module_goss_polynomial(n/q, K, name) ** q
    for j in range(0, n.log(q).n().floor() + 1):
        pol += (X*_carlitz_module_goss_polynomial(n - q ** j, K, name)
                     *(1/product_of_monic_polynomials(j, K)))
    return pol
