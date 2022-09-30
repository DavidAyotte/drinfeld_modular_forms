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

@cache
def Aexpansion(C, k, n, max_deg, prec, name='t'):
    r"""
    Return the A-expansion `\sum a^(k-n) G_n (t_a)` by using polynomials of
    degree less than or equal to ``max_deg`` and up to precision ``prec``.

    EXAMPLES::

        sage: from drinfeld_modules.all import *
        sage: A.<T> = GF(3)['T']
        sage: C = CarlitzModule(A)
        sage: Aexpansion(C, 4, 1, 3, 10)
        t + t^5 + (2*T^3 + T)*t^7 + t^9 + (T^3 + 2*T)*t^11 + O(t^13)
    """
    q = C.q()
    F = C.field_of_constants()
    A = C.base_polynomial_ring()
    if prec >= q ** (max_deg + 1):
        print("Warning: Cannot ensure the required precision. Choose ``max_deg``, so that ``prec <= q^(max_deg + 1)``.")
    characteristic = F.characteristic()
    if n > characteristic ** ((k-n).valuation(characteristic)):
        print("Warning: The computed A-expansion is not a modular form.")
    gn = C.goss_polynomial(n)
    normed_goss_pol = gn.coefficients()[0].inverse_of_unit() * gn
    ans = 0
    for j in range(0, max_deg + 1):
        for a in A.monics(of_degree=j):
            ans = ans + a ** (k-n)*normed_goss_pol(C.ta(a, prec, name))
    return ans

def compute_delta(C, max_deg, prec, name='t'):
    q = C.q()
    return Aexpansion(C, q ** 2 - 1, q - 1, max_deg, prec)
