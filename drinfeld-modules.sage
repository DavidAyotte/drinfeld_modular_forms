
r"""
This sage module defines a class named DrinfeldModule

TODO: quick explanation of Drinfeld modules here.

EXAMPLES::

The class constructor must take in inputs a FunctionField objects over a finite field::

    sage: K.<T> = FunctionField(GF(3))
    sage: D = DrinfeldModule(K, [1,1]); D
    Drinfeld module of rank 2 over the Rational function field in T over Finite Field of size 3 defined by T |--> X^9 + X^3 + T*X

It is possible to work directly with the Carlitz module. In that case, there is no need in passing the coefficients to the class constructor::

    sage: K.<T> = FunctionField(GF(5))
    sage: C = CarlitzModule(K); C
    Drinfeld module of rank 1 over the Rational function field in T over Finite Field of size 5 defined by T |--> X^5 + T*X

The call method is implemented to compute the image of a given Drinfeld module at a polynomial with `F_q`-coefficients::

    sage: K.<T> = FunctionField(GF(3))
    sage: D = DrinfeldModule(K, [1,1])
    sage: D(T)                                                                                                                                   
    X^9 + X^3 + T*X
    sage: D(T^2+1)                                                                                                                            
    X^81 + 2*X^27 + (T^9 + T + 1)*X^9 + (T^3 + T)*X^3 + (T^2 + 1)*X

.. WARNING::

    Computing the action at a polynomial of degree >= 2 can be computationnally very demanding.



AUTHORS:

- DAVID AYOTTE (2021): initial version, https://github.com/DavidAyotte/drinfeld-modular-forms

This code is based on the original code written by Alex Petrov (see the file "petrov/AlexPetrov-original-code-drinfeld-modular-forms.sage")

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


class AdditivePolynomials:
    #TODO: implement this class
    pass

class DrinfeldModule:
    r"""
    This class defines a Drinfeld module.

    INPUT:

    - ``K`` -- A function field over a finite field.
    - ``coef`` -- A list [a1, a2,..., an] corresponding to to the additive polynomial T*X + a1*X^q + a2*X^(q^2) + ... + an*X^(q^n). 

    EXAMPLES:

        sage: K.<T> = FunctionField(GF(3))
        sage: D = DrinfeldModule(K, [1,2]); D
        Drinfeld module of rank 2 over the Rational function field in T over Finite Field of size 3 defined by T |--> 2*X^9 + X^3 + T*X
        sage: D(T)
        2*X^9 + X^3 + T*X
        sage: D(T^2+1)
        X^81 + X^27 + (2*T^9 + 2*T + 1)*X^9 + (T^3 + T)*X^3 + (T^2 + 1)*X
    """
    def __init__(self, K, coef):
        if K.constant_base_field().cardinality() == Infinity:
            raise ValueError("The constant base field must be a finite field")
        self.function_field = K
        self.coefficients = coef

    def ring_of_regular_functions(self):
        r"""
        Return the subring of function regular at infinity of the defining function field.

        OUTPUT: The polynomial ring F_q[T] where F_q is the finite field of q elements.

        EXAMPLES:

            sage: K.<T> = FunctionField(GF(5, 'a'))
            sage: C = CarlitzModule(K)
            sage: C.ring_of_regular_functions()
            Univariate Polynomial Ring in T over Finite Field of size 5

        Changing the name for the generator of the function field will also change the name of the generator::

            sage: K.<U> = FunctionField(GF(13^2, 'b'))
            sage: D = DrinfeldModule(K, [1,1])
            sage: D.ring_of_regular_functions()
            Univariate Polynomial Ring in U over Finite Field in b of size 13^2

        TESTS::

            sage: K.<T> = FunctionField(GF(5, 'a'))
            sage: C = CarlitzModule(K)
            sage: C.ring_of_regular_functions()
            Univariate Polynomial Ring in T over Finite Field of size 5
            sage: K.<U> = FunctionField(GF(13^2, 'b'))
            sage: D = DrinfeldModule(K, [1,1])
            sage: D.ring_of_regular_functions()
            Univariate Polynomial Ring in U over Finite Field in b of size 13^2
        """
        return PolynomialRing(self.function_field.constant_base_field(), self.function_field.gen())

    def symbolic_additive_polynomials(self):
        r"""
        Return the polynomial ring A[X] where A is the ring of functions regular at infinity. Here X is to be seen as the frobenius X: w |--> w^p.

        OUTPUT: The polynomial ring F_q[T][X] where F_q is the finite field of q-elements. Note that this method is for representation purpose and
        does not return the real ring of additive polynomials. Please use the call method of the class to compute the action.

        EXAMPLES:
        
            sage: K.<T> = FunctionField(GF(5, 'a'))
            sage: D = DrinfeldModule(K, [1,1])
            sage: D.symbolic_additive_polynomials()
            Univariate Polynomial Ring in X over Univariate Polynomial Ring in T over Finite Field of size 5
            sage: C = CarlitzModule(K)
            sage: C.symbolic_additive_polynomials()
            Univariate Polynomial Ring in X over Univariate Polynomial Ring in T over Finite Field of size 5

        TESTS::

            sage: K.<T> = FunctionField(GF(5, 'a'))
            sage: D = DrinfeldModule(K, [1,1])
            sage: D.symbolic_additive_polynomials()
            Univariate Polynomial Ring in X over Univariate Polynomial Ring in T over Finite Field of size 5
            sage: C = CarlitzModule(K)
            sage: C.symbolic_additive_polynomials()
            Univariate Polynomial Ring in X over Univariate Polynomial Ring in T over Finite Field of size 5

        .. TODO::

            Implement the real ring of additive polynomial (with composition law)
        """
        return PolynomialRing(self.ring_of_regular_functions(), 'X')

    def rank(self):
        r"""
        Return the rank of the Drinfeld module.

        OUTPUT: an integer >= 1 corresponding to the rank of the Drinfeld module. Let D a Drinfeld module over F_q[T]. Then the rank correspond
        to the degree in X^q of the additive polynomial D(T).

        EXAMPLES:

            sage: K.<T> = FunctionField(GF(11, 'a'))
            sage: C = CarlitzModule(K); C
            Drinfeld module of rank 1 over the Rational function field in T over Finite Field of size 11 defined by T |--> X^11 + T*X
            sage: C.rank()
            1
            sage: D2 = DrinfeldModule(K, [1,1])
            sage: D2.rank()
            2
            sage: D3 = DrinfeldModule(K, [1,1,1])
            sage: D3.rank()
            3
            sage: D = DrinfeldModule(K, [0,0,0,0,1])
            sage: D.rank()
            5
        
        TESTS::

            sage: K.<T> = FunctionField(GF(11, 'a'))
            sage: C = CarlitzModule(K); C
            Drinfeld module of rank 1 over the Rational function field in T over Finite Field of size 11 defined by T |--> X^11 + T*X
            sage: C.rank()
            1
            sage: D2 = DrinfeldModule(K, [1,1])
            sage: D2.rank()
            2
            sage: D3 = DrinfeldModule(K, [1,1,1])
            sage: D3.rank()
            3
            sage: D = DrinfeldModule(K, [0,0,0,0,1])
            sage: D.rank()
            5
        """
        return len(self.coefficients)

    def characteristic(self):
        r"""
        This method return the caracteristic of the function field.

        OUTPUT: An integer corresponding to the characteristic of the base function field, that is the characteristic of the constants
        base field F_q.

        EXAMPLES:

            sage: K.<T> = FunctionField(GF(7^3, 'a'))
            sage: C = CarlitzModule(K)
            sage: C.characteristic()
            7
            sage: K.characteristic()
            7
        """
        return self.function_field.characteristic()

    def prime_power_of_constant_field(self):
        r"""
        This method return the power of the prime characteristic of the constant field.

        OUTPUT: An integer n>1 corresponding such that if p is the characteristic of the defining function field, then p^n is the 
        size of the constant field.

        EXAMPLES:

            sage: K.<T> = FunctionField(GF(7^10, 'a'))
            sage: D = DrinfeldModule(K, [1,1,1])
            sage: D.prime_power_of_constant_field()
            10
            sage: K.constant_field()
            Finite Field in a of size 7^10
        """
        F = self.function_field.constant_base_field()
        return F.cardinality().ord(self.characteristic())

    def frobenius(self):
        return K.frobenius_endomorphism()^(self.prime_power_of_constant_field())

    def q(self):
        return self.characteristic()^(self.prime_power_of_constant_field())

    def action_of_T(self, exp=1):
        r"""
        This method return the action of the Drinfeld module at T^exp (T being the generator of the function field).

        INPUT:

        - ``exp`` -- integer (default: ``1``); must be non-negative. Correspond to the exponent of T.

        OUTPUT: A polynomial of two variables corresponding to the image of T^exp at the Drinfeld module.

        EXAMPLES:

            sage: K.<T> = FunctionField(GF(3, 'a'))
            sage: D = DrinfeldModule(K, [1,1])
            sage: D.action_of_T()
            X^9 + X^3 + T*X
            sage: D.action_of_T(2)
            X^81 + 2*X^27 + (T^9 + T + 1)*X^9 + (T^3 + T)*X^3 + T^2*X
            sage: D.action_of_T(-1)
            ---------------------------------------------------------------------------
            ValueError                                Traceback (most recent call last)
            <ipython-input-11-68b1cb80f9fd> in <module>
            ----> 1 D.action_of_T(-Integer(1))

            <string> in action_of_T(self, exp)

            ValueError: The exponent of T must be a non-negative integer

        .. WARNING::

            This method is computationnally very demanding even for small primes (e.g. primes > 7) with exponents >= 2.

        .. TODO:: 

            Investigate a more efficient algorithm for this method
        """
        if (exp < 0):
            raise ValueError("The exponent of T must be a non-negative integer")
        R = self.symbolic_additive_polynomials()
        Rbase = R.base_ring()
        X = R.gen()
        T = Rbase.gen()
        if exp == 0:
            return X
        p = self.characteristic()
        T_action = T*X
        for j,c in enumerate(self.coefficients):
            T_action += c*X^(p^((j+1)*(self.prime_power_of_constant_field())))
        if exp == 1:
            return T_action
        T_exp_action = T_action
        for i in range(1, exp):
            T_exp_action = T_exp_action.subs(T_action)
        return T_exp_action

    def action(self, a):
        r"""
        Compute the action of the Drinfeld module at a, an element of the function field.

        INPUT: 

        - ``a`` - a polynomial with coefficient in F_q

        OUTPUT: A polynomial of two variables corresponding to the image of `a` at the Drinfeld module.

        EXAMPLES:

            sage: K.<T> = FunctionField(GF(3, 'a'))
            sage: C = CarlitzModule(K)
            sage: C.action(T^2+1)
            X^9 + (T^3 + T)*X^3 + (T^2 + 1)*X
            sage: C.action(T^2+T+1)
            X^9 + (T^3 + T + 1)*X^3 + (T^2 + T + 1)*X

        This method is the same as using the call method of the Drinfeld module::

            sage: K.<T> = FunctionField(GF(3, 'a'))
            sage: C = CarlitzModule(K)
            sage: C(T)
            X^3 + T*X
            sage: C(T^2+1)
            X^9 + (T^3 + T)*X^3 + (T^2 + 1)*X

        
        """
        A = self.ring_of_regular_functions()
        a = A(a)
        action = 0
        coeffs = a.list()
        for power, coeff in enumerate(coeffs):
            action += coeff*self.action_of_T(exp=power)
        return action

    def j_invariant(self):
        if self.rank() != 2:
            raise ValueError("The j-invariant is only defined for rank 2 Drinfeld Module (current rank: %s)"%(self.rank()))
        g = self.coefficients[0]
        delta = self.coefficients[1]
        q = self.characteristic()^self.prime_power_of_constant_field()
        return g^(q-1)/delta

        
    def __repr__(self):
        return "Drinfeld module of rank %s over the %s defined by T |--> %s"%(self.rank(), self.function_field, str(self.action_of_T()))

    def __call__(self, a):
        return self.action(a)

class CarlitzModule(DrinfeldModule):
    #TODO add documentation
    def __init__(self, K):
        DrinfeldModule.__init__(self, K, [1])

    def Br(self, n):
        if n < 0:
            raise ValueError("The value of n must be a non-negative integer")
        T = self.function_field.gen()
        p = self.characteristic()
        q = p^self.prime_power_of_constant_field()
        return T^(q^n)-T

    def monic_pol_products(self, deg):
        if deg < 0:
            raise ValueError("deg must be a non-negative integer")
        if deg == 0:
            return self.function_field(1)
        prod = 1
        for j in range(1, deg+1):
            prod *= self.Br(j)^(self.q()^(deg-j))
        return prod

    def LCM_monic_pol(self, deg):
        if deg < 0:
            raise ValueError("deg must be non-negative integer")
        if deg == 0:
            return self.function_field(1)
        prod = 1
        for i in range(1, deg+1):
            prod *= self.Br(i)
        return prod

    def action(self, a):
        X = self.symbolic_additive_polynomials().gen()
        T = self.function_field.gen()
        f=a*X
        d = len(a.list())
        aux = (a^self.q() - a)/(self.Br(1))
        for i in range(1, d+1):
            f += aux*X^(self.q()^i)
            aux = (aux^self.q() - aux)/self.Br(i + 1)
        return f

    #TODO: add more methods relative to CarlitzModule