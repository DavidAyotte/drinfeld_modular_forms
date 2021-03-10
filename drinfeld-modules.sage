
#TODO: add documentation

class DrinfeldModule:
    def __init__(self, K, coef):
        self.function_field = K
        self.coefficients = coef

    def ring_of_regular_functions(self):
        r"""
        Return the subring of function regular at infinity of the defining function field.

        OUTPUT: The polynomial ring F_q[T] where F_q is the finite field of q elements.

        EXAMPLES:

            sage: K.<T> = FunctionField(GF(5, 'a'));
            sage: C = CarlitzModule(K);
            sage: C.ring_of_regular_functions()
            Univariate Polynomial Ring in T over Finite Field of size 5

        Changing the name for the generator of the function field will also change the name of the generator::

            sage: K.<U> = FunctionField(GF(13^2, 'b'));
            sage: D = DrinfeldModule(K, [1,1]);
            sage: D.ring_of_regular_functions()
            Univariate Polynomial Ring in U over Finite Field in b of size 13^2

        TESTS::

            sage: K.<T> = FunctionField(GF(5, 'a'));
            sage: C = CarlitzModule(K);
            sage: C.ring_of_regular_functions()
            Univariate Polynomial Ring in T over Finite Field of size 5
            sage: K.<U> = FunctionField(GF(13^2, 'b'));
            sage: D = DrinfeldModule(K, [1,1]);
            sage: D.ring_of_regular_functions()
            Univariate Polynomial Ring in U over Finite Field in b of size 13^2
        """
        return PolynomialRing(self.function_field.constant_base_field(), self.function_field.gen())

    def symbolic_additive_polynomials(self):
        return PolynomialRing(self.ring_of_regular_functions(), 'X')

    def rank(self):
        return len(self.coefficients)

    def characteristic(self):
        return self.function_field.constant_field().characteristic()

    def prime_power_of_constant_field(self):
        F = self.function_field.constant_base_field()
        return F.cardinality().ord(self.characteristic())

    def frobenius(self):
        return K.frobenius_endomorphism()^(self.prime_power_of_constant_field())

    def action_of_T(self, exp=1):
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
        A = self.ring_of_regular_functions()
        a = A(a)
        action = 0
        coeffs = a.list()
        for power, coeff in enumerate(coeffs):
            action += coeff*self.action_of_T(exp=power)
        return action
        
    def __repr__(self):
        return "Drinfeld module of rank %s for the %s defined by T |--> %s"%(self.rank(), self.function_field, str(self.action_of_T()))

    def __call__(self, a):
        return self.action(a)

class CarlitzModule(DrinfeldModule):
    def __init__(self, K):
        DrinfeldModule.__init__(self, K, [1])



    


