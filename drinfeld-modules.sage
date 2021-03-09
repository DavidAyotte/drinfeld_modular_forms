
#TODO: add documentation

class DrinfeldModule:
    def __init__(self, K, coef):
        self.function_field = K
        self.coefficients = coef

    def ring_of_regular_functions(self):
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

def CarlitzModule(K):
    return DrinfeldModule(K, [1])



    


