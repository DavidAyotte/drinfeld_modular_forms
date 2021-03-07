
#TODO: add documentation

class DrinfeldModule:
    def __init__(self, K, coef):
        self.function_field = K
        self.coefficients = coef

    def ring_of_regular_functions(self):
        return PolynomialRing(self.function_field.constant_base_field(), self.function_field.gen())

    def symbolic_additive_polynomials(self):
        return PolynomialRing(self.ring_of_regular_functions(), 'u')

    def rank(self):
        return len(self.coefficients)

    def prime_power_of_constant_field(self):
        F = self.function_field.constant_field()
        p = F.characteristic()
        return F.cardinality().ord(p)

    def frobenius(self):
        return K.frobenius_endomorphism()^constant_field_power(K)

    def action(self, a):
        #TODO: study the coefficients of a to define the action correctly
        tau = self.frobenius()
        print(tau)
        coef = self.coefficients
        gen = self.function_field.gen()
        symbolic_Frob = self.symbolic_additive_polynomials().gen()
        print(coef[0])
        f = a*symbolic_Frob
        for i in [0..len(coef)-1]:
            tau_power = tau^(i+1)
            print(tau_power)
            f += coef[i]*tau_power(gen)
        return f

    


