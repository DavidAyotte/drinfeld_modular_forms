
def regular_ring_of_functions(K):
    return PolynomialRing(K.constant_base_field(), 'T')

def constant_field_power(K):
    F = K.constant_field()
    p = F.characteristic()
    return F.cardinality().ord(p)

class DrinfeldModule:
    def __init__(self, K, coef):
        self.function_field = K
        self.coefficients = coef

    def generator_tau_polynomial(self):
        pass

    def rank(self):
        return len(self.coefficients)

    def frobenius(self):
        return K.frobenius_endomorphism()^constant_field_power(K)

    def generator_action(self):
        tau = self.frobenius()
        print(tau)
        coef = self.coefficients
        gen = K.gen()
        print(coef[0])
        f = gen
        for i in [0..len(coef)-1]:
            tau_power = tau^(i+1)
            print(tau_power)
            f += coef[i]*tau_power(gen)
        return f

    


