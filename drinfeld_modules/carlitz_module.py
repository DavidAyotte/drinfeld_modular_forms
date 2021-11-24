from drinfeld_modules import DrinfeldModule

from sage.rings.all import FiniteField

class CarlitzModule(DrinfeldModule):
    def __init__(self, polynomial_base_ring, name='𝜏'):
        one = polynomial_base_ring.one()
        DrinfeldModule.__init__(self, one, name=name)

    def _repr_(self):
        return "Carlitz module over %s defined by:\n%s |--> %s" % (self.ring_of_constants(), self.base_ring().gen(), self.operator_polynomial())
