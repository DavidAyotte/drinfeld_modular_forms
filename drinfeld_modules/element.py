from sage.structure.element import Element

class DrinfeldModularFormsRingElement(Element):

    def __init__(parent, polynomial):
        self.polynomial = polynomial

        Element.__init__(self, parent)
