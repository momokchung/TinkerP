import atoms
import energi
import field
import empole3

def analysis(Energi, Atoms, Field):
    if Field.use_mpole:
        empole3.empole3(Energi, Atoms, Field)
    Energi.esum = Energi.em
