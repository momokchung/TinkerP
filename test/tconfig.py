import atoms
import field
import getcart
import mechanic

def setup(keyfile, xyzfile):
    at = atoms.Atoms(xyzfile)
    ff = field.Field(keyfile)

    getcart.getcart(ff,at)
    mechanic.mechanic(ff, at)
    return at, ff