import os
import numpy as np
import unittest
import tconfig
import analysis
import atoms
import field
import energi
import utils
import glob

path = os.path.abspath(os.path.dirname(__file__))

class Test_empole3(unittest.TestCase):

    def test_empole3a(self):
        files = glob.glob(f"{path}/testfiles/empole3/*.xyz")
        files.sort()
        files = [f.split(".xyz")[0] for f in files]
        for f in files:
            key = f + ".key"
            xyz = f + ".xyz"
            out = f + ".out"

            at,ff = tconfig.setup(key, xyz)

            en = energi.Energi()

            analysis.analysis(en, at, ff)

            outLines = utils.readClean(out)

            for line in outLines:
                if "Atomic Multipoles " in line:
                    split = line.split()
                    em = float(split[2])
                    nem = int(split[3])
            
            eps = 1e-6
            self.assertTrue(abs(em - en.em) < eps)
            self.assertEqual(nem, en.nem)
