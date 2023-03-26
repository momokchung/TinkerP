import numpy as np
import os
import utils
from sizes import *

class Field:
    def __init__(self, keyfile):
        self.keyfile = keyfile
        self.readKeyLines()

    def readKeyLines(self):
        keyfile = self.keyfile
        self.keyLines = utils.readClean(keyfile)

    def field(self):
        self.initialize()
        self.getprm()
        self.getkey()

    def initialize(self):
        self.use_mpole = True

    def getprm(self):
        keyfile = self.keyfile
        keyLines = self.keyLines
        prmfile = ""
        for line in keyLines:
            first = line.split()[0].upper()
            if first == "PARAMETER" or first == "PARAMETERS":
                prmfile = line.split()[1]
        if ".prm" not in prmfile.split("/")[-1]:
            prmfile += ".prm"
        if not os.path.isfile(prmfile):
            cwd = os.path.dirname(os.path.realpath(keyfile))
            prmfile = cwd+"/"+prmfile
        self.prmfile = prmfile
        self.setprm()
        self.initprm()
        self.readprm()

    def setprm(self):
        prmfile = self.prmfile
        keyfile = self.keyfile
        self.prmLines = utils.readClean(prmfile)
        prmLines = self.prmLines
        keyLines = self.keyLines
        maxnmp = 0
        for line in prmLines:
            first = line.split()[0].upper()
            if first == "MULTIPOLE":
                maxnmp += 1
        for line in keyLines:
            first = line.split()[0].upper()
            if first == "MULTIPOLE":
                maxnmp += 1
        self.maxnmp = maxnmp
        # allocate atomic multipole forcefield parameters
        self.multip = np.zeros((maxnmp, 13))
        self.mpaxis = np.array([" "*8]*maxnmp)
        self.kmp = np.array([" "*16]*maxnmp)

    def initprm(self):
        # set default control parameters for atomic multipole terms
        self.pentyp = "GORDON1"
        self.m2scale = 0.
        self.m3scale = 0.
        self.m4scale = 1.
        self.m5scale = 1.
        self.use_chgpen = False
        # allocate parameters
        self.atmcls = np.zeros(maxtyp, dtype=int)
        self.atmnum = np.zeros(maxtyp, dtype=int)
        self.ligand = np.zeros(maxtyp, dtype=int)
        self.weight = np.zeros(maxtyp)
        self.symbol = np.array([" "*3]*maxtyp)
        self.describe = np.array([" "*24]*maxtyp)
        self.cpele = np.zeros(maxclass)
        self.cpalp = np.zeros(maxclass)

    def readprm(self):
        nmp = 0
        prmfile = self.prmfile
        prmLines = self.prmLines
        nprm = len(prmLines)
        iprm = 0
        while (iprm < nprm):
            line = prmLines[iprm]
            self.prmkey(line)
            split = line.split()
            first = split[0].upper()
            if first == "ATOM":
                if "\"" in line:
                    split = line.split("\"")
                    split0 = split[0].split()
                    split1 = [split[1]]
                    split2 = split[2].split()
                    split = split0 + split1 + split2
                elif "\'" in line:
                    split = line.split("\'")
                    split0 = split[0].split()
                    split1 = [split[1]]
                    split2 = split[2].split()
                    split = split0 + split1 + split2
                ia = 0
                clss = 0
                atn = 0
                wght = 0.
                lig = 0
                ia = int(split[1])
                clss = int(split[2])
                if clss == 0:
                    clss = ia
                self.atmcls[ia] = clss
                if ia != 0:
                    self.symbol[ia] = split[3]
                    self.describe[ia] = split[4]
                    self.atmnum[ia] = int(split[5])
                    self.weight[ia] = float(split[6])
                    self.ligand[ia] = int(split[7])
            elif first == "MULTIPOLE":
                ia = 0
                ib = 0
                ic = 0
                id = 0
                axt = "Z-then-X"
                pl = np.zeros(13)
                split = line.split()[1:]
                lenSplit = len(split)
                if lenSplit == 5:
                    ia = int(split[0])
                    ib = int(split[1])
                    ic = int(split[2])
                    id = int(split[3])
                    pl[0] = float(split[4])
                elif lenSplit == 4:
                    ia = int(split[0])
                    ib = int(split[1])
                    ic = int(split[2])
                    pl[0] = float(split[3])
                elif lenSplit == 3:
                    ia = int(split[0])
                    ib = int(split[1])
                    pl[0] = float(split[2])
                elif lenSplit == 2:
                    ia = int(split[0])
                    pl[0] = float(split[1])
                iprm += 1
                line = prmLines[iprm]
                split = line.split()
                split = [float(s) for s in split]
                pl[1] = split[0]
                pl[2] = split[1]
                pl[3] = split[2]
                iprm += 1
                line = prmLines[iprm]
                split = line.split()
                split = [float(s) for s in split]
                pl[4] = split[0]
                iprm += 1
                line = prmLines[iprm]
                split = line.split()
                split = [float(s) for s in split]
                pl[7] = split[0]
                pl[8] = split[1]
                iprm += 1
                line = prmLines[iprm]
                split = line.split()
                split = [float(s) for s in split]
                pl[10] = split[0]
                pl[11] = split[1]
                pl[12] = split[2]
                if ib == 0:
                    axt = "None"
                if ib != 0 and ic == 0:
                    axt = "Z-Only"
                if ib < 0 or ic < 0:
                    axt = "Bisector"
                if ic < 0 and id < 0:
                    axt = "Z-Bisect"
                if max(ib,ic,id) < 0:
                    axt = "3-Fold"
                ib = abs(ib)
                ic = abs(ic)
                id = abs(id)
                pa = str(ia).zfill(4)
                pb = str(ib).zfill(4)
                pc = str(ic).zfill(4)
                pd = str(id).zfill(4)
                self.kmp[nmp] = pa+pb+pc+pd
                self.mpaxis[nmp] = axt
                self.multip[nmp,0] = pl[0]
                self.multip[nmp,1] = pl[1]
                self.multip[nmp,2] = pl[2]
                self.multip[nmp,3] = pl[3]
                self.multip[nmp,4] = pl[4]
                self.multip[nmp,5] = pl[7]
                self.multip[nmp,6] = pl[10]
                self.multip[nmp,7] = pl[7]
                self.multip[nmp,8] = pl[8]
                self.multip[nmp,9] = pl[11]
                self.multip[nmp,10] = pl[10]
                self.multip[nmp,11] = pl[11]
                self.multip[nmp,12] = pl[12]
                nmp += 1
            elif first == "CHGPEN":
                ia = 0
                pel = 0.
                pal = 0.
                line = prmLines[iprm]
                split = line.split()[1:]
                ia = int(split[0])
                pel = float(split[1])
                pal = float(split[2])
                if ia != 0:
                    self.cpele[ia] = abs(pel)
                    self.cpalp[ia] = pal
            iprm += 1

    def prmkey(self, line):
        split = line.split()
        first = split[0].upper()
        # set control parameters for atomic multipole potentials
        if first == "PENETRATION":
            self.pentyp = split[1].upper()
        elif first == "MPOLE-12-SCALE":
            self.m2scale = float(split[1])
            if self.m2scale > 1.:
                self.m2scale = 1. / self.m2scale
        elif first == "MPOLE-13-SCALE":
            self.m3scale = float(split[1])
            if self.m3scale > 1.:
                self.m3scale = 1. / self.m3scale
        elif first == "MPOLE-14-SCALE":
            self.m4scale = float(split[1])
            if self.m4scale > 1.:
                self.m4scale = 1. / self.m3scale
        elif first == "MPOLE-15-SCALE":
            self.m5scale = float(split[1])
            if self.m5scale > 1.:
                self.m5scale = 1. / self.m5scale

    def getkey(self):
        keyfile = self.keyfile
        keyLines = self.keyLines
        for line in keyLines:
            self.prmkey(line)

    def basefile(self):
        self.control()

    def control(self):
        digits = 4
        keyLines = self.keyLines
        for line in keyLines:
            split = line.split()
            if split[0].upper() == "DIGITS":
                digits = int(split[1])
        self.digits = digits
