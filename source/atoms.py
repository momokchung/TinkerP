import numpy as np
import utils
from sizes import *
from units import *
import field

class Atoms:
    def __init__(self, xyzfile):
        self.xyzfile = xyzfile

    def readxyz(self):
        xyzfile = self.xyzfile
        xyzLines = utils.readClean(xyzfile)
        n = 0
        n = int(xyzLines[0].split()[0])
        if n != int(xyzLines[-1].split()[0]):
            raise Exception("Number of atoms in title and xyz files inconsistent.")
        self.n = n
        #initialize coordinates and connectivities for each atom
        self.tag = -1*np.ones(n, dtype=int)
        self.name = np.array([" "*3]*n)
        self.x = np.zeros(n)
        self.y = np.zeros(n)
        self.z = np.zeros(n)
        self.type = np.zeros(n, dtype=int)
        self.n12 = np.zeros(n, dtype=int)
        self.i12 = -1*np.ones((n,maxval), dtype=int)
        # read the coordinates and connectivities for each atom
        lines = xyzLines[1].split()
        start = 2
        for line in lines:
            if not utils.isfloat(line):
                start = 1
        for i,line in enumerate(xyzLines[start:]):
            split = line.split()
            self.tag[i] = int(split[0])-1
            self.name[i] = split[1]
            self.x[i] = float(split[2])
            self.y[i] = float(split[3])
            self.z[i] = float(split[4])
            self.type[i] = int(split[5])
            l = len(split[6:])
            for j in range(l):
                self.i12[i,j] = int(split[j+6])-1
        # for each atom, count and sort its attached atoms
        for i in range(n):
            j = 0
            for k in self.i12[i]:
                if k == -1:
                    break
                j += 1
            self.n12[i] = j
            i12sub = self.i12[i][:j]
            i12sub.sort()

    def attach(self):
        maxn13 = 3 * maxval
        maxn14 = 9 * maxval
        maxn15 = 27 * maxval
        n = self.n
        self.n13 = np.zeros(n, dtype=int)
        self.n14 = np.zeros(n, dtype=int)
        self.n15 = np.zeros(n, dtype=int)
        self.i13 = -1*np.ones((n, maxn13), dtype=int)
        self.i14 = -1*np.ones((n, maxn14), dtype=int)
        self.i15 = -1*np.ones((n, maxn15), dtype=int)
        # loop over all atoms finding all the 1-3 relationships
        for i in range(n):
            for j in range(self.n12[i]):
                jj = self.i12[i,j]
                for k in range(self.n12[jj]):
                    kk = self.i12[jj,k]
                    if kk == i:
                        continue
                    found = False
                    for m in range(self.n12[i]):
                        if kk == self.i12[i,m]:
                            found = True
                            break
                    if found:
                        continue
                    self.n13[i] += 1
                    self.i13[i,self.n13[i]-1] = kk
            i13sub = self.i13[i][:self.n13[i]]
            i13sub.sort()
        # loop over all atoms finding all the 1-4 relationships
        for i in range(n):
            for j in range(self.n13[i]):
                jj = self.i13[i,j]
                for k in range(self.n12[jj]):
                    kk = self.i12[jj,k]
                    if kk == i:
                        goto_30 = True
                        break
                    for m in range(self.n12[i]):
                        if kk == self.i12[i,m]:
                            goto_30 = True
                            break
                    else:
                        for m in range(self.n13[i]):
                            if kk == self.i13[i,m]:
                                goto_30 = True
                                break
                        else:
                            self.n14[i] += 1
                            self.i14[i] [self.n14[i]-1]= kk
                if goto_30:
                    goto_30 = False
                    continue
            goto_30 = False
            i14sub = self.i14[i][:self.n14[i]]
            i14sub.sort()
        # loop over all atoms finding all the 1-5 relationships
        for i in range(n):
            for j in range(self.n14[i]):
                jj = self.i14[i,j]
                for k in range(self.n12[jj]):
                    kk = self.i12[jj,k]
                    if kk == i:
                        goto_50 = True
                        break
                    for m in range(self.n12[i]):
                        if kk == self.i12[i,m]:
                            goto_50 = True
                            break
                    else:
                        for m in range(self.n13[i]):
                            if kk == self.i13[i,m]:
                                goto_50 = True
                                break
                        else:
                            for m in range(self.n14[i]):
                                if kk == self.i14[i,m]:
                                    goto_50 = True
                                    break
                            else:
                                self.n15[i] += 1
                                self.i15[i, self.n15[i]-1] = kk
                if goto_50:
                    goto_50 = False
                    continue
                else:
                    goto_50 = False
            goto_50 = False
            i15sub = self.i15[i][:self.n15[i]]
            i15sub.sort()

    def katom(self, Field):
        n = self.n
        self.clss = np.zeros(maxatm, dtype=int)
        self.atomic = np.zeros(maxatm, dtype=int)
        self.valence = np.zeros(maxatm, dtype=int)
        self.mass = np.zeros(maxatm)
        self.name = np.array([" "*3]*maxatm)
        self.story = np.array([" "*24]*maxatm)
        # transfer atom type values to individual atoms
        for i in range(n):
            k = self.type[i]
            if k == 0:
                self.clss[i] = 0
                self.atomic[i] = 0
                self.mass[i] = 0.
                self.valence[i] = 0
                self.story[i] = "Undefined Atom Type     "
            else:
                if Field.symbol[k] != "   ":
                    self.name[i] = Field.symbol[k]
                self.clss[i] = Field.atmcls[k]
                self.atomic[i] = Field.atmnum[k]
                self.mass[i] = Field.weight[k]
                self.valence[i] = Field.ligand[k]
                self.story[i] = Field.describe[k]

    def kmpole(self,Field):
        self.maxpole = 13
        n = self.n
        self.ipole = -1*np.ones(n, dtype=int)
        self.polsiz = np.zeros(n, dtype=int)
        self.pollist = -1*np.ones(n, dtype=int)
        self.zaxis = -1*np.ones(n, dtype=int)
        self.xaxis = -1*np.ones(n, dtype=int)
        self.yaxis = -1*np.ones(n, dtype=int)
        self.pole = np.zeros((n,self.maxpole))
        self.rpole = np.zeros((n,self.maxpole))
        self.mono0 = np.zeros(n)
        self.polaxe = np.array([" "*8]*n)
        # count the number of existing multipole parameters
        blank = "                "
        maxnmp = Field.maxnmp
        nmp = maxnmp
        for i in range(maxnmp-1, -1, -1):
            if Field.kmp[i] == blank:
                nmp = i - 1
        mpt = np.zeros(maxnmp)
        mpz = np.zeros(maxnmp)
        mpx = np.zeros(maxnmp)
        mpy = np.zeros(maxnmp)
        # store the atom types associated with each parameter
        for i in range(nmp):
            mpt[i] = int(Field.kmp[i][:4])
            mpz[i] = int(Field.kmp[i][4:8])
            mpx[i] = int(Field.kmp[i][8:12])
            mpy[i] = int(Field.kmp[i][12:])
        # assign multipole parameters via only 1-2 connected atoms
        n = self.n
        goto140 = False
        for i in range(n):
            it = self.type[i]
            for imp in range(nmp):
                if it == mpt[imp]:
                    ztyp = mpz[imp]
                    xtyp = mpx[imp]
                    ytyp = mpy[imp]
                    for j in range(self.n12[i]):
                        ji = self.i12[i,j]
                        jt = self.type[ji]
                        if jt == ztyp:
                            for k in range(self.n12[i]):
                                ki = self.i12[i,k]
                                kt = self.type[ki]
                                if kt == xtyp and ki != ji:
                                    if ytyp == 0:
                                        self.pollist[i] = i
                                        self.zaxis[i] = ji
                                        self.xaxis[i] = ki
                                        self.polaxe[i] = Field.mpaxis[imp]
                                        for m in range(13):
                                            self.pole[i,m] = Field.multip[imp,m]
                                        goto140 = True
                                        break
                                    for l in range(self.n12[i]):
                                        li = self.i12[i,l]
                                        lt = self.type[li]
                                        if lt == ytyp and li != ji and li != ki:
                                            self.pollist[i] = i
                                            self.zaxis[i] = ji
                                            self.xaxis[i] = ki
                                            self.yaxis[i] = li
                                            self.polaxe[i] = Field.mpaxis[imp]
                                            for m in range(13):
                                                self.pole[i,m] = Field.multip[imp,m]
                                            goto140 = True
                                            break
                                    if goto140:
                                        break
                            if goto140:
                                break
                    if goto140:
                        break
            if goto140:
                goto140 = False
                continue
        # assign multipole parameters via 1-2 and 1-3 connected atoms
            for imp in range(nmp):
                if it == mpt[imp]:
                    ztyp = mpz[imp]
                    xtyp = mpx[imp]
                    ytyp = mpy[imp]
                    for j in range(self.n12[i]):
                        ji = self.i12[i,j]
                        jt = self.type[ji]
                        if jt == ztyp:
                            for k in range(self.n13[i]):
                                ki = self.i13[i,k]
                                kt = self.type[ki]
                                path = False
                                for m in range(self.n12[ki]):
                                    if self.i12[ki,m] == ji:
                                        path = True
                                if kt == xtyp and path:
                                    if ytyp == 0:
                                        self.pollist[i] = i
                                        self.zaxis[i] = ji
                                        self.xaxis[i] = ki
                                        self.polaxe[i] = Field.mpaxis[imp]
                                        for m in range(13):
                                            self.pole[i,m] = Field.multip[imp,m]
                                        goto140 = True
                                        break
                                    for l in range(self.n13[i]):
                                        li = self.i13[i,l]
                                        lt = self.type[li]
                                        path = False
                                        for m in range(self.n12[li]):
                                            if self.i12[li,m] == ji:
                                                path = True
                                        if lt == ytyp and li != ki and path:
                                            self.pollist[i] = i
                                            self.zaxis[i] = ji
                                            self.xaxis[i] = ki
                                            self.yaxis[i] = li
                                            self.polaxe[i] = Field.mpaxis[imp]
                                            for m in range(13):
                                                self.pole[i,m] = Field.multip[imp,m]
                                            goto140 = True
                                            break
                                    if goto140:
                                        break
                            if goto140:
                                break
                    if goto140:
                        break
            if goto140:
                goto140 = False
                continue
        # assign multipole parameters via only a z-defining atom
            for imp in range(nmp):
                if it == mpt[imp]:
                    ztyp = mpz[imp]
                    xtyp = mpx[imp]
                    ytyp = mpy[imp]
                    for j in range(self.n12[i]):
                        ji = self.i12[i,j]
                        jt = self.type[ji]
                        if jt == ztyp:
                            if xtyp == 0:
                                self.pollist[i] = i
                                self.zaxis[i] = ji
                                self.polaxe[i] = Field.mpaxis[imp]
                                for m in range(13):
                                    self.pole[i,m] = Field.multip[imp,m]
                                goto140 = True
                                break
                    if goto140:
                        break
            if goto140:
                goto140 = False
                continue
        # assign multipole parameters via no connected atoms
            for imp in range(nmp):
                if it == mpt[imp]:
                    ztyp = mpz[imp]
                    xtyp = mpx[imp]
                    ytyp = mpy[imp]
                    if ztyp == 0:
                        self.pollist[i] = i
                        self.polaxe[i] = Field.mpaxis[imp]
                        for m in range(13):
                            self.pole[i,m] = Field.multip[imp,m]
                        goto140 = True
                        break
            if goto140:
                goto140 = False
                continue
        #  convert the dipole and quadrupole moments to Angstroms, quadrupole divided by 3 for use as traceless values
        for i in range(n):
            for k in range(1,4):
                self.pole[i,k] = self.pole[i,k] * bohr
            for k in range(4,13):
                self.pole[i,k] = self.pole[i,k] * bohr**2 / 3.
        # get the order of the multipole expansion at each site
        n = self.n
        npole = n
        polmax = 0
        for i in range(n):
            size = 0
            for k in range(self.maxpole):
                if self.pole[i,k] != 0:
                    size = max(k+1,size)
            if size > 4:
                size = 13
            elif size > 1:
                size = 4
            self.polsiz[i] = size
            polmax = max(polmax,size)
        # initialize charge penetration parameters
        n = self.n
        self.pcore = np.zeros(n)
        self.pval = np.zeros(n)
        self.pval0 = np.zeros(n)
        self.palpha = np.zeros(n)
        # assign the charge penetration charge and alpha parameters 
        ncp = 0
        for i in range(n):
            self.pcore[i] = 0.
            self.pval[i] = self.pole[i,0]
            self.pval0[i] = self.pval[i]
            self.palpha[i] = 0.
            ic = self.clss[i]
            if ic != 0:
                self.pcore[i] = Field.cpele[ic]
                self.pval[i] = self.pole[i,0] - Field.cpele[ic]
                self.pval0[i] = self.pval[i]
                self.palpha[i] = Field.cpalp[ic]
        # remove zero or undefined electrostatic sites from the list
        if Field.use_mpole:
            npole = 0
            ncp = 0
            for i in range(n):
                if self.polsiz[i] != 0:
                    self.ipole[npole] = i
                    self.pollist[i] = npole
                    self.mono0[i] = self.pole[i,0]
                    if self.palpha[i] != 0.:
                        ncp = ncp + 1
                    npole = npole + 1
        self.npole = npole
        self.ncp = ncp
        # turn off atomic multipole potentials if not used
        if npole == 0:
            Field.use_mpole = False
        if ncp != 0:
            Field.use_chgpen = True
