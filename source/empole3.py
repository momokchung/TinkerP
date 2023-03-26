import atoms
import energi
import field
import numpy as np
import chkpole
import rotpole
import math
import damping
from units import *

def empole3(Energi, Atoms, Field):
    empole3a(Energi, Atoms, Field)

def empole3a(Energi, Atoms, Field):
    if Atoms.npole == 0:
        return
    # zero out total atomic multipole energy and partitioning
    nem = 0
    em = 0
    n = Atoms.n
    npole = Atoms.npole
    use_chgpen = Field.use_chgpen
    pentyp = Field.pentyp
    aem = np.zeros(n)
    dmpi = np.zeros(9)
    dmpk = np.zeros(9)
    dmpik = np.zeros(9)
    # check the sign of multipole components at chiral sites
    chkpole.chkpole()
    # rotate the multipole components into the global frame
    rotpole.rotpole(Atoms, 'MPOLE')
    # set conversion factor, cutoff and switching coefficients
    f = electric / dielec
    mode = "MPOLE"
    m2scale = Field.m2scale
    m3scale = Field.m3scale
    m4scale = Field.m4scale
    m5scale = Field.m5scale
    mscale = np.ones(n)
    # calculate the multipole interaction energy term
    for ii in range(npole-1):
        i = Atoms.ipole[ii]
        iz = Atoms.zaxis[i]
        ix = Atoms.xaxis[i]
        iy = abs(Atoms.yaxis[i])
        xi = Atoms.x[i]
        yi = Atoms.y[i]
        zi = Atoms.z[i]
        ci = Atoms.rpole[i,0]
        dix = Atoms.rpole[i,1]
        diy = Atoms.rpole[i,2]
        diz = Atoms.rpole[i,3]
        qixx = Atoms.rpole[i,4]
        qixy = Atoms.rpole[i,5]
        qixz = Atoms.rpole[i,6]
        qiyy = Atoms.rpole[i,8]
        qiyz = Atoms.rpole[i,9]
        qizz = Atoms.rpole[i,12]
        if use_chgpen:
            corei = Atoms.pcore[i]
            vali = Atoms.pval[i]
            alphai = Atoms.palpha[i]
        # set exclusion coefficients for connected atoms
        for j in range(Atoms.n12[i]):
            mscale[Atoms.i12[i,j]] = m2scale
        for j in range(Atoms.n13[i]):
            mscale[Atoms.i13[i,j]] = m3scale
        for j in range(Atoms.n14[i]):
            mscale[Atoms.i14[i,j]] = m4scale
        for j in range(Atoms.n15[i]):
            mscale[Atoms.i15[i,j]] = m5scale
        # evaluate all sites within the cutoff distance
        for kk in range(ii+1, npole):
            k = Atoms.ipole[kk]
            kz = Atoms.zaxis[k]
            kx = Atoms.xaxis[k]
            ky = abs(Atoms.yaxis[k])
            proceed = True
            if proceed:
                xr = Atoms.x[k] - xi
                yr = Atoms.y[k] - yi
                zr = Atoms.z[k] - zi
                r2 = xr*xr + yr* yr + zr*zr
                r = math.sqrt(r2)
                ck = Atoms.rpole[k,0]
                dkx = Atoms.rpole[k,1]
                dky = Atoms.rpole[k,2]
                dkz = Atoms.rpole[k,3]
                qkxx = Atoms.rpole[k,4]
                qkxy = Atoms.rpole[k,5]
                qkxz = Atoms.rpole[k,6]
                qkyy = Atoms.rpole[k,8]
                qkyz = Atoms.rpole[k,9]
                qkzz = Atoms.rpole[k,12]
                # intermediates involving moments and separation distance
                dir_ = dix*xr + diy*yr + diz*zr
                qix = qixx*xr + qixy*yr + qixz*zr
                qiy = qixy*xr + qiyy*yr + qiyz*zr
                qiz = qixz*xr + qiyz*yr + qizz*zr
                qir = qix*xr + qiy*yr + qiz*zr
                dkr = dkx*xr + dky*yr + dkz*zr
                qkx = qkxx*xr + qkxy*yr + qkxz*zr
                qky = qkxy*xr + qkyy*yr + qkyz*zr
                qkz = qkxz*xr + qkyz*yr + qkzz*zr
                qkr = qkx*xr + qky*yr + qkz*zr
                dik = dix*dkx + diy*dky + diz*dkz
                qik = qix*qkx + qiy*qky + qiz*qkz
                diqk = dix*qkx + diy*qky + diz*qkz
                dkqi = dkx*qix + dky*qiy + dkz*qiz
                qiqk = 2.*(qixy*qkxy+qixz*qkxz+qiyz*qkyz) \
                         + qixx*qkxx+qiyy*qkyy+qizz*qkzz
                # get reciprocal distance terms for this interaction
                rr1 = f * mscale[k] / r
                rr3 = rr1 / r2
                rr5 = 3. * rr3 / r2
                rr7 = 5. * rr5 / r2
                rr9 = 7. * rr7 / r2
                # find damped multipole intermediates and energy value
                if use_chgpen:
                    corek = Atoms.pcore[k]
                    valk = Atoms.pval[k]
                    alphak = Atoms.palpha[k]
                    term1 = corei*corek
                    term1i = corek*vali
                    term2i = corek*dir_
                    term3i = corek*qir
                    term1k = corei*valk
                    term2k = -corei*dkr
                    term3k = corei*qkr
                    term1ik = vali*valk
                    term2ik = valk*dir_ - vali*dkr + dik
                    term3ik = vali*qkr + valk*qir - dir_*dkr \
                                 + 2.*(dkqi-diqk+qiqk)
                    term4ik = dir_*qkr - dkr*qir - 4.*qik
                    term5ik = qir*qkr

                    damping.damppole(r,9,alphai,alphak,dmpi,dmpk,dmpik,pentyp)
                    rr1i = dmpi[0]*rr1
                    rr3i = dmpi[2]*rr3
                    rr5i = dmpi[4]*rr5
                    rr1k = dmpk[0]*rr1
                    rr3k = dmpk[2]*rr3
                    rr5k = dmpk[4]*rr5
                    rr1ik = dmpik[0]*rr1
                    rr3ik = dmpik[2]*rr3
                    rr5ik = dmpik[4]*rr5
                    rr7ik = dmpik[6]*rr7
                    rr9ik = dmpik[8]*rr9
                    e = term1*rr1 + term1i*rr1i \
                          + term1k*rr1k + term1ik*rr1ik \
                          + term2i*rr3i + term2k*rr3k \
                          + term2ik*rr3ik + term3i*rr5i \
                          + term3k*rr5k + term3ik*rr5ik \
                          + term4ik*rr7ik + term5ik*rr9ik
                # find standard multipole intermediates and energy value
                else:
                    term1 = ci*ck
                    term2 = ck*dir_ - ci*dkr + dik
                    term3 = ci*qkr + ck*qir - dir_*dkr \
                               + 2.*(dkqi-diqk+qiqk)
                    term4 = dir_*qkr - dkr*qir - 4.*qik
                    term5 = qir*qkr
                    e = term1*rr1 + term2*rr3 + term3*rr5 \
                           + term4*rr7 + term5*rr9
                # increment the overall multipole energy components
                if e != 0.:
                    nem = nem + 1
                    em = em + e
                    aem[i] = aem[i] + 0.5*e
                    aem[k] = aem[k] + 0.5*e
        # reset exclusion coefficients for connected atoms
        for j in range(Atoms.n12[i]):
           mscale[Atoms.i12[i,j]] = 1.
        for j in range(Atoms.n13[i]):
           mscale[Atoms.i13[i,j]] = 1.
        for j in range(Atoms.n14[i]):
           mscale[Atoms.i14[i,j]] = 1.
        for j in range(Atoms.n15[i]):
           mscale[Atoms.i15[i,j]] = 1.
    Energi.em = em
    Energi.nem = nem