import atoms
import numpy as np
import math
from mathconstants import *

def rotpole(Atoms, poltyp):
    # rotate local multipoles to global frame at each site
    a = np.identity((3))
    if poltyp.upper() == "MPOLE":
        for i in range(Atoms.n):
            if Atoms.pollist[i] != -1:
                planar = rotmat(Atoms,i,a)
                rotsite(Atoms,i,a,planar)

def rotmat(Atoms,i,a):
    # get coordinates and frame definition for multipole site
    xi = Atoms.x[i]
    yi = Atoms.y[i]
    zi = Atoms.z[i]
    iz = Atoms.zaxis[i]
    ix = Atoms.xaxis[i]
    iy = abs(Atoms.yaxis[i])
    axetyp = Atoms.polaxe[i]
    planar = False
    # get Z-Only rotation matrix elements for z-axis only
    if axetyp == "Z-Only":
        dx = Atoms.x[iz] - xi
        dy = Atoms.y[iz] - yi
        dz = Atoms.z[iz] - zi
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        a[0,2] = dx / r
        a[1,2] = dy / r
        a[2,2] = dz / r
        dx = 1.
        dy = 0.
        dz = 0.
        dot = a[0,2]
        eps = 0.707
        if abs(dot) > eps:
            dx = 0.
            dy = 1.
            dot = a[1,2]
        dx = dx - dot*a[0,2]
        dy = dy - dot*a[1,2]
        dz = dz - dot*a[2,2]
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        a[0,0] = dx / r
        a[1,0] = dy / r
        a[2,0] = dz / r
    # get Z-then-X rotation matrix elements for z- and x-axes
    elif axetyp == "Z-then-X":
        dx = Atoms.x[iz] - xi
        dy = Atoms.y[iz] - yi
        dz = Atoms.z[iz] - zi
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        a[0,2] = dx / r
        a[1,2] = dy / r
        a[2,2] = dz / r
        dx = Atoms.x[ix] - xi
        dy = Atoms.y[ix] - yi
        dz = Atoms.z[ix] - zi
        dot = dx*a[0,2] + dy*a[1,2] + dz*a[2,2]
        dx = dx - dot*a[0,2]
        dy = dy - dot*a[1,2]
        dz = dz - dot*a[2,2]
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        a[0,0] = dx / r
        a[1,0] = dy / r
        a[2,0] = dz / r
    # get Bisector rotation matrix elements for z- and x-axes
    elif axetyp == "Bisector":
        dx = Atoms.x[iz] - xi
        dy = Atoms.y[iz] - yi
        dz = Atoms.z[iz] - zi
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        dx1 = dx / r
        dy1 = dy / r
        dz1 = dz / r
        dx = Atoms.x[ix] - xi
        dy = Atoms.y[ix] - yi
        dz = Atoms.z[ix] - zi
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        dx2 = dx / r
        dy2 = dy / r
        dz2 = dz / r
        dx = dx1 + dx2
        dy = dy1 + dy2
        dz = dz1 + dz2
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        a[0,2] = dx / r
        a[1,2] = dy / r
        a[2,2] = dz / r
        dot = dx2*a[0,2] + dy2*a[1,2] + dz2*a[2,2]
        dx = dx2 - dot*a[0,2]
        dy = dy2 - dot*a[1,2]
        dz = dz2 - dot*a[2,2]
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        a[0,0] = dx / r
        a[1,0] = dy / r
        a[2,0] = dz / r
    # get Z-Bisect rotation matrix elements for z- and x-axes;
    # use alternate x-axis if central atom is close to planar
    elif axetyp == "Z-Bisect":
        dx = Atoms.x[iz] - xi
        dy = Atoms.y[iz] - yi
        dz = Atoms.z[iz] - zi
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        a[0,2] = dx / r
        a[1,2] = dy / r
        a[2,2] = dz / r
        dx = Atoms.x[ix] - xi
        dy = Atoms.y[ix] - yi
        dz = Atoms.z[ix] - zi
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        dx1 = dx / r
        dy1 = dy / r
        dz1 = dz / r
        dx = Atoms.x[iy] - xi
        dy = Atoms.y[iy] - yi
        dz = Atoms.z[iy] - zi
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        dx2 = dx / r
        dy2 = dy / r
        dz2 = dz / r
        dx = dx1 + dx2
        dy = dy1 + dy2
        dz = dz1 + dz2
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        dx = dx / r
        dy = dy / r
        dz = dz / r
        dot = dx*a[0,2] + dy*a[1,2] + dz*a[2,2]
        angle = 180. - radian*math.acos(dot)
        # eps = 15.
        eps = 0.
        if angle < eps:
            planar = True
            dx = dy1*dz2 - dz1*dy2
            dy = dz1*dx2 - dx1*dz2
            dz = dx1*dy2 - dy1*dx2
            dot = dx*a[0,2] + dy*a[1,2] + dz*a[2,2]
            if dot < 0.:
                dx = -dx
                dy = -dy
                dz = -dz
                dot = -dot
        dx = dx - dot*a[0,2]
        dy = dy - dot*a[1,2]
        dz = dz - dot*a[2,2]
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        a[0,0] = dx / r
        a[1,0] = dy / r
        a[2,0] = dz / r
    # get 3-Fold rotation matrix elements for z- and x-axes;
    # use alternate z-axis if central atom is close to planar
    elif axetyp == "3-Fold":
        dx = Atoms.x[iz] - xi
        dy = Atoms.y[iz] - yi
        dz = Atoms.z[iz] - zi
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        dx1 = dx / r
        dy1 = dy / r
        dz1 = dz / r
        dx = Atoms.x[ix] - xi
        dy = Atoms.y[ix] - yi
        dz = Atoms.z[ix] - zi
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        dx2 = dx / r
        dy2 = dy / r
        dz2 = dz / r
        dx = Atoms.x[iy] - xi
        dy = Atoms.y[iy] - yi
        dz = Atoms.z[iy] - zi
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        dx3 = dx / r
        dy3 = dy / r
        dz3 = dz / r
        dx = dx1 + dx2 + dx3
        dy = dy1 + dy2 + dy3
        dz = dz1 + dz2 + dz3
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        # eps = 0.15
        eps = 0.
        if r < eps:
            planar = True
            dx2 = Atoms.x[ix] - Atoms.x[iz]
            dy2 = Atoms.y[ix] - Atoms.y[iz]
            dz2 = Atoms.z[ix] - Atoms.z[iz]
            dx3 = Atoms.x[iy] - Atoms.x[iz]
            dy3 = Atoms.y[iy] - Atoms.y[iz]
            dz3 = Atoms.z[iy] - Atoms.z[iz]
            dx4 = dy2*dz3 - dz2*dy3
            dy4 = dz2*dx3 - dx2*dz3
            dz4 = dx2*dy3 - dy2*dx3
            dot = dx4*dx + dy4*dy + dz4*dz
            if dot > 0.:
                dx = dx4
                dy = dy4
                dz = dz4
            else:
                dx = -dx4
                dy = -dy4
                dz = -dz4
            r = math.sqrt(dx*dx + dy*dy + dz*dz)
        a[0,2] = dx / r
        a[1,2] = dy / r
        a[2,2] = dz / r
        dot = dx1*a[0,2] + dy1*a[1,2] + dz1*a[2,2]
        dx = dx1 - dot*a[0,2]
        dy = dy1 - dot*a[1,2]
        dz = dz1 - dot*a[2,2]
        r = np.sqrt(dx*dx + dy*dy + dz*dz)
        a[0,0] = dx / r
        a[1,0] = dy / r
        a[2,0] = dz / r
    # finally, find rotation matrix elements for the y-axis
    a[0,1] = a[2,0]*a[1,2] - a[1,0]*a[2,2]
    a[1,1] = a[0,0]*a[2,2] - a[2,0]*a[0,2]
    a[2,1] = a[1,0]*a[0,2] - a[0,0]*a[1,2]
    return planar

def rotsite(Atoms,ii,a,planar):
    inpole = Atoms.pole[ii]
    outpole = Atoms.rpole[ii]
    # copy input multipoles and modify at planar sites
    spole = np.copy(inpole)
    if planar:
        axetyp = Atoms.polaxe[ii]
        if axetyp == "Z-Bisect":
            spole[1] = 0.
            spole[6] = 0.
            spole[10] = 0.
            spole[4] = 0.5 * (spole[4]+spole[8])
            spole[8] = spole[4]
        elif axetyp == "3-Fold":
            spole[1:] = 0
    # monopoles are the same in any coordinate frame
    outpole[0] = spole[0]
    # rotate input dipoles to final coordinate frame
    for i in range(1,4):
        outpole[i] = 0.
        for j in range(1,4):
            outpole[i] += spole[j]*a[i-1,j-1]
    # rotate input quadrupoles to final coordinate frame
    mp = np.zeros((3,3))
    rp = np.zeros((3,3))
    k = 4
    for i in range(3):
        for j in range(3):
            mp[i,j] = spole[k]
            rp[i,j] = 0.
            k += 1
    for i in range(3):
        for j in range(3):
            if j < i:
                rp[i,j] = rp[j,i]
            else:
                for k in range(3):
                    for m in range(3):
                        rp[i,j] += a[i,k]*a[j,m]*mp[k,m]
    k = 4
    for i in range(3):
        for j in range(3):
            outpole[k] = rp[i,j]
            k += 1
    return
