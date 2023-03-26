#!/usr/bin/env python3
import argparse
import os
import mechanic
import field
import atoms
import getcart
import initial
import energi
import analysis

def enrgyze(Energi, Atoms, Field):
    analysis.analysis(Energi, Atoms, Field)
    # print out the total potential energy of the system
    prec = "{:16.4f}"
    space = 8
    if Field.digits >= 6:
        prec = "{:18.6f}"
        space = 6
    if Field.digits >= 8:
        prec = "{:20.8f}"
        space = 4
    fstr = "\n Total Potential Energy :" + " "*space + prec + " Kcal/mole"
    print(fstr.format(Energi.esum))

def partyze(Energi, Field):
    fstr = "\n Energy Component Breakdown :           Kcal/mole        Interactions\n"
    print(fstr)
    prec = "{:16.4f}"
    space = 5
    if Field.digits >= 6:
        prec = "{:18.6f}"
        space = 3
    if Field.digits >= 8:
        prec = "{:20.8f}"
        space = 1
    form1 = " "*space + prec + "{:17d}"
    if Field.use_mpole and (Energi.nem != 0 or Energi.em != 0.):
        fstr = " Atomic Multipoles" + " "*10 + form1
        print(fstr.format(Energi.em, Energi.nem))

if __name__ == "__main__":
    initial.initial()
    # Initialize the parser
    parser = argparse.ArgumentParser(
        description="Tinker Python"
    )

    # Add the parameters positional/optional
    parser.add_argument("xyzfile", help="Cartesian Coordinate File Name")
    parser.add_argument("analysis", help="Total Potential Energy and its Components [E]")
    parser.add_argument("-k", "--key")

    # Parse the arguments
    args = parser.parse_args()

    path = os.path.abspath(os.path.dirname(__file__))

    # Open xyz file
    if os.path.isfile(args.xyzfile):
        xyzfile = args.xyzfile
    else:
        xyzfile = f"{path}/{args.xyzfile}"

    # Open key file
    if args.key == None:
        fn = xyzfile.split(".")[0]
        keyfile = f"{fn}.key"
        if not os.path.isfile(keyfile):
            keyfile = "tinker.key"
    else:
        if os.path.isfile(args.key):
            keyfile = args.key
        else:
            keyfile = f"{path}/{args.key}"

    at = atoms.Atoms(xyzfile)
    ff = field.Field(keyfile)

    getcart.getcart(ff,at)
    mechanic.mechanic(ff, at)

    doenergy = False
    if args.analysis.upper() == "E":
        doenergy = True
    
    en = energi.Energi()

    if doenergy:
        enrgyze(en, at, ff)

    if doenergy:
        partyze(en,ff)
