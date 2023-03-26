import numpy as np

def readClean(file):
    with open(file) as f:
        lines = f.read().splitlines()

    lines = [i for i in lines if i]
    lines = [i.strip() for i in lines]
    lines = [i for i in lines if i[0] != '#']

    return lines

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
