import os
import sys

path = os.path.abspath(os.path.dirname(__file__))
parent = os.path.abspath(os.path.join(path, os.pardir))
source = os.path.join(parent,'source')
test = os.path.join(parent,'test')
sys.path.append(source)
sys.path.append(test)
