#!/home/delta137/anaconda/bin/python
# Methode de shooting, pour obtenir la solution
# de kink, voir Rajaraman chapitre 2 pour les eq diff
# (2.25) et la solution analytique (2.28)

import numpy
import math
import pylab

x = numpy.linspace(1, 10, 100)
h = x[1] - x[0]  # stepsize
y1=x/10
