import random
import numpy

for x in range (0, 100):
	P=random.random()
	s = numpy.random.binomial(1, P, 4)
	print s
