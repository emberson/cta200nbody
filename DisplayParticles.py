import numpy as np
import matplotlib.pyplot as pl
import os
import time

def display(ngrid, ncube, path):
	delay = 0.001
	files = os.listdir(path)
	numFiles = len(files)
	pl.ion()
	pl.figure(1)

	for i in xrange(1, numFiles + 1):
		x, y, z, vx, vy, vz = np.genfromtxt(path + "TimeStamp" + str(i) + ".txt", dtype = float, unpack = True)
		pl.scatter(x, y, s = 3)
		pl.axes().set_xlim((0., float(ngrid * ncube)))
		pl.axes().set_ylim((0., float(ngrid * ncube)))
		pl.xlabel("X Position")
		pl.ylabel("Y Position")
		pl.title("Jd Is Our Favourite Best")
		pl.draw()
		time.sleep(delay)
		pl.clf()

	pl.close(1)

