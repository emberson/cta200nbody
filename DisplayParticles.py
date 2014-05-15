import numpy as np
from matplotlib import pyplot as pl
#import matplotlib.pyplot as pl
import os
import time

def display(ncube, ngrid):
	# Time delay before displaying new plot.
	delay = 0.001

	path = 'TimeStamps/'

	# Determines the number of files over which to loop.
	files = os.listdir(path)
	numFiles = len(files)

	# Prepares the plot to be displayed.
	#pl.ion()
	#pl.figure(1)

	# Loops over all files and extracts the component-based data.
	for i in xrange(1, numFiles + 1):
		x, y, z, vx, vy, vz, t  = np.genfromtxt(path + "TimeStamp" + str(i) + ".txt", dtype = float, unpack = True)
		pl.figure(1)
		pl.clf()
		# Plots x positions against y positions to get an xy-plane slice. 
		pl.scatter(x, y, s = 3)
		pl.axes().set_xlim((0., float(ngrid * ncube)))
		pl.axes().set_ylim((0., float(ngrid * ncube)))
		pl.xlabel("X Position")
		pl.ylabel("Y Position")
		pl.title("JD Is Our Favourite Best")
		pl.savefig(path + "TimeStamp" + str(i) + ".png")

		# Draws the plot and then delays the next iteration.
		#pl.draw()
		#time.sleep(delay)
		

	# Closes the plot at the very end.
	pl.close(1)

display(4,2)
