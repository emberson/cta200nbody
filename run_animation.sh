# Arguments: ncube, ngrid, first animation frame number, last animation frame
#            number, animation frame step, path to data, prefix to data files
#
# *if animation starts at 1 and animation step is 1, then set last animation 
#  frame to 'timesteps' to animation all data

##python animationmaster.py 2 4 0 200 1 njones/cta200nbody/TimeStamps/ TimeStamp
python animationmaster.py 90 1 100 529 5 dpereira/HaloPos/ HaloPos
