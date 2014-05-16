import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import argparse

# Get program arguments
parser=argparse.ArgumentParser()
parser.add_argument('ncube',type=int)
parser.add_argument('ngrid',type=int)
parser.add_argument('anim_start',type=int)
parser.add_argument('timesteps',type=int)
parser.add_argument('anim_step',type=int)
parser.add_argument('path')
parser.add_argument('prefix')
args=parser.parse_args()

# Set arguments to parameter values
ncube=args.ncube
ngrid=args.ngrid
anim_start=args.anim_start
timesteps=args.timesteps
anim_step=args.anim_step
path='/mnt/scratch-lustre/'+args.path
prefix=args.prefix

# Do check
if (anim_start < 1):
    anim_start=1

################################### STUFF ######################################
def _blit_draw(self, artists, bg_cache):
    # Handles blitted drawing, which renders only the artists given instead
    # of the entire figure.
    updated_ax = []
    for a in artists:
        # If we haven't cached the background for this axes object, do
        # so now. This might not always be reliable, but it's an attempt
        # to automate the process.
        if a.axes not in bg_cache:
            # bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.bbox)
            # change here
            bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.figure.bbox)
        a.axes.draw_artist(a)
        updated_ax.append(a.axes)

    # After rendering all the needed artists, blit each axes individually.
    for ax in set(updated_ax):
        # and here
        # ax.figure.canvas.blit(ax.bbox)
        ax.figure.canvas.blit(ax.figure.bbox)

animation.Animation._blit_draw = _blit_draw
################################################################################


# START HERE
# Open figure and add subplots for each 2d projection
fig=plt.figure(num=1,figsize=(15,6))
ax1=fig.add_subplot(1,3,1,xlim=(0,ncube*ngrid),ylim=(0,ncube*ngrid),aspect='equal',xlabel='X',ylabel='Y')
img1, = ax1.plot([],[],'ko',markersize=2)
ax2=fig.add_subplot(1,3,2,xlim=(0,ncube*ngrid),ylim=(0,ncube*ngrid),aspect='equal',xlabel='X',ylabel='Z')
img2, = ax2.plot([],[],'ko',markersize=2)
ax3=fig.add_subplot(1,3,3,xlim=(0,ncube*ngrid),ylim=(0,ncube*ngrid),aspect='equal',xlabel='Y',ylabel='Z')
img3, = ax3.plot([],[],'ko',markersize=2)

# Setup a 'string' object to hold the changing title
time_text=ax2.text(.5,1.05,'',transform=ax2.transAxes,va='center')

# Initialize the blank axes and title
def init():
    time_text.set_text('')
    img1.set_data([0], [0])
    img2.set_data([0], [0])
    img3.set_data([0], [0])
    return time_text, img1, img2, img3

# Animate the title and scatter plot using the simulation's output data 
def animate(i,*args):
    time_text.set_text('dt='+str(i+1))
##    data=np.loadtxt(path+prefix+str(i+1)+'.txt')

    data=np.fromfile(path+prefix+str(i+1)+'.bin',dtype='float32')
    data=data.reshape(data.shape[0]/3,3)
##    img1.set_data(data[:,0],data[:,1]) # xy projection 
##    img2.set_data(data[:,0],data[:,2]) # xz projection
##    img3.set_data(data[:,1],data[:,2]) # yz projection

    index=np.arange(0,data.shape[0],100)

    img1.set_data(data[index,0],data[index,1])
    img2.set_data(data[index,0],data[index,2])
    img3.set_data(data[index,1],data[index,2])
    
    return time_text, img1, img2, img3

# Animate
anim = animation.FuncAnimation(fig,animate,init_func=init,frames=range(anim_start-1,timesteps,anim_step),interval=100,blit=True)

# Save animation as a mp4 (hopefully)000
#anim.save('animation.mp4', fps=30)

plt.show()
