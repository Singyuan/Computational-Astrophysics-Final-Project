import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# read the text file
firstcol = np.loadtxt("anidata.txt", usecols=(0,), dtype=int) # read first column
N_In = firstcol[0] # x axis not include ghost
U = np.loadtxt("anidata.txt",skiprows=1) # read anthoer row
totalrow = len(U) # count the total length of data
frame = int(totalrow/N_In)

# x axis
x = U[0:N_In, 1]

# initialize animation
def init():
   line_d.set_xdata( x )
   line_u.set_xdata( x )
   line_p.set_xdata( x )
   return line_d, line_u, line_p

# update animation
def update( i ):
   oldline = i*N_In
   newline = (i+1)*N_In
   if ( newline> totalrow+1):
      anim.event_source.stop()

   t = U[oldline, 0]
   d = U[oldline:newline, 2]
   u = U[oldline:newline, 3]
   P = U[oldline:newline, 4]

   line_d.set_ydata( d )
   line_u.set_ydata( u )
   line_p.set_ydata( P )
   ax[0].set_title( 't = %3.4f' % (t) )

   return line_d, line_u, line_p


fig, ax = plt.subplots( 3, 1, sharex=True, sharey=False, dpi=140 )
fig.subplots_adjust( hspace=0.1, wspace=0.0 )
#fig.set_size_inches( 6.4, 12.8 )
line_d, = ax[0].plot( [], [], 'r-o', ls='-', markeredgecolor='k', markersize=3 )
line_u, = ax[1].plot( [], [], 'b-o', ls='-', markeredgecolor='k', markersize=3 )
line_p, = ax[2].plot( [], [], 'g-o', ls='-', markeredgecolor='k', markersize=3 )


anim = animation.FuncAnimation( fig, func=update, init_func=init, frames=frame, interval=10, repeat=False )
ax[2].set_xlabel( 'x' )
ax[0].set_ylabel( 'Density' )
ax[1].set_ylabel( 'Velocity' )
ax[2].set_ylabel( 'Pressure' )

ax[0].set_xlim( 0.0, 1.0 )
ax[0].set_ylim( +0.0, 1.2 )
ax[1].set_ylim( +0.0, 1.5 )
ax[2].set_ylim( +0.0, 1.2 )
plt.show()