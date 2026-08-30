# Python qbo.py
# By G. K. Vallis
# Should work in both python 2 and python 3.

"""
1-d QBO, following the model of Plumb (1977)
Uses 3rd order Adam-Bashforth time differencing (a useful default method)

Momentum forcing from the lower boundary
by two gravity waves of opposite phase speed.

Grid goes from j=0 at z = 0 to j = nz-1 at z = ztop.
Boundary conditions are ub[0] = 0 and ub[nz-1] = ub[nz-2]
The pythonic way has the top at ub[nz-1]
A useful exercise would be to make a movie of the output.
"""


# Default parameters that work are amplitude = 1 and time = 50, or 0.1 and 500
from __future__ import (division, print_function)
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from numba import jit
import time as tm

# 
# This is just a function that brings the plot-figure to the top
# if you are using a qt backend. 
def show_plot(figure_id=None):
    import matplotlib.pyplot as plt
    if figure_id is not None:
        fig = plt.figure(num=figure_id)
    else:
        fig = plt.gcf()
    pass
    fig.canvas.draw()
    if 'qt' in plt.get_backend().lower():
        plt.show(block=False)
        plt.pause(1e-7)
        fig.canvas.manager.window.activateWindow()
        fig.canvas.manager.window.raise_()
    else:
        plt.show()
        plt.pause(1e-7)
    pass
#

# Here are a few cosmetic details
fn = 18
fontb = {'size': fn}
mpl.rc('font', **fontb)
mpl.rc('text', usetex = False)


# Now the real stuff:
# Define  constants and parameters
# A few physical (but nondimensional) values first
ztop = 4
famp = 1                                # forcing amplitude (arbitrary-ish
#   If you change famp you will need to change the plotting times.
ztop = 4                                # top boundary level
gwdamp = 1.                             # GW damping rate, usually unity
diffusv = 0.02                          # mean flow diffusion, e.g., 0.02
diffusv = 0.21*famp                     # diffusion scales with forcing amplit
k1 = 1 ; k2 = 1                         # zonal wavenumbers for the two waves
c1 = 1 ; c2 = -1                        # phase speeds for the two waves
#
rossby = 0                              # 0 for Plumb, 1 for Lindzen-Holton

# Now some numerical values
# note that z=0 and z=ztop correspond to j = 0 and j = nz-1
nz = 100                                # number of vertical gridpoints
nzm1 = nz - 1                           # 100 is good number to play with
dz = 1.*ztop/nz                         # grid size
dz2 = dz**2
rdz2 = 1./dz2
dt = 0.001/famp                         # Scale timestep with amplitude
# time_stop = 40                        # stopping time
time_stop = 41/famp                     # Maybe make it amplitude dependent
# as a rule of thumb, 20000 steps takes a few seconds on a laptop, with jit
# showing the movie of ub in realtime slows it down a lot
nsteps = int(time_stop/dt)
# nsteps = 20000
# nplot1 is used to make a simple movie. 
# Set nplot1 = 200 for a movie (slow), set nplot1 = 200000 for no movie (fast)
nplot1 = 500                          # Plot u-movie evey nplot1 (e.g 100)steps
nplot2 = 50                           # plot contour every so many time steps
nplot2len = np.int(nsteps/nplot2)     # size of plot array
time = 0                              # initial time


# Now initialize some of the basic arrays for later use
z = np.linspace(0, ztop, nz)            # vertical prediction levels
zplot = np.linspace(0, ztop, nz)        # plotting levels in vertical
ubplot = np.zeros((nz, nplot2len))      # time series for mean wind
timeplot = np.zeros(nplot2len)
G1 = np.zeros((nz, nplot2))             # time series for wave 1 flux
G2 = np.zeros((nz, nplot2))             # time series for wave 2 flux
F_gw1z = np.zeros(nz)                   # vector for wave 1 forcing at time n
F_gw2z = np.zeros(nz)                   # vector for wave 2 forcing at time n
f0 = np.zeros(nz)                       # vector for total forcing at time n
fm1 = np.zeros(nz)                      # forcing at n-1
fm2 = np.zeros(nz)                      # at time n-1

ub0 = 0                                 # mean wind at z = 0
ub = 0.1*np.sin(np.pi*z/4)              # initial mean wind profile
# ub = z*0.1                            # or try a shear


############################################################################


                                     # the jit decorator below is optional
                                     # numba/jit is a just-in-time LLVM compiler.
@jit                                 # It provides a big speedup (> x20)
def rhs_mean(Am, gwdamp, k1, k2, c1, c2, ub, f0, F_gw1z, F_gw2z,
             diffusv, dz, dz2, nz):
    """ rhs of mean wind equation """
    # first compute vertical dependence of momentum flux divergence
    # compute the eastward and westward wave forcing
    if rossby == 1: gwdamp = 1.
    gw1z = gwdamp/(k1*(ub-c1)**2)
    gw2z = gwdamp/(k2*(ub-c2)**2)
    gw10 = gwdamp/(k1*(ub0-c1)**2)
    gw20 = gwdamp/(k2*(ub0-c2)**2)
#

# Now, if rossby = 1, replace G2 with the expression for Rossby waves
    if rossby == 1:
        beta = 32;
        amplit = 1;
        k2 = 3 ;
        gw20 = amplit*gwdamp/(k2*(ub0-c2)**2)*(1/(k2*(ub0 - c2)) - 1);
        gpart = amplit*gwdamp/(k2*(ub-c2)**2) ;
        rpart = beta/(k2**2*(ub - c2)) - 1 ;
        gw2z = gpart*rpart ;
    # now evaluate effect of gwave forcing
    for j in range(0, nz):
        F_gw1z[j] = famp*gw1z[j]*np.exp(-(0.5*(gw10+gw1z[j]) +
                                  np.sum(gw1z[1:j-1]))*dz)
        F_gw2z[j] = -famp*gw2z[j]*np.exp(-(0.5*(gw20+gw2z[j]) +
                                  np.sum(gw2z[1:j-1]))*dz)
    pass
    # now add wave plus mean flow diffusion terms to  get total forcing
    f0[1:nzm1] = F_gw1z[1:nzm1] + F_gw2z[1:nzm1]  \
        + diffusv*rdz2*(ub[2:nz] - 2.*ub[1:nz-1]+ub[0:nz-2])
    return (F_gw1z, F_gw2z, f0)
pass


def ustep(ic, ub, time, f0, fm1, fm2):
    """time stepping for mean wind equation
       advance ubar by one timestep, 3rd order AB
       The right-hand side is input as array f0
    """
    if ic == 0:
        # forward euler
        dt1 = dt
        dt2 = 0.0
        dt3 = 0.0
    elif ic == 1:
        # AB2 at step 2
        dt1 = 1.5*dt
        dt2 = -0.5*dt
        dt3 = 0.0
    else:
        # AB3 from step 3 on
        dt1 = 23./12.*dt
        dt2 = -16./12.*dt
        dt3 = 5./12.*dt
    ub[1:nz-1] = ub[1:nz-1] + \
        (dt1*f0[1:nz-1] + dt2*fm1[1:nz-1] + dt3*fm2[1:nz-1])
    ub[nz-1] = ub[nz-2]      # top boundary condition of no gradient
    # ub[0] = 0  # not needed as it is not evolved.
    fm2[:] = fm1[:]
    fm1[:] = f0[:]
    time = time + dt
    return ub, fm1, fm2, time
pass

t_start = tm.clock()
ifig = 1
count = 0
# Begin time stepping
iplot1 = 0
iplot2 = 0
print ('nsteps, nplot1, nplot2', nsteps, nplot1, nplot2)
for ic in range(0,nsteps):
    count = count + 1
    F_gw1z, F_gw2z, f0 = rhs_mean(famp, gwdamp, k1, k2, c1, c2, ub, f0,
                                  F_gw1z, F_gw2z, diffusv, dz, dz2, nz)
    ub, fm1, fm2, time = ustep(ic, ub, time, f0, fm1, fm2)
    # now a movie
    if  ic % nplot1 == 0:
        iplot1 += 1
        fig1 = plt.figure(ifig)  # make a movie of wind
        if iplot1 == 1:
            show_plot()      # this brings the window to the front
        plt.clf()
        ax1 = fig1.add_subplot(1, 1, 1)
        ax1.plot(ub, zplot)
        plt.xlim((-1.0, +1.0))
        ax1.set_xlabel('zonal wind')
        ax1.set_ylabel('Height')
        plotlabel = "Time = %1.2f" % time
        ax1.set_title(plotlabel)
        plt.draw()
        plt.pause(1e-7)
    if  (ic % nplot2 == 0):
        ubplot[:,iplot2] = ub[:]
        timeplot[iplot2] = time
        iplot2 += 1
pass
t_end = tm.clock()
print ('CPU time = ', t_end - t_start)

# Now do some plotting. Some of the parameters here are hardwired for a 
# forcing amplitude of one. 
ifig += 1
fig2 = plt.figure(ifig)
plt.clf()
iplstart = 0
plt.set_cmap('Spectral_r')
ax1 = fig2.add_axes([0.1, 0.1, 0.8, 0.5])
ax1.set_xlabel ('Time',size='18')
ax1.set_ylabel ('Height', size='18')
ax1.set_title('Zonal mean zonal wind',size=18)
pcc = ax1.contourf(timeplot[iplstart:], z, ubplot[:,iplstart:])
ax1.contour(timeplot[iplstart:], z, ubplot[:,iplstart:], colors='k',linestyles='solid')
ax1.contour(timeplot[iplstart:], z, ubplot[:,iplstart:], levels=[0], linewidths=2)
ax2 = fig2.add_axes([0.92, 0.1, 0.02, 0.5])
ax1.set_xticks([20, 25, 30, 35, 40])
ax1.set_xlim((19.5,40.5))
ax1.set_yticks([0, 1, 2, 3., 4])
ax1.set_ylim((0, 3.2))
fig2.colorbar(pcc, cax=ax2)
plt.savefig("qbo-contour.pdf")
show_plot()

ifig += 1
fig3 = plt.figure(ifig)
plt.clf()
plint = 35
for ipl in range(1,7):
    itpl = 400+(ipl-1)*plint
    timepl = timeplot[itpl]
    ax = fig3.add_subplot(2, 3, ipl)
    ax.plot(ubplot[:,itpl], zplot[:],linewidth=2)
    ax.set_xticks([-1, 0, 1.])
    ax.set_xlim((-1.0, +1.0))
    # ax.set_xlim((-0.8, +0.8))
    ax.set_yticks([0, 1, 2, 3.])
    ax.set_ylim((0, +3.1))
    plotlabel = "t = %1.2f"  % timepl
    ax.set_title(plotlabel,size=12)
    if ipl == 1:
        ax.set_ylabel('Height', size=14)
    elif ipl ==4:
        ax.set_ylabel('Height', size=14)
    if ipl > 3:
        ax.set_xlabel('u', style='italic', size=14)
    pass
pass
plt.savefig("qbo-period.pdf")  # saves a PDF plot
show_plot()

ifig += 1
fig = plt.figure(ifig)
plt.clf()
ax = fig.add_subplot(1, 1, 1)
ax.plot(timeplot, ubplot[50,:], 'k--', label = 'z = 2', linewidth=2)
ax.plot(timeplot, ubplot[13,:], 'b-', label = 'z = 0.5', linewidth=2)
ax.set_yticks([-0.8, -0.4, 0, 0.4, 0.8 ])
ax.set_ylim((-0.8, 0.8))
ax.set_xticks(np.arange(0, 41, 10))
ax.set_xlim((0, 40))
ax.tick_params(labelsize=18)
ax.set_xlabel ('Time',size='20')
ax.set_ylabel ('Zonal wind', size='20')
ax.legend(frameon = True, prop={'size':20,'style':'italic'},loc=3)
show_plot()
plt.show()
plt.savefig("qbowind-time.pdf")  # saves a PDF plot
