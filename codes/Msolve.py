
"""
 A solver for the one-dimensional M-equation in thermocline dynamics.
 Uses scipy.integrate_bvp

 This code can be used to reproduce figs 20.12 and 20.13 in Vallis (2017)
 The code solves the boundary value problem:
     w w_{zzz}  = \kappa w_{zzzz}  in 0 > z > -1
 with bounday conditions on w and w_{zz} at z = 0, -1.
 To understand this code you need two things:
     1. To understand tha equation itself and the boundary condtions.
        This is described in  Vallis (2017), the section 20.6 on the
        internal thermocline.
     2. To understand the functionality of the scipy function solve_bvp.

 For the record, the actual figures in Vallis were obtained
 with a Fortran code that differenced and then solved the nonlinear BV
 equations 'by hand'. The code here is shorter and simpler, but uses a
 scipy library function and so is more of a black box.

 Author: G. K. Vallis
"""

from __future__ import (division, print_function)
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
# from mystuff import *
from scipy.integrate import solve_bvp
import time as tm

# Here are a few cosmetic details:
fn = 16
fontb = {'size': fn}
mpl.rc('font', **fontb)
mpl.rc('text', usetex = False)

# Now set the physical paramters
kappa = 4e-4    # vertical diffusivity = 4e-4 in book
# kappa = 1e-6
kr = 1./kappa
# kr = 10000
# Following are the boundary conditions.
# t for temperature, w for vertical celocity
ttop =  10.0 ; tbot = 0
wtop = -1; wbot =  0


def fun(x, y):
    """
    This is the routine that defines the equation to be solved.
    We solve w_{zzzz} = kr w w_{zzz}
    """
    return np.vstack((y[1], y[2], y[3], kr*y[0]*y[3] ))
pass

def bc(ya, yb):
    """
    THis is the routine for the boundary conditions
    the values ttop, tbot, wtop, wbot are given.
    """
    return np.array([ya[0]-wbot, yb[0] - wtop, ya[2] - tbot, yb[2] + ttop])
pass

# Now define the initial grid.
n_init = 200    #  may need to go much higher if kappa is very small
x = np.linspace(-1, 0, n_init)
y_a = np.zeros((4, x.size))

# now find the solution
res_a = solve_bvp(fun, bc, x, y_a)

# Now extract the solution for plotting
n_plot = 200   # 200 is enough for most plots
x_plot = np.linspace(-1, 0, n_plot)
y_plot_a = res_a.sol(x_plot)[0]
y_plot_b = res_a.sol(x_plot)[2]
y_plot_c = res_a.sol(x_plot)[3]

# Now the plotting.
# First get a quick look at everything
plt.figure(1)
plt.clf()
# Note arbitrary factors in these plots
plt.plot(5*y_plot_a, x_plot, label = '5 * Vertical velocity')
plt.plot(-y_plot_b, x_plot, label='Temperature')
plt.plot(-y_plot_c*0.1, x_plot, label='0.1 * dT/dz')
plt.legend(loc='best')
plt.xlabel("Temperature and Vertical velocity")
plt.ylabel("Depth, z")
plt.pause(1.e-11)

# Now make a figure similar to Fig. 20.12 in Vallis (2017)
# If you change parameters you may need to change the plot axis limits.
plt.figure(2)
plt.clf()
plt.subplot(1,3,1)
plt.plot(y_plot_a, x_plot)
plt.title('Vertical velocity')
plt.ylabel("Depth, z")
plt.xticks(np.arange(-1.5,1.1,0.5))
plt.xlim(-1.2,0.2)
plt.subplot(1,3,2)
plt.plot(-y_plot_b, x_plot)
plt.title('Temperature')
plt.xticks(np.arange(-5,20,5))
plt.yticks([])
plt.xlim(-1,12)
plt.subplot(1,3,3)
plt.plot(-y_plot_c, x_plot)
xmax = np.max(-y_plot_c)
plt.xlim(-0.2*xmax,1.2*xmax)
plt.locator_params(axis='x', tight = True, nbins=6)
plt.yticks([])
plt.title('dT/dz')
plt.pause(1.e-11)
plt.savefig('int_therm.pdf')
plt.show()
