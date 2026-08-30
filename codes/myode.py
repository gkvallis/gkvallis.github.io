"""
This is a simple ode solver, for problems like
the Lorenz equations.  You can easily write your own functions for
different equation sets.
You can use Runge Kutta or Adams-Bashforth 3rd order.
However, the timestep is fixed --- there is no adaptive adjusting.

The code is not a black box; the user should it to their own needs!
It is not especially well documented, but it is simple enough just to follow
code itself. 

A version of this code was used in section 11.4 of the AOFD book. 

by G. K. Vallis
"""

# This version is set up for the Lorenz equations with RK4, or Adams Bashforth.
# Comment out the appropriate line in the timestepper.

from __future__ import (print_function, division)
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

def show_plot(figure_id=None):
    """
    Brings the plot to the foreground if using qt backend
    """
#    import matplotlib.pyplot as plt
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


def rk_4(tm, h, y, f):
    """
    Runge-Kutta integrator
    h is timestep
    tm is the time
    f is function to be integrated
    y are the variables
    """
    k1 = h * f(tm, y)
    k2 = h * f(tm + 0.5*h, y + 0.5*k1)
    k3 = h * f(tm + 0.5*h, y + 0.5*k2)
    k4 = h * f(tm + h, y + k3)
    return tm + h, y + (k1 + 2*(k2 + k3) + k4)/6.0
pass


def adbash3(tm, ic, dt, y, func, fm1, fm2):
    """ Adams Bashforth 3rd order
    Routine steps forward one time step
    """
    # f0 = np.zeros_like(y)
    if ic == 0:
        # forward euler
        dt1 = 1.
        dt2 = 0.0
        dt3 = 0.0
    elif ic == 1:
    #else:
        # AB2 at step 2
        dt1 = 1.5
        dt2 = -0.5
        dt3 = 0.0
    else:
        # AB3 from step 3 onwards
        dt1 = 23./12.
        dt2 = -16./12.
        dt3 = 5./12.
    f0 = func(tm, y)
    y = y + dt*(dt1*f0 + dt2*fm1 + dt3*fm2)
    fm2 = fm1
    fm1 = f0
    return tm+dt, y , fm1, fm2
pass



def Lorenz(t, state):
    """
    3 component Lorenz model, with parameters from a function
    """
    x, y, z = state
    sigma, beta, rho = get_params(state)
    dt_x = sigma*(y-x)
    dt_y = x*(rho-z)-y
    dt_z = x*y-beta*z
    return np.array([dt_x, dt_y, dt_z])

# @jit
def get_params(state):
    """
    returns parameters for ode function
    """
    params = np.array([10., 8./3., 28.])
    return params

time_start = time.clock()
#global gparams
t = 0
dt = 1./100.0
End_time = 40

N_steps = int(End_time/dt)
nplot = 1
N_plot = int(End_time/(dt*nplot))
print ('Nsteps, Nplot', N_steps, N_plot)


# This is the initial conditions for the system of ODES [x,y,z]
state = np.array([1.0, 1., 1.])
plot_array_x = np.zeros((N_plot))
plot_array_y = np.zeros((N_plot))
plot_array_z = np.zeros((N_plot))
time_array = np.zeros((N_plot))
fm1 = 0
fm2 = 0

print('%10f %10f %10f %10f' % (t, state[0], state[1], state[2]))

time_array[0] = t
plot_array_x[0] = state[0]
plot_array_y[0] = state[1]
plot_array_z[0] = state[2]


j = 0
nprint = 200

# following loop is the timestep integrator.
for i in range(N_steps):
    #t, state = rk4(t, dt, state, params, funLorenz)
    if (i % nplot == 0):
        time_array[j] = t
        plot_array_x[j] = state[0]
        plot_array_y[j] = state[1]
        plot_array_z[j] = state[2]
        j = j+1
        # print 'j =', j
    if (i % nprint == 0):
        print('%10f %10f %10f %10f' % (t, state[0], state[1], state[2]))
    pass
    # t, state = rk_4(t, dt, state, Lorenz)
    t, state, fm1, fm2 = adbash3(t, i, dt, state, Lorenz, fm1, fm2)
pass
pass
time_end = time.clock()
print ('CPU time = ', time_end - time_start)

plt.close('all')
plt.figure(1)
plt.plot(time_array, plot_array_x)
plt.xlabel('Time')
plt.ylabel('X')

plt.figure(2)
plt.subplot(1,2,1)
plt.plot(plot_array_z, plot_array_x,label="butterflies",color="#008DFF")
plt.legend() #(p1, ["butterfly"])
plt.savefig('Lorenz2d.pdf')
show_plot()

fig = plt.figure(3)
ax = fig.gca(projection='3d')
ax.plot(plot_array_x, plot_array_y, plot_array_z,color="#2294D2")
ax.set_xlabel("X", style='italic')
ax.set_ylabel("Y", style='italic')
ax.set_zlabel("Z", style='italic')
plt.savefig('Lorenz3d.pdf')
show_plot()
# show_plot(1)
# show_plot(3)
# plt.show()

