"""
This code calculates the Gill response analytically.
Author: G. K. Vallis
""" 

def myheaviside(x):
    return (x > 0)*1
myhv = myheaviside


def Hv(x):
    return np.heaviside(x, 0.5)


# Now the functions for the Gill problem

def QGill(x, y, A):
    return  A * np.cos(k*x) * Hv(Lf - abs(x))
    # return G
#
    

def Q0(x, y, A):
    q00 = -alk*(alpha*np.cos(k*x) + \
                        k * (np.sin(k*x) + np.exp(-alpha*(x + Lf))) )
    q00 = q00*Hv(Lf - np.abs(x))  # *Hv(x+Lf)
    q01 = -alk*k*(1 + np.exp(-2*alpha*Lf))  \
                                   * np.exp(alpha*(Lf - x) ) 
    q01 = q01 * Hv(x-Lf) 
    q0 = (q00 + q01)*A 
    return q0
#

def Q2(x, y, A):
    q21 = alk3*(-3*alpha*np.cos(k*x) + k * (np.sin(k*x) - \
                              np.exp(3*alpha*(x - Lf))))
    q21 = q21*Hv(Lf - abs(x))
    q22 = -alk3*k*(1 + np.exp(-6*alpha*Lf))*np.exp(3*alpha*(x+Lf))
    q22 = q22*Hv(-Lf -x)
    q2 = (q21 + q22)*A
    return q2
#

def Phi(x, y):
    """ This is the Gill solution for pressure """
    phi = 0.5 * (Q0(x, y, A) + Q2(x, y, A)*(y**2+1) ) * np.exp(-(y**2)/4)
    return phi
#

def uvel(x, y):
    u = 0.5*(Q0(x, y, A) + Q2(x, y, A)*(y**2 - 3))*np.exp(-y**2/4)
    return u

def vvel(x, y):
    v = y*(1.0*QGill(x, y, A) +  4.*alpha * Q2(x, y, A))*np.exp(-(y**2)/4)
    return v
#

def wvel(x, y):
    w = 0.5 * (2*QGill(x, y, A) + alpha*Q0(x, y, A) + \
             alpha * Q2(x, y, A) * (1 + y**2) ) * np.exp(-y**2/4) 
    return w
#



import numpy as np
import matplotlib as mpl
# from showplots import show_plot, draw_plot, show_plots
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
# from scipy import ndimage




# NUMERICAL and DOMAIN PARAMETERS
ifig = 0
Ly = 5.0     # Domain size y
Lx = 15.0    # Domain size x
A = 1        # amplitude of Gill forcing. Constant for now
ifig = 0
eps = np.spacing(1)  # small number
ny = 201          # Number of grid points in y
nx = 501          # Number of grid points in x
dy = Ly/(ny-1)    # Grid step in y (space)
dx = Lx/(nx-1)    # Grid step in x (space), different because periodic
rdy  = 1./(dy*2.)
rdx  = 1./(dx*2.)
rdy2 = 1./(dy**2)
rdx2 = 1./(dx**2)



# Parameters of Gill Problem:
Lf = 2.0     # forcing scale in Gill problem
alpha = 0.1  # Friction parameter 
# Convenience Physical Parameters
k = np.pi/(2*Lf)   # a wavenumber in Gill problem
alk = 1./(alpha**2 + k**2)
alk3 = 1./(9*alpha**2 + k**2) 


#  Initialize various fields 
myzeros = np.zeros((nx,ny))
v = np.zeros((nx,ny))
u = np.zeros((nx,ny))
w = np.zeros((nx,ny))
pressure = np.zeros((nx, ny)) 
gsol = np.zeros((nx, ny))
rsol = np.zeros((nx, ny))
ones = np.ones((nx, ny))

    
# Now set up the grid:
x = np.linspace(-Lx, +Lx, nx)
y = np.linspace(-Ly, Ly, ny)
# xx, yy = np.meshgrid(x, y)
xx, yy = np.meshgrid(x, y, indexing = 'ij')


# Nww find analytic solution
gsol[:, :] = 0.5 * Q0(xx, yy, A)*np.exp(-(y**2)/4)
rsol[:, :] = 0.5 * Q2(xx, yy, A)*(y**2+1)*np.exp(-(y**2)/4)
v[:, :] = vvel(xx, yy)
u[:, :] = uvel(xx, yy)
w[:, :] = wvel(xx, yy)
pressure[:, :] =  Phi(xx, yy)

# Now a couple of very simple plots. 
# You, dear reader, can easily add more, and prettier ones, very easily.
plt.close('all') 
ifig=1
pfig1 = plt.figure(ifig, figsize=(7, 7))
pfig1.canvas.manager.window.move(200, 100)
ax1 = pfig1.add_subplot(2,1,1)
plevels = np.linspace(-1.8,0,7)
ax1.contour(xx, yy, pressure, plevels) 
ax1.set_xlim(-10,15) 
ax2 = pfig1.add_subplot(2,1,2)
ax2.plot(x, pressure[:, ny//2]) 
ax2.set_xlim(-10,15) 
#
ifig += ifig
pfig2 = plt.figure(ifig, figsize=(7, 7))
pfig2.canvas.manager.window.move(900, 100)
ax1 = pfig2.add_subplot(2,1,1)
ax1.contour(xx,yy,gsol, plevels)  
ax1.set_xlim(-10,15) 
# ax1.plot(x, gsol[:, ny//2])  
ax2 = pfig2.add_subplot(2,1,2)
ax2.contour(xx,yy,rsol, plevels) 
ax2.set_xlim(-10,15) 
#ax2.plot(x, rsol[:, ny//2]) 
# show_plots(1, 2)  
plt.pause(1.e-5)
