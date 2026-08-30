#  Advection-condensation-diffusion model. 
# GK Vallis, 2015.

"""
This is perhaps the simplest example of a 2D advection-diffusion-condensation model, 
as described in chapter 18 of AOFD.
In the example below, a single clockwise streamfunction advects water vapour. 
The water vapour souce is the bottom surface, which  is saturated.
The temperature falls linearly from 26 C to - 50C at the top.
The saturated mixing ratio is calculated using a Clausius-Clapeyron relation.
The water vpour condenses when it becomes saturated, in a 'fast condensation' limit.
The advection is a simple, flux-form, scheme. In this code the streamfunction is assumed constant
in time.
"""


def diffuse_wv(lapphi,phi,Ky,Kz,rdy2,rdz2):
   """
   This function
   evaluates the derivatives in the Laplacian.
   """
   lapphi[1:-1, 1:-1] = \
        Ky*(phi[2:, 1:-1] - 2*phi[1:-1, 1:-1] + phi[:-2,1:-1])*rdz2 + \
        Kz*(phi[1:-1, 2:] - 2*phi[1:-1, 1:-1] + phi[1:-1,:-2])*rdy2
   return lapphi
#


def advect_wv(advphi,v,w,phi,rdy,rdz):
    """
    Advects the variable phi by v and w
    This provides the positive advective derivative in flux form
    """
    advphi[1:-1, 1:-1] = \
     (phi[2:, 1:-1]*w[2:, 1:-1] - phi[:-2, 1:-1]*w[:-2, 1:-1])*rdz + \
     (phi[1:-1, 2:]*v[1:-1, 2:] - phi[1:-1, :-2]*v[1:-1, :-2])*rdy
    return advphi
#


def stream(v, w, psi):
    """
    Calculates velocities from stream function
    convention is v = - \pp psi z,  w = \pp psi y
    """
    w[:, 1:-1] = (psi[:, 2:] - psi[:, 0:-2]) * rdy
    v[1:-1, :] = - (psi[2:, :] - psi[0:-2, :]) * rdz
    return v, w
#


def CC(esat, T):
    """ Clausis Clapeyron Tetens style """
    esat[:, :] = 0.3619*np.exp(17.67*T[:, :]/(T[:, :]+243.3))
    return esat
#


def condense(Sink,wvsat,wv,myzeros,rtau):
    """ Provides a condensation sink for water vapour """
    Sink[:, :] = (wvsat[:, :] - wv[:, :])*rtau
    Sink[:, :] = np.minimum(Sink, myzeros)
    return Sink
#



def show_plot(figure_id=None):
    """
    Draws the figure and raises it to the foreground
    Does not raise the figure if not a qt backend.
    """
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



import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# from mystuff import *

mpl.rcParams.update({'figure.autolayout': True})
#font1 = {'family': 'Myriad Pro',
#         'weight': 'normal',
#         'style': 'normal',
#         'size': 14}

# plt.close('all')
# plt.close(1)
# plt.clf()



#  PHYSICAL PARAMETERS
Ky = 0.005    # Diffusion coefficient in y direction
Kz = 0.005    # Diffusion coefficient in z direction
Ly = 1.0     # Domain size y
Lz = 1.0     # Domain size z
Ltime = 1.   # Integration time  (may be reset later - see nt)
So = 1.0     # Source term

#  NUMERICAL PARAMETERS

ifig = 0      #  figure number
ny = 100       # Number of grid points in x (increase for publication)
nz = 102       # Number of grid points in y
dt = 0.0003      # Grid step (time)
dy = Ly/(ny-1)   # Grid step in y (space)
dz = Ly/(nz-1)   # Grid step in z (space)
nt = Ltime/dt    #  number of timesteps set physically
nt = 15000       #  nmber of timesteps set numerically
rdy  = 1./(dy*2.)
rdz  = 1./(dz*2.)
rdy2 = 1./(dy**2)
rdz2 = 1./(dz**2)
tau = 2*dt
rtau = 1./tau
eps = np.spacing(1)  #  a very small number

yy = np.linspace(0,Ly,ny)   #  y grid
zz = np.linspace(0,Lz,nz)   #  z grid
yv, zv = np.meshgrid(yy, zz)

#  Initialize fields to some value
myzeros = np.zeros((nz,ny))
wv_new = np.zeros((nz,ny))
wv_old = np.zeros((nz,ny))
S = np.zeros((nz,ny))
rh = np.zeros((nz,ny))
T = np.zeros((nz,ny))
wvsat = np.zeros((nz,ny))
lapwv= np.zeros((nz,ny))
wv = np.zeros((nz,ny))
v = np.zeros((nz,ny))
w = np.zeros((nz,ny))
diffwv = np.zeros((nz,ny))
advectwv = np.zeros((nz,ny))
est = np.zeros((nz,ny))
psi = np.zeros((nz,ny))    #  streamfunction

#  Set up a non-zero streamfunction
# psi[:,:] = 2*np.sin(yv[:,:]*3*np.pi)*np.sin(zv[:,:]*np.pi) \
#            * np.exp(-6.*yv[:,:]) *np.exp(-3*zv[:,:])\
#            + 0.0*np.sin(yv[:,:]*1 *np.pi)*np.sin(zv[:,:]*np.pi) \
#            + 1.0*np.sin(yv[:,:]*1 *np.pi)*np.sin(zv[:,:]*np.pi) \
#            * np.exp(-1.*yv[:,:])*np.exp(-3*zv[:,:])

# This is a single cell. 
psi[:,:] = np.sin(yv[:,:]*np.pi)*np.sin(zv[:,:]*np.pi)

(v,w) = stream(v,w,psi)    #  get velocities from streamfunction

#  Now set temperature, in Celsius
#  Set equatorial and polar temperatures, and interpolate linearly
tequator = 26.
tpole = 26
ttrop  = -50
#  tpole = 26
T[0,0] = tequator
T[-1,0] = ttrop
T[0,-1] = tpole
T[-1,-1] = ttrop
#  Now interpolate horizontally along ground and tropopause
T[0,:] = (T[0,-1] - T[0,0])*yy[:] + T[0,0]
T[-1,:] = (T[-1,-1] - T[-1,0])*yy[:] + T[-1,0]
#  Now fill in by vertical interpolation

for j in range(0,nz):
  T[j,:] = (T[-1,:] - T[0,:])*zz[j] + T[0,:]
#

wvsat = CC(est,T)

plt.set_cmap('Blues')
if 1:    #  plot or not the basic set up
    ifig += 1
    plt.clf()
    plt.figure(ifig)
    plt.subplot(2,2,1)
    plt.title('streamfunction')
    plt.contourf(yy,zz,psi)
    plt.colorbar()
    plt.subplot(2,2,3)
    # plt.contourf(yv[1:-1,1:-1],zv[1:-1,1:-1],v[1:-1,1:-1])
    plt.contourf(yy[1:-1],zz[1:-1],v[1:-1,1:-1])
    plt.title('v')
    plt.colorbar()
    plt.subplot(2,2,4)
    plt.contourf(yv[1:-1,1:-1],zv[1:-1,1:-1],w[1:-1,1:-1])
    plt.title('W')
    plt.colorbar()
    plt.subplot(2,2,2)
    plt.contourf(yv,zv,T)
    plt.colorbar()
    plt.title('Temperature')
    show_plot()
    plt.savefig('wvsetup.pdf')
pass


#  set boundary conditions on water vapour)

wv[0,:]=0.
wv[-1,:]=0.
wv[:,0]=0
wv[:,-1]=0

#  Apply a no-flux boundary condition
#  This needs to be in the main loop
wv[:,0]  = wv[:,1]
wv[:,-1] = wv[:,-2]
wv[:,0]  = wv[:,1]
wv[-1,:] = wv[-2,:]
wv[0,:] = wvsat[0,:]  # except at the bottom


ifig += 1
plt.figure(ifig)
plt.clf()
#if 'qt' in plt.get_backend().lower():
#    plt.get_current_fig_manager().window.raise_()
#pass
show_plot()

#  Now enter main time stepping loop
for n in range(0, nt):
    diffwv = diffuse_wv(diffwv, wv_old, Ky, Kz, rdy2, rdz2)
    advectwv = advect_wv(advectwv,v,w,wv,rdy,rdz)
    S = condense(S,wvsat,wv_old,myzeros,rtau)
    wv_new[1:-1,1:-1] = wv_old[1:-1,1:-1] + \
       (0.*S[1:-1,1:-1]+ diffwv[1:-1,1:-1] - advectwv[1:-1,1:-1])*dt
    wv_old[:,:] = wv[:,:]   #  slightly clunky alogithm
    wv_new[1:-1,1:-1] = np.minimum(wvsat[1:-1,1:-1],wv_new[1:-1,1:-1])
    wv[:,:] = wv_new[:,:]

#  Now set a no-flux boundary condition This needs to be in the main loop
    #  wv[:,0]  = wv[:,1]
    wv[:,-1] = wv[:,-2]
    wv[:,0]  = wv[:,1]
    wv[-1,:] = wv[-2,:]
    wv[0,:] = wvsat[0,:]  # Dirichlet BC at the bottom

    #  Evaluate relative humidity
    rh[:,:] = wv[:,:]/wvsat[:,:]
    logwv = np.log10(abs(wv)+ eps)

    # Plot every nplot time steps
    nplot = 500
    if (n % nplot == 0):
        plt.subplot(2,2,1)
        plotlabel = "Water vapour, t = %1.2f" % (n * dt)
        #  plt.pcolormesh(yy,zz,wv, shading='flat')
        plt.contourf(yy,zz,logwv,shading='flat')
        plt.contour(yv,zv,psi,colors='r')
        plt.title(plotlabel)
        plt.axis('image')
        plt.draw()
        #if 'qt' in plt.get_backend().lower():
        #    QtGui.qApp.processEvents()
        plt.subplot(2, 2, 2)
        plt.pause(0.1)
        plotlabel = "Relativive humidity, t = %1.2f" % (n * dt)
        #  plt.pcolormesh(yy,zz,wv, shading='flat')
        plt.contourf(yy, zz, rh, shading='flat')
        plt.contour(yv, zv, psi, colors='r')
        plt.title(plotlabel)
        plt.axis('image')
        # plt.show()
        # plt.pause(0.1)
        show_plot()
        print ('Water vapour, n = ', n)
pass

show_plot()
print ('Loop finished, nt = ', nt) #  end looping

ifig += 1
plt.figure(ifig)
plt.clf()
plt.subplot(2,2,1)
plotlabel = "Log wv at t = %1.2f" %(n * dt)
# plt.pcolormesh(yy,zz,wv, shading='flat')
plt.contourf(yy,zz,logwv,shading='flat')
plt.colorbar()
plt.title(plotlabel)
# plt.axis('image')
plt.draw()
# if 'qt' in plt.get_backend().lower():
#       QtGui.qApp.processEvents()

plt.subplot(2,2,2)
# plt.pcolormesh(yy,zz,wv, shading='flat')
plt.contourf(yv[1:-1,1:-1],zv[1:-1,1:-1],rh[1:-1,1:-1],shading='flat')
plt.title(plotlabel)
# plt.axis('image')
plt.draw()
plt.colorbar()
plt.title("Relative humidity.")

plt.subplot(2,2,3)
# plt.pcolormesh(yy,zz,wv, shading='flat')
plt.contourf(yv[1:-1,1:-1],zv[1:-1,1:-1],T[1:-1,1:-1],shading='flat')
plt.title(plotlabel)
#  plt.axis('image')
plt.colorbar()
plt.draw()
plt.title("Temperature")

plt.subplot(2,2,4)
plt.title("Streamfunction")
plt.contourf(yv,zv,psi,shading='flat')
#  plt.axis('image')
plt.colorbar()
show_plot()

ifig += 1
ifig=3
plt.set_cmap('Blues')
plt.figure(ifig)
plt.clf()
plt.subplot(2,2,1)
rhlevs = np.arange(0,1.2,0.2)
plotlabel = "$\mathcal{H}$ and $\psi$"
plotlabel = "Relative humidity"
#  plt.pcolormesh(yv[1:-1,1:-1],zv[1:-1,1:-1],rh[1:-1,1:-1],shading='flat')
plt.contourf(yv[1:-1,1:-1],zv[1:-1,1:-1],rh[1:-1,1:-1],levels=rhlevs, \
       shading='flat')
plt.colorbar(ticks=np.arange(0,1.2,0.2))
plt.xticks(np.arange(0.0,1.2,0.5))
plt.yticks(np.arange(0.0,1.2,0.5))
#  cbar = plt.colorbar(cax, ticks=[-1, 0, 1])
#  cbar.ax.set_yticklabels(['< -1', '0', '> 1'])#   vertically oriented colorbar
#  plt.contour(yv[1:-1,1:-1],zv[1:-1,1:-1],psi[1:-1,1:-1],colors='b')
plt.contour(yv,zv,psi,colors='r',linewidth=0.5)
#  plt.contour(yy,zz,'psi', 'r-')
plt.title(plotlabel)
plt.xlabel('y',style='italic')
plt.ylabel('z',style='italic')
#  plt.axis('image')
#
qmin = np.nanmin(wv)
qmax = np.nanmax(wv)
#  qmin = 3.e-6
wvlev = np.logspace(np.log10(qmin),np.log10(qmax),num=9)
plt.subplot(2,2,2)
plotlabel = "Water Vapour"
#  plt.contourf(yy,zz,logwv,shading='flat')
#  plt.pcolormesh(yy,zz,logwv,shading='flat')
plt.contourf(yy,zz,wv,shading='flat', \
    locator=mpl.ticker.LogLocator(),levels=wvlev)
#  plt.colorbar()
cbar=plt.colorbar(format='%.1e')
cbar.ax.tick_params(labelsize=12)
plt.contour(yv,zv,psi,colors='r',linewidth=0.5)
plt.xticks(np.arange(0.0,1.2,0.5))
plt.yticks(np.arange(0.0,1.2,0.5))
plt.title(plotlabel)
plt.xlabel('y',style='italic')
#  plt.ylabel('z',style='italic')
#  plt.axis('image')
show_plot()


plt.savefig('wvcell.pdf')



plt.show()
