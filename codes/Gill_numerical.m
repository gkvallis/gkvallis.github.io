% Script to solve Gill problem.

% Original code by C. Bretherton
% with modifications by G. Vallis

% We use an FFT 
% with periodic BC in x and finite differencing in y.
% with BC v = 0 at N, S walls.

% You can set whatever heating or cooling you wish if you modify the 
% code appropriately. 

clear all
% mymap = cbrewer('div','Spectral',16);
%mymap = cbrewer('div', 'RdBu',8);
%colormap(mymap)

colormap('parula')
% mymap = brighten(0.4) ; 



%-------------User-defined parameters follow...--------------

 zonalcomp = 0; % No zonal compensation of heating (Original Gill problem)
%  zonalcomp = 1; % Gill problem with zonall compensated heating
%  nodiss = 1; % Inviscid, no thermal diffusion  
 nodiss = 0; % Thermal diffusivity = Rayleigh damping (Original Gill problem)

%  Set Rayleigh friction (nondimensionalized in units of sqrt(beta*c))
 if(nodiss)     % that is, if nodiss = 1
   a = 0.0001;  % viscosity, a = 0 would produce divide by zero.
 else
   a = 0.15;  
 end

%  Define other physical parameters
 H = 1;     % Layer depth
 g = 1;     % gravity
 beta = 1;  % df/dy

%  Domain size and number of mesh points
  % 30 and 10 are canonical values for lx and ly
 lx = 50; % periodic domain width -lx/2 < x < lx/2
 ly = 20; % Rigid walls at +/- ly/2
 hires=1 ;
 nx = hires*256; % Number of gridpoints in x
 ny = hires*120; % Number of y gridpoints in -ly/2 <= y < ly/2

 
%  Parameters defining plot window in x,y and stride (in gridpoints) for
%  velocity vector plotting

 xmin = -lx/2;
 xmax = lx/2;
 xmin = -15 ;
 xmax = 15 ; 
 ymin = -5;
 ymax = 5;
 
 stride = 3*hires;

%----------End of basic parameters---------  

 dx = lx/nx;
 x = -lx/2 + dx*(0:(nx-1));
 dy = ly/ny;
 y = -ly/2 + dy*(0:ny);
 [X,Y] = meshgrid(x,y);
 myzero = 0.*x  ;

% Shallow water phase speed
c = sqrt(g*H);

%   Define mass source M(x,y)
% for iy = 1:2
iy = 1 ;
 sx = 2;    % Mass source half-width in x
 sy = 1;    % Mass source half-width in y (1 is default)
 %sy = 0.1;  % Use 0.1 this for line source
 x0 = -2;    % Central x of mass source
 x1 = +6;
 y0 = iy-1;    % Central y of mass source.
 y1 = y0;

 kh = pi/(2*sx);
 phase = kh*(X-x0);
 phase(X-x0>sx) = pi/2;
 phase(X-x0<-sx) = pi/2;
 F = cos(phase);  % source 1
 
 phase1 = kh*(X-x1);
 phase1(X-x1>sx) = pi/2;
 phase1(X-x1<-sx) = pi/2;
 F1 = cos(phase1);   % source 2
 
 %F = 1;  % Use this for a line source
 
 % Use the following for a separate heat source and a sink
 M0 = - F.*exp(-(Y-y0).^2/(2*sy^2)); 
 M1 = +F1.*exp(-(Y-y1).^2/(2*sy^2));
 M = M0 + M1 ; 
 
% Use the following for a zonal temperature gradient
%  M = A*(Y -x1 - x0)
 Mag = -0.4  ;  % or -1 or  -0.7 ;
 sxx = 15 ;
 ZM = Mag/(x0 - x1) * (2*X - (x0 + x1)).*exp(-(Y-y0).^2/(2*sy^2));;
 for i = 1:ny
     for j = 1:nx
         if X(i,j) > x1
             ssy = 1 + 0.2*abs(X(i,j) - x1);
             ssy = sy;
             ZM(i,j) = -Mag*exp(-(X(i,j)-x1)^2/sxx^2)...
                 *exp(-(Y(i,j)-y0).^2/(2*ssy^2));
         end
         if X(i,j) < x0
             ssy = 1 + 0.2*abs(X(i,j) - x0);
             ssy = sy;
             ZM(i,j) =  Mag*exp(-(X(i,j)-x0)^2/sxx^2)...
                 *exp(-(Y(i,j)-y0).^2/(2*ssy^2));
         end
     end
 end
M = ZM;

% Following for the canonical Gill problem
igill = 1
if(igill) 
 iy = 1
 sx = 2;    % Mass source half-width in x
 sy = 1;    % Mass source half-width in y (1 is default)
 %sy = 0.1;  % Use 0.1 this for line source
 x0 = 0;    % Central x of mass source
 y0 = iy-1;    % Central y of mass source.

 kh = pi/(2*sx);
 phase = kh*(X-x0);
 phase(X-x0>sx) = pi/2;
 phase(X-x0<-sx) = pi/2;
 F = cos(phase);   % Use this for the Gill model
 %F = 1;  % Use this for a line source
 M =  F.*exp(-(Y-y0).^2/(2*sy^2));
end
% done with Gill problem 

% now do some stuff for all the forcing
M(1,:) = 0;                % Zero the heating at top and bottom bdries.
 M(ny+1,:) = 0;
 if (zonalcomp)   % if zonalcomp = 1
   M = M - (mean(M.')).'*ones(1,nx); % Correct zonal avg heating to be zero
 end
 
  % M = M - (mean(M.')).'*ones(1,nx);; % compensate anyway
   
   
 dMdy = -(Y-y0).*M/sy^2;
 d2Mdy2 = ((Y-y0).^2/sy^2 - 1).*M/sy^2;
 Mhat = (fft(M.')).';
 dMdyhat = (fft(dMdy.')).';
 
%  Define waveno. matrix 

 kx = (2*pi/lx)*[0:(nx/2 - 1) (-nx/2):(-1)];
 KX = ones(ny+1,1)*kx;

%-------------------Now find the solution --------------------------------

%  Define thermal diffusivity b
b = a;

%  Define v source term Sv = (a*d/dy - beta*y*d/dx)M/H
 Svhat = (a*dMdyhat - beta*1i*KX.*Y.*Mhat)/H;

%  Solve 
%  (-b(a^2 + beta^2y^2)/c^2 + a*del^2 + beta*d/dx)v = Sv,
%  or
%  a*d2vhat/dy2 + (-b(a^2 + beta^2y^2)/c^2 -a*k^2 + i*k*beta)vhat = Svhat
%
%  where Sv = (a*d/dy - beta*y*d/dx)M/H is the same source term
%  as in WTG (but this time we don't remove wavenumber zero).
%  This is done as a loop over wavenumbers, using only the interior 
%  y-gridpoints 2:ny (since v = 0 at the boundaries)

 vhat = zeros(ny+1,nx);
 d1 = a/dy^2;
 for i = 1:nx
   k = kx(i);
   d0 = -2*d1 - a*k^2 + 1i*k*beta;
   e = ones(ny-1,1);
   Av = spdiags([d1*e d0*e-b*(a^2 + (beta*(y(2:ny))').^2)/c^2 d1*e],...
                -1:1,ny-1,ny-1);
   r = Svhat(2:ny,i);
   vhat(2:ny,i) = Av\r;
 end
 v = real((ifft(vhat.')).');

%  Calculate phi from
%     (b/c^2)phi + du/dx = M/H - dv/dy 
%  and
%     a*u = -dphi/dx + beta*y*v
%  Eliminating u between these equations,
%     (a*b/c^2 - d2/dx2)phi = aM/H - a dv/dy - beta*y*dv/dx
%  whose FFT in x diagnoses phi

 dvdyhat = zeros(ny+1,nx);
 dvdyhat(2:ny,:) = (vhat(3:(ny+1),:) - vhat(1:(ny-1),:))/(2*dy);
 dvdyhat(1,:) = (vhat(2,:) - vhat(1,:))/dy;
 dvdyhat(ny+1,:) = (vhat(ny+1,:) - vhat(ny,:))/dy;
 phihat = (a*Mhat/H - a*dvdyhat - 1i*beta*Y.*KX.*vhat)./(a*b/c^2 + KX.^2);
 phi = real((ifft(phihat.')).');

 D = M/H - b*phi/c^2;


%  Calculate vorticity zeta from divergence and v:
%    a*zeta + beta*(y*D + v) = 0

 zeta = (beta/a)*(-Y.*D - v);

%  u calculated using div eqn: du/dx + dv/dy = M/H-b*phi/(c^2)

 dvdy = zeros(ny+1,nx);
 dvdy(2:ny,:) = (v(3:ny+1,:)  - v(1:(ny-1),:))/(2*dy);

 uhat = (Mhat/H - b*phihat/c^2-(fft(dvdy.')).');
%   uhat=(Mhat/H-(fft(dvdy.')).');				  

 uhat(:,2:nx) = uhat(:,2:nx)./(1i*KX(:,2:nx));

%  The k=0 components are indeterminate; for these go to zonally
%  averaged vorticity equation dudyhat(:,1) = -zetahat(:,1), with BC that
%  the meridional average of uhat(:,1) should equal zero. 

 zetahat = (fft(zeta.')).';
 dudyhath = -0.5*(zetahat(1:ny,1) + zetahat(2:(ny+1),1)); 
 uhat(:,1) = [0; dy*cumsum(dudyhath)];
 uhatmean = mean([uhat(2:ny,1);0.5*(uhat(1,1)+uhat(ny+1,1))]);
 uhat(:,1) = uhat(:,1) - uhatmean;

 ucompare = real((ifft(uhat.')).');
	    

%  u calculated using x momentum eqn: -beta.y.v=-d(phi)/dx-a.u 
		   dphidx=zeros(ny+1,nx);  

 dphidx(:,2:nx-1)=(phi(:,3:nx)-phi(:,1:nx-2))/(2*dx);
 dphidx(:,1)=(phi(:,2)-phi(:,nx))/(2*dx);
 dphidx(:,nx)=(phi(:,1)-phi(:,nx-1))/(2*dx);
 u=(-beta*Y.*v+dphidx)/(-a);

 % if (~nodiss)




%-------Various plots--------
ifig = 0 ;
SQ = 0.3*max(abs(u(ny/2+1,:)));
ifig = ifig+1 
figure(ifig)
hold off
subplot(2,1,2)
hold off
plot(x,-M(ny/2+1,:),'linewidth',2)
hold on
plot(x,u(ny/2+1,:)/4,'r--','linewidth',2)
plot(x,myzero,'k-','linewidth',0.5)
xlim([xmin xmax])
ylim([-1.2 1.2])
xlabel('x/L_{eq}')
ylabel('Heating')
pbaspect([10,2,1])   
hold off
% colormap(mymap)
    subplot(2,1,1)
    hold off
%    Plot Gill divergence
   cDmax = max(max(D))*1.1  ;
   cint = (-1:0.2:1)*cDmax ;
   % contourf(x,y,-D,cint,'k-')
   pcolor(x,y,D) ; shading interp ; % colormap(mymap)
   axis equal
   axis([xmin xmax ymin ymax])
   xlabel('x/L_{d}')
   ylabel('y/L_{d}')
   % text(0.75*xmax+0.25*xmin,0.85*ymax+0.15*ymin,...
     %   ['a = b = ' num2str(a) 'c/L_{eq}'])
   title('Velocity and divergence','fontsize',12,'FontWeight','normal')

%    Plot velocity vectors and vorticity
hold on
   %quiver(x(1:stride:nx),y(1:stride:ny),...
   %       u(1:stride:ny,1:stride:nx),v(1:stride:ny,1:stride:nx))
   
   quiver(x(1:stride:nx),y(1:stride:ny),...
          u(1:stride:ny,1:stride:nx),v(1:stride:ny,1:stride:nx),SQ,'k')
      
   %      quiver(x(1:stride:nx),y(1:stride:ny),...)
   hold off
  % colormap(mymap)
  

  
  ifig = ifig +1 
   figure(ifig)
   subplot(2,1,1)
   cphimax = max(max(abs(phi))) ;
   %   cphimax = 2;      % fixed contours
   cpos = (0.1:0.2:1.9)*cphimax;
   contour(x,y,phi,cpos,'k-')
   hold on
   contourf(x,y,phi,cpos)
   hold on
   contour(x,y,phi,-cpos,'k--')
   axis equal
   axis([xmin xmax ymin ymax])
   xlabel('x/R_{d}')
   ylabel('y/R_{d}')
   title('Gill geopotential')
   hold off

% Comment out the printing
% export_fig -transparent 'Gillish.pdf'
% print2eps 'Gillish.eps'
   
fprintf 'finished \n' 
return   
