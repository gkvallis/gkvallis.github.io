clear all
% This calculates the linear response to an isolated mountain, like 
% a Gaussian hump or the Witch of Agnesi. 

% Code by G. K. Vallis

% Change the values of the parameters wide and hydro to see different 
% solutions.
% Set the parameter agnesi to set Witch of Agnesi or Gaussian.

% set(0,'DefaultAxesFontName','Myriad Pro') % comment out if needed
% set(gca, 'FontName', 'Myriad Pro')

% Set resolution and domain
nmax = 128 ;
xmax = nmax;
zmax = 64;
kmax = nmax; 
lx = 2*pi*(1-1/xmax);


% Choose whether hump is thin or wide
% giving a shallow or deep solution respectively. 
wide = 1   %  set wide = 1 for wide mountain, otherwise 0

% For hydrostatic flow over a hump set hydro = 1
% regardless of setting of shallow or deep
hydro = 1

% Choose Witch of Agnesi (agnesi = 1) or Gaussian hump (agnesi = 0).
agnesi= 1

% Common parameters for deep or shallow:
    U = 1;
    N = 5. *1.25;
    lz = 2 ;
    h0 = 1.2*U/N ;
if(wide)
    a = 2.*U/N;  % long waves 
else
    a = U/(2*N); % short waves
end

% Hydrostatic case:
if(hydro)
    U = 1. ;
    N = 40. ;
    a = 10.*U/N ;
    h0 = 0.7*U/N ;
    lz = lx/30 ;
end


% --------------------------

% Now find solution:
x = linspace(0,lx,xmax);
z = linspace(0.,lz,zmax);
%z2 = linspace(-0.5,lz,zmax);
[X,Z] = meshgrid(x,z);
unit = linspace(-0.22,-0.22,xmax);
unit2 = linspace(-0.2,-0.2,xmax);
zero = linspace(0,0,xmax);

hgauss = h0*exp(-(x-pi).^2/a^2);
hgauss = hgauss  - mean(hgauss);

hagn = h0*a^2./(a^2 + (x - 0.5*lx).^2) ;
hagn = hagn  - mean(hagn);

if(agnesi)
    h = hagn;    % Witch of Agnesi solution
else
    h = hgauss;  % Gaussian hump
end

hsize = size(h);
 
% For a contour plot:
for p = 1:xmax
    for q = 1:zmax
        if Z(q,p) > h(p)
            H(q,p) = 0 ;
        else
            H(q,p) = 0.2  ;
        end
    end
end


% Now calculate FFT for the mountain
hf = fft(h) ;
% check for inverse
h2 = ifft(hf);

for j=1:nmax/2
    k = (j -1) ;
    dhf(j) = i*k*hf(j) ;  
end
for j=nmax/2+1:nmax
    k= j - (nmax+1) ; 
    dhf(j) = i*k*hf(j) ;
end


% Now calculate solutions for each Fourier mode, part I
 for j=1:nmax/2
     k = j -1 ;
     % Fist calculate the vertical wavenumner
     m2 = (N^2/U^2 - k^2);
     mr = -real(sqrt(m2));
     mi = -imag(sqrt(m2));
     m(j) = mr + i*mi;
     %
     w1(j) = -i*U*k*hf(j); 
     u1(j) = i*U*m(j)*hf(j);
     p1(j) = -i*U^2*m(j)*hf(j);
     psi1(j) = U*hf(j);
 % Now calculate the vertical variation for each wavenumber
   for q=1:zmax
       wkz(j,q) = w1(j)*exp(-i*m(j)*z(q));
       ukz(j,q) = u1(j)*exp(-i*m(j)*z(q));
       pkz(j,q) = p1(j)*exp(-i*m(j)*z(q));
       psikz(j,q) = psi1(j)*exp(-i*m(j)*z(q));
   end
 end
 % Now calculate solutions for each Fourier mode, part II
 for j=nmax/2 + 1:nmax
     k = j - (nmax+1);
     % First calculate the vertical wavenumner
     m2 = (N^2/U^2 - k^2);
     mr = -real(sqrt(m2));
     mi = imag(sqrt(m2));
     m(j) = mr + i*mi;
     %
     w1(j) = i*U*k*hf(j); 
     u1(j) = -i*U*m(j)*hf(j);
     p1(j) = i*U^2*m(j)*hf(j);
     psi1(j) = U*hf(j);
 % Now calculate the vertical variation for each wavenumber
   for q=1:zmax
       wkz(j,q) = w1(j)*exp(i*m(j)*z(q));
       ukz(j,q) = u1(j)*exp(i*m(j)*z(q));
       pkz(j,q) = p1(j)*exp(i*m(j)*z(q));
       psikz(j,q) = psi1(j)*exp(i*m(j)*z(q));
   end
 end
  
% Now calculate the physical space fields
% wxzsize=size(wxz)
 for q = 1:zmax
     wxz(q,1:xmax) = ifft(wkz(1:nmax,q));
     uxz(q,1:xmax) = ifft(ukz(1:nmax,q));
     pxz(q,1:xmax) = ifft(pkz(1:nmax,q));
     psixz(q,1:xmax) = ifft(psikz(1:nmax,q));
 end
 
 psiP = psixz(2:zmax,:);
 ZP    = Z(2:zmax,:);
 XP    = X(2:zmax,:);
 
  wxzsize=size(wxz);
 % Now reconstruct the z variation;
 
minusp = -pxz ;
psitot = psixz -1*U*Z;
psitotP = psiP - U*ZP;

labelfont=18 ;
axisfont=14 ;

% Here we arbitrarily plot minus pressure just to get nicer shading
figure(1)
subplot(2,1,1)
hold off
contourf(X,Z,-real(pxz),4,'k')  % use 5 or 6 for long waves
set(gca,'Fontsize',labelfont)
%contourf(XP,ZP,real(pxzP),5,'k')
colormap gray(128)
temp = colormap;
colormap2 = temp(96:128,:);
colormap(colormap2)
ylabel('Height,\it z','Fontsize',labelfont)
title('Pressure','Fontsize',labelfont,'Fontweight','normal')

hold on
subplot(2,1,2)
hold off
contour(X,Z,real(psitot),8,'k')  % 8 is good for long waves
set(gca,'Fontsize',labelfont)
hold on
hh = area(x,h,-0.2) ;
set(hh,'FaceColor',[0.3 0.3 0.3])
axis([0 lx -lz/50. lz])
xlabel('Horizontal distance,\it x','Fontsize',labelfont)
ylabel('Height,\it z','Fontsize',labelfont)
title('Streamfunction','Fontsize',labelfont,'Fontweight','normal')
%set(gca, 'FontName', 'Times')
set(gcf,'color','w');

% Uncomment one of the following if you wish, to print a PDF or EPS
%savefig('hump','pdf')
%savefig('humpe','eps')
%saveas(gcf,'hump3', 'pdf')
% export_fig  hump1.pdf
% export_fig  hump1.eps
print -depsc2 hump.eps
print -dpdf hump.pdf

return
