% Steps forward the linear shallow water model for Gill-like solutions
% By G. K. Vallis 2011 and onward.
% The code is not very clean and is not well commented. 
% But it is very straightforward. 
% This version steps forward the canonical Gill problem and plots it.

% First set some domain and grid parameters
clear all
cputime_begin=cputime
tic
imax=160;   % 160 is a good default
jmax=80;    % 80 is a good default
imax = 80; jmax = 40;  % for testing 
imaxd2 = imax/2;
jmaxd2 = jmax/2;
xmax = 20; 
ymax = 8;
xmin =-16;
ymin =-8;
dx = (xmax - xmin)/(imax-3)
dy = (ymax - ymin)/(jmax-1)
dx2 = 2.*dx;
dy2 = 2.*dy;
rdx2=1./dx2;
rdy2=1./dy2;

% Now set some physical and timestepping parameters
dt = 0.02;   % 0.02 works fof default imax, jmax
dt = 0.04
mmax = 10000;
ass = 0.02;
raw=0.5;
alphah = 0.15;
alphau = 0.15;
alphav = 0.15;
beta  = 0.5;
L = 2;
k = pi/(2*L);
% nu = 0.002;
 nu = 0.0 ;


for nt = 1:2
yn = nt
FQ = zeros(imax,jmax); FU = zeros(imax,jmax); FV = zeros(imax,jmax);
FU1 = zeros(imax,jmax); FQ2 = zeros(imax,jmax);

% Now some grid variables
for i = 1:imax
for j = 1:jmax
   x = (xmax-xmin)*(i-2)/(imax-3) + xmin;
   y = (ymax-ymin)*(j-1)/(jmax-1) + ymin; 
   expy = exp(-y^2/4);
   expx = exp(-x^2/2);
   xx(i,j) = x;
   yy(i,j) = y;
   y1(j) = y;
   f(i,j) = beta*yy(i,j);
% Evaluate forcing
  if x < -L 
     FQ(i,j) = 0;
  elseif x > L
     FQ(i,j) = 0;
  else
      FQ(i,j) = -y*cos(k*x)*expy - cos(k*x)*expy;  % Gill original
%     FU1(i,j) = cos(k*x)*expy;
%     FV(i,j) = 0.*cos(k*x)*expy;
  end
%      fqq  = -1.*exp(-(y-yn)^2/4)*exp(-x^2/2) - 0.*exp(-(y+2)^2/2)*exp(-x^2/2);
%     FU(i,j) = - exp(-(y)^2/4)*exp(-x^2/2);
     if nt == 1    
         FQ(i,j) = -cos(k*x)*expy;
         FQ2(i,30) = -1.;    % For line source
         FQ2(i,31) = -1.;
     elseif nt == 2
         FQ(i,j) = - y*cos(k*x)*expy ; 
         FQ2(i,35) = -1;
         FQ2(i,36) = -1;
     end
      if x < -L 
     FQ(i,j) = 0;
  elseif x > L
     FQ(i,j) = 0;
      end
        
%      if y < 0.9
%           FQ2(i,j) = 0;
%      elseif y > 1.1
%          FQ2(i,j) = 0;
%      end
%     FQ(i,j) =  FQ2(i,j) ; % + fqq;
%     FQ(i,j) = fqq; 
end
end

FQ(:,1) = 0.;
FQ(:,jmax)=0 ;

% Initialize variables
  u(1:imax,1:jmax) = 0.;
  v(1:imax,1:jmax) = 0.;
  h(1:imax,1:jmax) = 0.;
  uo(1:imax,1:jmax) = 0.;
  vo(1:imax,1:jmax) = 0.;
  ho(1:imax,1:jmax) = 0.;
  dhx1 = zeros(imax,jmax); dhx2 = zeros(imax,jmax); dhx = zeros(imax,jmax);
  dhy1 = zeros(imax,jmax); dhy2 = zeros(imax,jmax); dhy = zeros(imax,jmax);
  dux1 = zeros(imax,jmax); dux2 = zeros(imax,jmax); dux = zeros(imax,jmax);
  dvy1 = zeros(imax,jmax); dvy2 = zeros(imax,jmax); dvy = zeros(imax,jmax);
  urhs = zeros(imax,jmax); vrhs = zeros(imax,jmax); hrhs = zeros(imax,jmax);
  uass = zeros(imax,jmax); vass = zeros(imax,jmax); hass = zeros(imax,jmax);
  laplh = zeros(imax,jmax);

  % Now the timestepping
  % Do this in vector form for speed. 
for m = 1:mmax; 
    dhx1(2:imax-1,:) = diff(h(2:imax,:));
    dhx2(2:imax-1,:) = diff(h(1:imax-1,:));
    dhx(2:imax-1,:)  = (dhx1(2:imax-1,:) + dhx2(2:imax-1,:))*rdx2;
    
    dhy1(:,2:jmax-1) = diff(h(:,2:jmax),1,2);
    dhy2(:,2:jmax-1) = diff(h(:,1:jmax-1),1,2);
    dhy(:,2:jmax-1)  = (dhy1(:,2:jmax-1) + dhy2(:,2:jmax-1))*rdy2;
    
    
    dux1(2:imax-1,:) = diff(u(2:imax,:));
    dux2(2:imax-1,:) = diff(u(1:imax-1,:));
    dux(2:imax-1,:)  = (dux1(2:imax-1,:) + dux2(2:imax-1,:))*rdx2; 
    
    dvy1(:,2:jmax-1) = diff(v(:,2:jmax),1,2);
    dvy2(:,2:jmax-1) = diff(v(:,1:jmax-1),1,2);
    dvy(:,2:jmax-1)  = (dvy1(:,2:jmax-1) + dvy2(:,2:jmax-1))*rdy2;
    
    laplh = del2(ho,dx,dy);  % Lapacian
    
    urhs =  f(:,:).*v(:,:) - dhx(:,:) - alphau*uo(:,:) + FU(:,:);
    vrhs = -f(:,:).*u(:,:) - dhy - alphav*vo(:,:) + FV(:,:);
    hrhs = -(dux(:,:) + dvy(:,:)) + FQ(:,:) - alphah*ho(:,:) + nu*laplh;
 
    % Now step forward   
    un = uo + urhs*dt;
    vn = vo + vrhs*dt;
    hn = h + hrhs*dt;
    
  % Now the Raw filter
    uass = ass*(un + uo - 2*u);
    vass = ass*(vn + vo - 2*v);
    hass = ass*(hn + ho - 2*h);
    uo = u + raw*uass;
    vo = v + raw*vass;
    ho = h + raw*hass;

  % Now update the new fields
    u = un + (raw-1)*uass;
    v = vn + (raw-1)*vass; 
    h = hn + (raw-1)*hass;


  % Now do the periodic boundary values at the zonal ends
  % Periodic
    u(imax,:) = u(3,:);
    v(imax,:) = v(3,:);
    h(imax,:) = h(3,:);
    u(1,:) = u(imax-2,:);
    v(1,:) = v(imax-2,:);
    h(1,:) = h(imax-2,:);
    uo(imax,:) = uo(3,:);
    vo(imax,:) = vo(3,:);
    ho(imax,:) = ho(3,:);
    uo(1,:) = uo(imax-2,:);
    vo(1,:) = vo(imax-2,:);
    ho(1,:) = ho(imax-2,:);    
end

% Now calculate vertical velocity
w(:,:) = alphah*h(:,:) - FQ(:,:);
    

% Now decimate the field for quiver plots
incrx=imax/20;
incry=jmax/20;
ii=0; jj=0;
for i = 1:incrx:imax
    ii = ii+1;
    jj=0;
    for j=2:incry:jmax
        jj = jj+1;  
        xxd(ii,jj) = xx(i,j);
        yyd(ii,jj) = yy(i,j);
        ud(ii,jj) = u(i,j);
        vd(ii,jj) = v(i,j);
    end
end

% Now a smaller decimation to get round a Matlab quirk in dashed contours
incrx=2;
incry=1;
ii=0; jj=0;
for i = 1:incrx:imax
    ii = ii+1;
    jj=0;
    for j=2:incry:jmax
        jj = jj+1; 
        xxd2(ii,jj) = xx(i,j);
        yyd2(ii,jj) = yy(i,j);
        hd2(ii,jj) = h(i,j);
    end
end

for i=1:imax
 h(i,:) = smooth(h(i,:),3);
end

hh = zeros(imax,jmax);
 hh = -h; 
%Now plot some solutions
figure(nt)
hold off

%subplot(1,2,nt)
hold off
contourf(xx,yy,h,10,'k')
%contour(xxd2,yyd2,hd2,[0 -0.3 -0.6 -0.9 -1.2 -1.5 -1.8],'k')
%contour(xxd2,yyd2,hd2,[-2.25:0.3:-0.15],'k')
%contour(xx,yy,h,[-2.25:0.3:-0.15],'k')
%contourf(xx,yy,hh,[0.9:0.3:1.5],'k-')
hold on
%contour(xx,yy,hh,[0.3:0.3:0.6],'b')
%contourf(xx,yy,h,[0.5:0.5:5],'b')
%contour(xxd2,yyd2,hd2,[0.15:0.3:2.25],'b--','LineWidth',0.7)
%set(gcf,'renderer','painter');
hold on
quiver(xxd,yyd,ud,vd,0.65,'b')
title('Pressure and Velocity','fontsize',8)
axis([xmin xmax ymin ymax])
ylabel('\it y','fontsize',9)
set(gca,'fontsize',8)
xlabel('\it x','fontsize',9)
pbaspect([2 1 1])

colormap('gray')
colormap2 = (colormap+6)/7. ;
colormap(colormap2)

end



% Print if needs be:  (uses export_fig and print2eps)
% export_fig -transparent 'Gill.pdf'
% print2eps 'Gill.eps'

toc

return

% --------------------- Code stops here -----------------------------------


figure(2)
hold off
contour(xx,yy,w,[0.2 0.4 0.6 0.8 1.0 1.2],'k')
%contour(xx,yy,w,8)
hold on
contour(xx,yy,w, [0 0],'k--','LineWidth',0.7)
contour(xx,yy,w, [-0.1 -0.1],'k.') %,'LineWidth',1)
pbaspect([2 1 1])
% colormap('gray')
% colormap2 = (colormap+4)/5. ;
% colormap(colormap2)



cpu=cputime - cputime_begin 

figure(1)
hold off

hh = zeros(imax,jmax);
for i=1:imax
 hh(i,:) = -smooth(h(i,:),3);
end
subplot(3,2,nt)
hold off
%contourf(xx,yy,h,5,'k')
%contour(xxd2,yyd2,hd2,[0 -0.3 -0.6 -0.9 -1.2 -1.5 -1.8],'k')
%contour(xxd2,yyd2,hd2,[-2.25:0.3:-0.15],'k')
%contour(xx,yy,hh,[-2.25:0.3:-0.15],'k')
contourf(xx,yy,hh,[0.8:0.4:2.0],'k-')
hold on
contour(xx,yy,hh,[0.4 0.4],'b')
%contour(xxd2,yyd2,hd2,[0.15:0.3:2.25],'b--','LineWidth',0.7)
%set(gcf,'renderer','painter');
hold on
quiver(xxd,yyd,ud,vd,0.65,'b')
title('Pressure and Velocity','fontsize',8)
axis([xmin xmax ymin ymax])
axis([-10 15 -5 5])
ylabel('\it y','fontsize',9)
set(gca,'fontsize',8)
xlabel('\it x','fontsize',9)
pbaspect([2 1 1])

colormap('jet')
colormap2 = (colormap+4)/5. ;
colormap(colormap2)
print -depsc2 Gill_line1.eps


!skimm Gill_line1.pdf &

return


