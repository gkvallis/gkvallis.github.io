% This code calculates the analytic solution of the canonical Gill model
% and plots the results. 
% The code is not wll documented internally. However, it is straightforward
% to read, in connjuction with a description of the Gill solution, such as
% may be found in Chapter 8 of AOFD.

% Code by G. K. Vallis

clear all


% b = repmat(linspace(0,1,200),20,1);
% imshow(b,[],'InitialMagnification','fit')

% Set the resolution
imax=200;
jmax=200;
imaxd2 = imax/2;
jmaxd2 = jmax/2;
alpha = 0.1;
beta = 1;
L = 2;
k = pi/(2*L);
alk = 1./(alpha^2 + k^2);
alk3 = 1./(9*alpha^2 + k^2);
xmax = 15; 
ymax = 5;
xmin=-10;
ymin=-5;
for i = 1:imax
for j = 1:jmax
x = (xmax-xmin)*(i-1)/(imax-1) + xmin;
y = (ymax-ymin)*(j-1)/(jmax-1) + ymin;
expy = exp(-y^2/4);
xx(i,j) = x;
yy(i,j) = y;
y1(j) = y;
zer(j) = 0.0;
xe=xmax; yn=ymax;

% Now start the calculation.
% Evaluate forcing
if x < -L 
    FQ = 0;
elseif x > L
   FQ = 0;
else
   FQ = cos(k*x);
end

% Now evaluate q0
if x < -L 
    q0 = 0;
elseif x > L
    q0 = -alk*k*(1+exp(-2*alpha*L))*exp(alpha*(L-x));
else
    q0 = -alk*(alpha*cos(k*x) + k*(sin(k*x) + exp(-alpha*(x+L))));
end

% Now use q0 to evaluate the Kelvin wave solution
uk(i,j) = 0.5*q0*expy;
pk(i,j) = uk(i,j);
vk(i,j) = 0;
wk(i,j) = 0.5*(alpha*q0 + FQ)*expy;

% Now evaluate q2, which we will use for the Rossby solution.
if x > L 
    q2 = 0;
elseif x < -L
    q2 = -alk3*k*(1+exp(-6*alpha*L))*exp(3*alpha*(x+L));  
else
    q2 = alk3*(-3*alpha*cos(k*x) + k*(sin(k*x) - exp(3*alpha*(x-L))));
end

% Now use the q2 to reconstruct the Rossby wave solution
pr(i,j) = 0.5*q2*(1+y^2)*expy;
ur(i,j) = 0.5*q2*(y^2 -3)*expy;
vr(i,j) = (FQ + 4*alpha*q2)*y*expy;
wr(i,j) = 0.5*(FQ+alpha*q2*(1+y^2))*expy;

% Now construct the Total solution
ukr(i,j) = uk(i,j) + ur(i,j);
pkr(i,j) = pk(i,j) + pr(i,j);
vkr(i,j) = vk(i,j) + vr(i,j);
wkr(i,j) = wk(i,j) + wr(i,j);

% Now, for future reference, put the q's into arrays
qq0(i) = q0;
qq2(i) = q2;

end
end




% Now decimate the field for quiver plots
% there are better ways to do this...
incr=10;
ii=0; jj=0;
for i = 1:incr:imax
    ii = ii+1;
    jj=0;
    for j=1:incr:jmax
        jj = jj+1;  
        xxd(ii,jj) = xx(i,j);
        yyd(ii,jj) = yy(i,j);
        ukrd(ii,jj) = ukr(i,j);
        vkrd(ii,jj) = vkr(i,j);
        urd(ii,jj) = ur(i,j);
        vrd(ii,jj) = vr(i,j);
        ukd(ii,jj) = uk(i,j);
        vkd(ii,jj) = vk(i,j);
    end
end


%Now plot some solutions. First the pressure
figure(1)
hold off
contourf(xx,yy,pkr,[ -0.0 -0.3 -0.6 -0.9 -1.2 -1.5 -1.8],'k')
hold on
quiver(xxd,yyd,ukrd,vkrd,0)
colormap('bone')
colormap2 = (colormap+4)/5. ;
colormap(colormap2)

pbaspect([2 1 1])
% print -depsc2 Gill_anal.eps
% !epstopdf Gill_anal.eps ; skimm Gill_anal.pdf &


% Now plot a variety of fields.
figure(2)
hold off
subplot(3,2,1)
hold off
%contourf(xx,yy,pk,5,'k')
contourf(xx,yy,pk,[  -0.3 -0.6 -0.9 -1.2 -1.5 -1.8],'k')
%contourf(xx,yy,pk,[ -0.15 -0.3 -0.45 -0.6 -0.75 -0.9 -1.05 -1.2],'k')
hold on
quiver(xxd,yyd,ukd,vkd,0)
title('Kelvin pressure','fontsize',8,'fontweight','normal')
axis([xmin xmax ymin ymax])
ylabel('\it y','fontsize',9)
set(gca,'fontsize',8)
%xlabel('\it x','fontsize',8)


subplot(3,2,2)
hold off
%contourf(xx,yy,pr,5,'k')
contourf(xx,yy,pr,[  -0.3 -0.6 -0.9 -1.2 -1.5 -1.8],'k')
hold on
%contour(xx,yy,pr,[0 0],'k')
hold on
quiver(xxd,yyd,urd,vrd,0)
title('Rossby pressure','fontsize',8,'fontweight','normal') 
%    'FontName','Cambria')   doesn't work
set(gca,'fontsize',8)
%axis([xmin xmax ymin ymax])
% set(gca,'ytick',[ ])
% set(gca,'xtick',[ ])


subplot(3,2,3)
hold off
%contourf(xx,yy,pkr,5,'k')
contourf(xx,yy,pkr,[ -0.0 -0.3 -0.6 -0.9 -1.2 -1.5 -1.8])
hold on
quiver(xxd,yyd,ukrd,vkrd,0)
title('Total pressure','fontsize',8,'fontweight','normal')
%axis([xmin xmax ymin ymax])
ylabel('\it y','fontsize',9)
xlabel('\it x','fontsize',9)
set(gca,'fontsize',8)
%pbaspect([2,1,1])


%wkr = wkr - mean(mean(wkr)) ; 

subplot(3,2,4)
%set(gca,'FontName','Myriad Pro')
hold off
%contourf(xx,yy,-wkr,[-0.3 -0.6 -0.9 -1.2 -1.5],'k')
contourf(xx,yy,-wkr,[-0.2 -0.5 -0.8 -1.1 -1.4 -1.7],'k')
hold on
contour(xx,yy,-wkr, [0 0],'r--','LineWidth',0.7)
contour(xx,yy,-wkr, [0.1 0.1],'b') 
contour(xx,yy,-wkr, [0.1 0.2],'k'); colormap('white')
% colorbar
%set(gcf,'renderer','zbuffer');
quiver(xxd,yyd,ukrd,vkrd,0)
% title('Total vertical velocity','fontsize',8,'fontweight','normal')
title('Total vertical velocity','fontsize',8,'fontweight','normal',...
    'FontName','Myriad Pro')
%axis([xmin xmax ymin ymax])
%pbaspect([2,1,1])
xlabel('\it x','fontsize',9)
set(gca,'fontsize',8)



%colormap('summer')
%beta=.5
%brighten(beta)

% colormap('gray')
colormap('bone')
colormap2 = (colormap+4)/5. ;
colormap(colormap2)
%brighten(1)

% export_fig -transparent 'Gill-anal.pdf'
% print2eps ('Gill_anal.eps')
% print -dpdf Gill_anal.pdf
% /users/gkv/scripts/skimm Gill-anal.pdf
% export_fig 'Gillanal.eps'  % does not work. 

return
 
% Now plot the overturning streamfunctions
% In fact I donn't use this code. 

mmax=20;
mmaxd2 = mmax/2;
zmax = 1;
zmin=0.0
for i = 1:imax
for m = 1:mmax 
    x = (xmax-xmin)*(i-1)/(imax-1) + xmin;
    z = (zmax-zmin)*(m-1)/(mmax-1) + zmin;
    ssn = sin(pi*z/zmax);
    psiz(i,m) = (1/sqrt(pi))*(qq0(i) - qq2(i))*ssn;
    xx2(i,m) = x;
    zz2(i,m) = z;
end
end
for j=1:jmax
for m=1:mmax
    y = (ymax-ymin)*(j-1)/(jmax-1) + ymin;
    z = (zmax-zmin)*(m-1)/(mmax-1) + zmin;
    ssn = sin(pi*z/zmax);
    psim(j,m) = (-y/3.)*exp(-y^2/4)*ssn;
    yy2(j,m) = y;
    zz2(j,m) = z;
end
end
    
figure(3)
hold off
subplot(2,1,1)
contour(xx2,zz2,psiz,8,'b');
xlabel('\it x')
ylabel('\it z')
pbaspect([2 1 1])
axis([xmin xmax zmin zmax])

subplot(2,1,2)
contour(yy2,zz2,psim,8,'k');
xlabel('\it y') 
ylabel('\it z')
pbaspect([2 1 1])

% Print if you want to: 
 %print -dpdf Gill-anal.pdf
 %!epstopdf Gill-anal.eps

 
 
 
 % Now do the analytic solution for a line source
 
 fac3 = 3/(2*factorial(3)); fac4 = 5/(4*factorial(4)); fac5 = 3*7/(4*factorial(5)); 
 for i = 1:imax
 for j = 1:jmax
     x = xx(i,j); y = yy(i,j); 
     % First calculate the cylinder functions
     yc = y;
     poly = (1 - yc + 0.25*yc^2 - fac3*yc^3 + fac4*yc^4 - fac5*yc^5);
     expy = exp(-yc^2/4);
     dpoly = (-1 + 0.5*yc - 3*fac3*yc^2 + 4*fac4*yc^3 - 5*fac5*yc^4);
     dexpy = -0.5*yc*expy;
     exx(i,j) = expy;
     polyy(i,j) = poly;
     U(i,j)  =  expy*poly;
     Ud(i,j) =  expy*dpoly + dexpy*poly; 
     yc = -y;
     poly = (1 - yc + 0.25*yc^2 - fac3*yc^3 + fac4*yc^4 - fac5*yc^5);
     expy = exp(-yc^2/4);
     dpoly = (-1 + 0.5*yc - 3*fac3*yc^2 + 4*fac4*yc^3 - 5*fac5*yc^4);
     dexpy = -0.5*yc*expy;
     Um(i,j)  =  expy*poly;
     Umd(i,j) =  expy*dpoly + dexpy*poly; 
     yc = 1;
     poly = (1 - yc + 0.25*yc^2 - fac3*yc^3 + fac4*yc^4 - fac5*yc^5);
     expy = exp(-yc^2/4);
     dpoly = (-1 + 0.5*yc - 3*fac3*yc^2 + 4*fac4*yc^3 - 5*fac5*yc^4);
     dexpy = -0.5*yc*expy;
     U1 =  expy*poly;
     U1d =  expy*dpoly + dexpy*poly; 
     yc = -1;
     poly = (1 - yc + 0.25*yc^2 - fac3*yc^3 + fac4*yc^4 - fac5*yc^5);
     expy = exp(-yc^2/4);
     dpoly = (-1 + 0.5*yc - 3*fac3*yc^2 + 4*fac4*yc^3 - 5*fac5*yc^4);
     dexpy = -0.5*yc*expy;
     Um1 =  expy*poly;
     Um1d =  expy*dpoly + dexpy*poly; 
     yc = 0;
     poly = (1 - yc + 0.25*yc^2 - fac3*yc^3 + fac4*yc^4 - fac5*yc^5);
     expy = exp(-yc^2/4);
     dpoly = (-1 + 0.5*yc - 3*fac3*yc^2 + 4*fac4*yc^3 - 5*fac5*yc^4);
     dexpy = -0.5*yc*expy;
     U0 =  expy*poly;
     U0d =  expy*dpoly + dexpy*poly; 
     
     % Now solutions for y0 = 1;
     if y > 1
         v1(i,j) = alpha*U(i,j)/U1d;
         u1(i,j) = (y/2.)*U(i,j)/U1d;
         p1(i,j) = -Ud(i,j)/U1d;
     elseif y <= 1
         v1(i,j) = -alpha*Um(i,j)/Um1d;
         u1(i,j) = -(y/2.)*Um(i,j)/Um1d;
         p1(i,j) = -Umd(i,j)/Um1d;
     end
          % Now solutions for y0 =0;
     if y > 0
         v0(i,j) = alpha*U(i,j)/U0d;
         u0(i,j) = (y/2.)*U(i,j)/U0d;
         p0(i,j) = -Ud(i,j)/U0d;
     elseif y <= 0
         v0(i,j) = -alpha*Um(i,j)/U0d;
         u0(i,j) = -(y/2.)*Um(i,j)/U0d;
         p0(i,j) = -Umd(i,j)/U0d;
     end
     %
 end
 end
         
 % Now decimate the field for quiver plots
incr=10;
ii=0; jj=0;
for i = 1:incr:imax
    ii = ii+1;
    jj=0;
    for j=1:incr:jmax
        jj = jj+1;  
        xxd(ii,jj) = xx(i,j);
        yyd(ii,jj) = yy(i,j);
        u1d(ii,jj) = u1(i,j);
        v1d(ii,jj) = v1(i,j);
        u0d(ii,jj) = u0(i,j);
        v0d(ii,jj) = v0(i,j);
    end
end

figure(4)
hold off
contour(xx,yy,p1,5,'k')
hold on
quiver(xxd,yyd,u0d,v0d,0)
        
return

