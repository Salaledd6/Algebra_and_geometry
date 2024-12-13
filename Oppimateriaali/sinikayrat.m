%% sinikayra y=A*sin(w*t+phi)
clear
A=2%amplitudi
w=100*pi %kulmataajuus 
phi=130*pi/180 %vaihekulma


T=2*pi/w %jakso
%kuvaaja t=0...4T
t=0:T/100:4*T; 
y=A*sin(w*t+phi);

t0=0.3*T %tarkasteluhetki
y0=A*sin(w*t0+phi) %y-koordinaatti
x0=A*cos(w*t0+phi) %x-koordinaatti



figure(1)
plot(t,y,'b') 
hold
plot(t0,y0,'r.','markersize',20) 
plot([t0,t0],[0,y0],'r') 
plot([0,4*T],[0,0],'k') 
hold off
grid
xlabel('aika t')
ylabel('y = A sin(wt+\phi)')
ylim([-1.2*A,1.2*A])

figure(2)
plot([0,x0],[0,y0],'b','linewidth',2) %keppi
hold 
plot(x0,y0,'r.','markersize',20) 
kk=0:1:360; 
plot(A*cosd(kk),A*sind(kk),'k') 
plot([x0,x0],[0,y0],'r') 
plot([-1.2*A,1.2*A],[0,0],'k') 
plot([0,0],[-1.2*A,1.2*A],'k') 
hold off
axis([-1.2*A,1.2*A,-1.2*A,1.2*A]) 
axis square  
grid

%% s.9-10, RL-piiri
clear
R=0.5
L=0.005
U=5
w=100*pi
T=2*pi/w
t=0:T/100:4*T;
u=U*sin(w*t);
I=U/sqrt(R^2+(w*L)^2)
phi=atan(w*L/R)
i=I*sin(w*t-phi);

figure(1)
plot(t,u,'r','linewidth',1.5)
hold 
plot(t,i,'b','linewidth',1.5)
hold off 
grid
legend({'u','i'},'fontsize',14)

figure(2)
r=max([U,I])+1
plot([-r,r],[0,0],'k')
hold
plot([0,0],[-r,r],'k')
p1=plot([0,U],[0,0],'r','linewidth',3)
p2=plot([0,I*cos(-phi)],[0,I*sin(-phi)],'b','linewidth',3)
hold off
axis([-r,r,-r,r])
axis square
grid
legend([p1,p2],{'u','i'},'fontsize',14)

%% s.17-18, sinikayrien  y1=A1sin(wt) ja y2=A2sin(wt+phi) summa
clear
w=40*2*pi
A1=10
A2=4
phi=160*pi/180
%y1=A*sin(w*t)
%y2=A2*sin(w*t+phi)
%y=y1+y2=A*sin(w*t+theta)
A=sqrt((A1+A2*cos(phi))^2+(A2*sin(phi))^2)
theta=atan2(A2*sin(phi),A1+A2*cos(phi))
T=2*pi/w
t=0:T/100:4*T;
y1=A1*sin(w*t);
y2=A2*sin(w*t+phi);
y=A*sin(w*t+theta);

figure(1)
plot(t,y1,'b',t,y2,'g',t,y,'r','linewidth',1.2)
grid
legend({'y_1','y_2','y_1+y_2'},'fontsize',12)
title(['A_1 = ',num2str(A1),', A_2 = ',num2str(A2),', \phi = ',num2str(phi),...
       ', A = ',num2str(A),', \theta = ',num2str(theta)],'fontsize',12)


A1x=A1
A1y=0
A2x=A2*cos(phi)
A2y=A2*sin(phi)
Ax=A*cos(theta)
Ay=A*sin(theta)
figure(2)
plot([0,A1x],[0,0],'b','linewidth',2)
hold on
plot([0,A2x],[0,A2y],'g','linewidth',2)
plot([0,Ax],[0,Ay],'r','linewidth',2)
plot([A2x,A2x+A1x],[A2y,A2y+A1y],'b','linewidth',1)
plot([A1x,A1x+A2x],[A1y,A1y+A2y],'g','linewidth',1)
hold off
grid
axis equal
legend({'y_1','y_2','y_1+y_2'},'fontsize',12)
