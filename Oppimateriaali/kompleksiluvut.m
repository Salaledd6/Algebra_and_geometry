clear
%suorakulmainen muoto->kulmamuoto
z=3+4*i 
r=abs(z) %itseisarvo
phi=angle(z) %suuntakulma radiaaneina
%tai
x=real(z) %reaaliosa
y=imag(z) %imaginaariosa
[phi,r]=cart2pol(x,y)%cartesian to polar

%% kulmamuoto->suorakulmainen muoto
r=3
phi=125*pi/180
x=r*cos(phi)
y=r*sin(phi)
%tai
[x,y]=pol2cart(phi,r)%polar to cartesian
z=x+y*i

%% laskutoimitukset
clear
x=5
y=2
z=x+y*i
a=2
b=4
w=a+b*i
zpw=z+w
zmw=z-w
zkw=z*w
zjw=z/w
%%
subplot(221)
plot([0,x],[0,y],'r','linewidth',2)
hold
plot([0,a],[0,b],'b','linewidth',2)
plot([0,real(zpw)],[0,imag(zpw)],'g','linewidth',2)
hold off
grid
axis equal
legend('z','w','z+w')

subplot(222)
plot([0,x],[0,y],'r','linewidth',2)
hold
plot([0,a],[0,b],'b','linewidth',2)
plot([0,real(zmw)],[0,imag(zmw)],'g','linewidth',2)
hold off
grid
axis equal
legend('z','w','z-w')

subplot(223)
plot([0,x],[0,y],'r','linewidth',2)
hold
plot([0,a],[0,b],'b','linewidth',2)
plot([0,real(zkw)],[0,imag(zkw)],'g','linewidth',2)
hold off
grid
axis equal
legend('z','w','z*w')
   
subplot(224)
plot([0,x],[0,y],'r','linewidth',2)
hold
plot([0,a],[0,b],'b','linewidth',2)
plot([0,real(zjw)],[0,imag(zjw)],'g','linewidth',2)
hold off
grid
axis equal
legend('z','w','z/w')

%% sinikayrien yhteenlasku
clear
close all
A1=10
phi1=0
A2=4
phi2=60*pi/180

[x1,y1]=pol2cart(phi1,A1)
Y1=x1+y1*i
[x2,y2]=pol2cart(phi2,A2)
Y2=x2+y2*i

Y=Y1+Y2
A=abs(Y)
phi=angle(Y)
x=real(Y)
y=imag(Y)

figure(1)
plot([0,x1],[0,y1],'b','linewidth',2)
hold
plot([0,x2],[0,y2],'g','linewidth',2)
plot([0,x],[0,y],'r','linewidth',2)
plot([x1,x],[y1,y],'g','linewidth',1)
plot([x2,x],[y2,y],'b','linewidth',1)
hold off
grid
axis equal
legend('Y_1','Y_2','Y_1+Y_2')

%sinikayrat
w=20*2*pi
T=2*pi/w
t=0:T/100:4*T;
y1=A1*sin(w*t+phi1);
y2=A2*sin(w*t+phi2);
y=A*sin(w*t+phi);

figure(2)
plot(t,y1,'b','linewidth',1.0)
hold
plot(t,y2,'g','linewidth',1.0)
plot(t,y,'r','linewidth',1.0)
hold off
grid
legend('y_1','y_2','y_1+y_2')
xlabel('t')