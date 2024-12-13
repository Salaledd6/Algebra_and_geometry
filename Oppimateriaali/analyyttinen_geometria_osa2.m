% s.11-12, pisteiden [x1,y1] ja [x2,y2]
% kautta kulkeva suora y=kx+b
clear
close all
x1=-2
y1=3
x2=4
y2=1
k=(y2-y1)/(x2-x1)
b=-k*x1+y1
%vasen piste (sovitaan, etta x1<x2)
xv=x1-2
yv=k*xv+b
%oikea piste
xo=x2+2
yo=k*xo+b

plot([xv,xo],[yv,yo],'b','linewidth',2)
hold on
plot([x1,x2],[y1,y2],'r.','markersize',20)
hold off
grid
title(['k = ',num2str(k),', b = ',num2str(b)])

%% s.19, paraabeli 
clear
a=1.5
b=1
c=6
%huippu
x0=-b/(2*a)
y0=a*x0^2+b*x0+c
%x:n arvot kuvaajaa varten
L=2 %paraabelin leveys = 2L
xmin=x0-L
xmax=x0+L
x=xmin:L/1000:xmax; %alku:askel:loppu
%y:n arvot
y=a*x.^2+b*x+c;
plot(x,y,'b','linewidth',2)
hold
plot(x0,y0,'r.','markersize',20)
hold off
grid
axis equal

%% s.25, kolmen pisteen kautta kulkeva paraabeli

clear
close all
x1=-1
y1=2
x2=1
y2=0
x3=2
y3=4


% tapa 2

A=x1^2-x2^2
B=x1-x2
C=x1^2-x3^2
D=x1-x3
E=y1-y2
F=y1-y3

a=(D*E-B*F)/(A*D-B*C)
b=(A*F-C*E)/(A*D-B*C)
c=y1-a*x1^2-b*x1
%%
%kuvaaja
xmin=min([x1,x2,x3])-1
xmax=max([x1,x2,x3])+1
x=xmin:(xmax-xmin)/100:xmax; %alku:askel:loppu
%y:n arvot
y=a*x.^2+b*x+c;

plot(x,y,'linewidth',1.5)
hold
plot([x1,x2,x3],[y1,y2,y3],'r.','markersize',20)
hold off
grid

%% s.27 paraabelin ja suoran leikkauspisteet
clear
%suora y=kx+p
k=1
p=3
%paraabeli y=ax^2+bx+c
a=0.5
b=-2
c=2

D=(b-k)^2-4*a*(c-p)

if D>=0 %leikkaavat
    x1=(-(b-k)-sqrt(D))/(2*a)
    y1=k*x1+p
    x2=(-(b-k)+sqrt(D))/(2*a)
    y2=k*x2+p
    xmin=min([x1,x2])-2
    xmax=max([x1,x2])+2
    x=xmin:(xmax-xmin)/100:xmax;
    y=a*x.^2+b*x+c;
    ymin=k*xmin+p
    ymax=k*xmax+p
else
    x0=-b/(2*a);
    xmin=x0-3
    xmax=x0+3
    x=xmin:(xmax-xmin)/100:xmax;
    y=a*x.^2+b*x+c;
    ymin=k*xmin+p
    ymax=k*xmax+p
end

plot(x,y,'linewidth',1.5)
hold
plot([xmin,xmax],[ymin,ymax],'linewidth',1.5)
if D>=0
    plot([x1,x2],[y1,y2],'r.','markersize',20)
end
hold off
grid
title(['k = ',num2str(k),', p = ',num2str(p),', a = ',num2str(a),...
       ', b = ',num2str(b),', c = ',num2str(c)])


%% s.29, polttopiste ja heijastusominaisuus

clear
%close all
a=0.5%y=ax^2
L=2
x=-L:L/100:L;
y=a*x.^2;
%polttopiste
Fx=0
Fy=1/(4*a)
%heijastuspiste
xt=0.2
yt=a*xt^2
%tangentin kulmakerroin (s.31-32)
k=2*a*xt
plot(x,y,'linewidth',2)
hold on
plot([xt,xt,Fx],[yt+2,yt,Fy],'g','linewidth',2)%sade
plot([xt-1,xt+1],[yt-k*1,yt+k*1],'k','linewidth',2)%tangentti
plot(Fx,Fy,'r.','markersize',20)%polttopiste
plot(xt,yt,'b.','markersize',20)%heijastumispiste
hold off
axis equal
grid

%% s.35 tasaisesti kiihtyva liike
%s=v0*t+1/2*a*t^2

clear
v0=5
L=30
amin=-v0^2/(2*L)
a=amin+1
a=6.1
%L=v0*T+1/2*a*T^2
if a>0
    T=(-v0+sqrt(v0^2+2*a*L))/a;
elseif a>=amin
    T1=(-v0+sqrt(v0^2+2*a*L))/a;
    T2=(-v0-sqrt(v0^2+2*a*L))/a;
    T=T2;
else
     T=-2*v0/a;
end

t=0:T/100:T;
s=v0*t+1/2*a*t.^2;

plot(t,s,'linewidth',1.5)
hold
plot([0,T],[L,L],'k','linewidth',1.5)
if a>0
    plot(T,L,'r.','markersize',20)
elseif a>=amin
    plot([T1,T2],[L,L],'r.','markersize',20)
end
hold off
grid
xlim([0,T])
ylim([0,max([1.2*L,max(s)])])
if a>0
title(['L = ',num2str(L),', v_0 = ',num2str(v0),', a = ',num2str(a),', amin = ',num2str(amin),', T = ',num2str(T)])
elseif a>=amin
title(['L = ',num2str(L),', v_0 = ',num2str(v0),', a = ',num2str(a),', amin = ',num2str(amin),', T_1 = ',num2str(T1),',  T_2 = ',num2str(T2)])
else
title(['L = ',num2str(L),', v_0 = ',num2str(v0),', a = ',num2str(a),', amin = ',num2str(amin)])
end
xlabel('aika t')
ylabel('etaisyys s')

%% s.45 vaakasuora lentomatka
clear
h=5
v0=25
g=9.81
alfa=0:0.1:90;
xmax=v0^2/g*cosd(alfa).*(sind(alfa)+sqrt(sind(alfa).^2+2*g*h/v0^2));

%suurimmillaan, kun alfa=alfa0
alfa0=atand(v0/sqrt(v0^2+2*g*h))
%maksimiarvo
xmax0=v0/g*sqrt(v0^2+2*g*h)

figure(1)
plot(alfa,xmax,'linewidth',1.5)
hold
plot(alfa0,xmax0,'r.','markersize',20)
hold off
grid
xlabel('lahtokulma \alpha')
ylabel('vaakasora lentomatka')
title(['xmax = ',num2str(xmax0),', kun \alpha = ',num2str(alfa0)])

%% s.51, maali
clear
L=10
H=6
g=9.81
minv0=sqrt(g*L^2/(sqrt(H^2+L^2)-H))
%%
v0=minv0+1
if v0>=minv0
k=g*L^2/v0^2+H
C=sqrt(L^2+H^2)
beta=asind(k/C)
theta=atan2d(H,-L)
alfa1=1/2*(180+beta-theta)
alfa2=1/2*(360-beta-theta)

a1=-g/(2*(v0*cosd(alfa1))^2)
b1=tand(alfa1)
a2=-g/(2*(v0*cosd(alfa2))^2)
b2=tand(alfa2)

x=0:L/100:L;
y1=a1*x.^2+b1*x;
y2=a2*x.^2+b2*x;

plot(x,y1,'linewidth',1.5)
hold on
plot(x,y2,'linewidth',1.5)
plot(L,H,'r.','markersize',20)
hold off
grid
axis equal
title(['v0 = ',num2str(v0),', min v0 = ',num2str(minv0),...
    ', \alpha_1 = ',num2str(alfa1),', \alpha_2 = ',num2str(alfa2)])
end
%% s.59, ympyran ja suoran leikkauspiste
clear
close all %sulkee kuvaikkunat
k=0.5
b=1.5
x0=5
y0=7
r=5

%tapa2
A=k^2+1
B=2*b*k-2*k*y0-2*x0
C=b^2-2*b*y0+x0^2+y0^2-r^2

D=B^2-4*A*C

%ympyran pisteet
t=0:360;
xymp=x0+r*cosd(t);
yymp=y0+r*sind(t);

plot(xymp,yymp,'b','linewidth',2)
hold on
plot(x0,y0,'b.','markersize',20)
xv=x0-r-1 %suoran vasen piste kuvaa varten
xo=x0+r+1 %suoran oikea piste kuvaa varten
yv=k*xv+b
yo=k*xo+b

plot([xv,xo],[yv,yo],'k','linewidth',2)
if D>=0 %leikkaavat
x1=(-B-sqrt(D))/(2*A)
y1=k*x1+b
x2=(-B+sqrt(D))/(2*A)
y2=k*x2+b
plot([x1,x2],[y1,y2],'r.','markersize',20)
end
hold off
axis equal
grid



%% s.61 paikannus
clear
%close all
%tukiasemat
%P1
x1=0
y1=0
%P2
x2=6
y2=-2
%P3
x3=3
y3=3
%testipiste P
x0=3
y0=1
%mittaustulokset

r1=sqrt((x0-x1)^2+(y0-y1)^2)-0.7 %PP1+virhe
r2=sqrt((x0-x2)^2+(y0-y2)^2)-0.4 %PP2+virhe
r3=sqrt((x0-x3)^2+(y0-y3)^2)-0.9 %PP3+virhe

%arvio P:n paikalle
%yhtaloparin
%ax+by=e
%cx+dy=f
%ratkaisu x,y

a=2*(x3-x1)
b=2*(y3-y1)
e=x3^2-x1^2+y3^2-y1^2+r1^2-r3^2
c=2*(x3-x2)
d=2*(y3-y2)
f=x3^2-x2^2+y3^2-y2^2+r2^2-r3^2

x=(d*e-b*f)/(a*d-b*c)
y=(a*f-c*e)/(a*d-b*c)

%ympyrat
t=0:1:360;
co=cosd(t);
si=sind(t);

ymp1x=x1+r1*co;
ymp1y=y1+r1*si;
ymp2x=x2+r2*co;
ymp2y=y2+r2*si;
ymp3x=x3+r3*co;
ymp3y=y3+r3*si;

%ympyraparien maaraamat suorat
%ympyrat 1 ja 2
L=sqrt((x1-x2)^2+(y1-y2)^2)%keskipisteiden valinen etaisyys
L1=L/2+(r1^2-r2^2)/(2*L) %suoran etaisyys keskipisteesta [x1,y1]
th=atan2d(y2-y1,x2-x1) %keskipisteiden yhdysjanan suuntakulma
px=x1+L1*cosd(th) %suoran piste yhdysjanalla
py=y1+L1*sind(th)
R=max([r1,r2]) %suoran pituus  = 2R
p12ax=px+R*cosd(th+90); %suoran paatepisteet
p12ay=py+R*sind(th+90);
p12bx=px+R*cosd(th-90);
p12by=py+R*sind(th-90);
%ympyrat 1 ja 3
L=sqrt((x1-x3)^2+(y1-y3)^2)
L1=L/2+(r1^2-r3^2)/(2*L)
th=atan2d(y3-y1,x3-x1)
px=x1+L1*cosd(th)
py=y1+L1*sind(th)
R=max([r1,r3])
p13ax=px+R*cosd(th+90);
p13ay=py+R*sind(th+90);
p13bx=px+R*cosd(th-90);
p13by=py+R*sind(th-90);
%ympyrat 2 ja 3
L=sqrt((x2-x3)^2+(y2-y3)^2)
L1=L/2+(r2^2-r3^2)/(2*L)
th=atan2d(y3-y2,x3-x2)
px=x2+L1*cosd(th)
py=y2+L1*sind(th)
R=max([r2,r3])
p23ax=px+R*cosd(th+90);
p23ay=py+R*sind(th+90);
p23bx=px+R*cosd(th-90);
p23by=py+R*sind(th-90);



plot(ymp1x,ymp1y,'r','linewidth',2)
hold
plot(ymp2x,ymp2y,'g','linewidth',2)
plot(ymp3x,ymp3y,'b','linewidth',2)
plot(x1,y1,'r.','markersize',20)
plot(x2,y2,'g.','markersize',20)
plot(x3,y3,'b.','markersize',20)
p1=plot(x,y,'k.','markersize',20) %arvio
p2=plot(x0,y0,'m.','markersize',20) %todellinen
plot([p12ax,p12bx],[p12ay,p12by],'k')
plot([p13ax,p13bx],[p13ay,p13by],'k')
plot([p23ax,p23bx],[p23ay,p23by],'k')
hold off
grid
axis equal
legend([p1,p2],{'arvio','P'},'fontsize',12)

%% s.74, ellipsin heijastusominaisuus
clear
a=3
b=2
F1=sqrt(a^2-b^2)
F2=-F1
%parametrimuoto (s.71-72)
t=0:1:360;
x=a*cosd(t);
y=b*sind(t);

%ellipsin piste
k=75
Px=x(k)
Py=y(k)

%PF1+PF2=2a ?
PF1=sqrt((Px-F1)^2+(Py-0)^2)
PF2=sqrt((Px-F2)^2+(Py-0)^2)
['PF1+PF2 = ',num2str(PF1+PF2)]
['2a = ',num2str(2*a)]
%%
%ellipsin tangentti pisteessa P
%kulmakerroin (ellipsi_ja_hyperbeli.pdf)
k=-Px*Py/(a^2-Px^2)
L=2 %tangentin leveys = 2L
xv=Px-L
yv=Py-k*L
xo=Px+L
yo=Py+k*L



plot(x,y,'linewidth',1.5)
hold on
plot([-(a+1),a+1],[0,0],'k')
plot([0,0],[-(a+1),a+1],'k')
plot([F1,F2],[0,0],'r.','markersize',20)
plot([F1,Px,F2],[0,Py,0],'r','linewidth',1.3)
plot([xv,xo],[yv,yo],'k','linewidth',1.3)
plot(Px,Py,'b.','markersize',20)
hold off
grid
axis([-(a+1),a+1,-(a+1),a+1])
axis square



%% s.77- hyperbeli
clear
close
a=5
b=2
t=-2:0.01:2;
x=a*cosh(t);
y=b*sinh(t);

%polttopisteet
F1x=-sqrt(a^2+b^2)
F1y=0
F2x=sqrt(a^2+b^2)
F2y=0

%hyperbelin piste

k=100
Px=x(k);
Py=y(k);

%tarkastus
%PF1-PF2=2a ? (s.89)
PF1=sqrt((Px-F1x)^2+(Py-F1y)^2)
PF2=sqrt((Px-F2x)^2+(Py-F2y)^2)
PF1-PF2
2*a


plot(x,y,'b','linewidth',2) %oikea haara
hold
plot(-x,y,'b','linewidth',2) %vasen haara
r=max(x)
plot([-r,r],[0,0],'k')
plot([0,0],[-r,r],'k')
p1=plot([0,a],[0,0],'g','linewidth',3)
p2=plot([0,0],[0,b],'m','linewidth',3)
plot([a,a,-a,-a,a],[-b,b,b,-b,-b],'k')
plot([-r,r],[-b/a*r,b/a*r],'c')
plot([-r,r],[b/a*r,-b/a*r],'c')
p3=plot([F1x,F2x],[F1y,F2y],'r.','markersize',20)
p4=plot(Px,Py,'k.','markersize',20)
hold off
grid
axis([-r,r,-r,r])
axis square
title(['a = ',num2str(a),', b = ',num2str(b),', F_1,F_2 = \pm ',num2str(F1x),', P = [',num2str(Px),',',num2str(Py),']'])
legend([p1,p2,p3,p4],{'a','b','F_1,F_2','P'},'fontsize',11)

%% s.81- hyperboliset kosini ja sini ja kaanteiset area-funktiot
clear
t=-3:0.01:3;
x=cosh(t);
y=sinh(t);
c=1:0.01:10;
tc=acosh(c);
s=-10:0.01:10;
ts=asinh(s);

figure(1)
subplot(121)
plot(t,x,'linewidth',1.5)
grid
xlabel('t')
ylabel('cosh(t)')
xlim([-3,3])
ylim([0,10])

subplot(122)
plot(c,tc,'linewidth',1.5)
grid
xlabel('c')
ylabel('acosh(c)')
ylim([0,3])
xlim([0,10])


figure(2)
subplot(121)
plot(t,y,'linewidth',1.5)
grid
xlabel('t')
ylabel('sinh(t)')
xlim([-3,3])
ylim([-10,10])


subplot(122)
plot(s,ts,'linewidth',1.5)
grid
xlabel('s')
ylabel('asinh(s)')
ylim([-3,3])
xlim([-10,10])

%% s.91
clear
%close
L=10 %polttopisteiden F1 ja F2 valinen etaisyys F1F2
%polttopisteet
F1x=-L/2
F1y=0
F2x=L/2
F2y=0
d=-7%PF1-PF2, pitaa olla valilla -L ... L

%puoliakselit
a=abs(d)/2 %abs = itseisarvo, absolute value
b=sqrt((L/2)^2-a^2)

%hyperbelin pisteet
t=-2:0.01:2;
x0=(F1x+F2x)/2 %keskipiste
if d>0
xa=x0+a*cosh(t); %oikea haara
xb=x0-a*cosh(t); %vasen haara
else
xa=x0-a*cosh(t); %vasen haara
xb=x0+a*cosh(t); %oikea haara
end
y=b*sinh(t);

%P
k=300
Px=xa(k)
Py=y(k)

%tarkastus, PF1-PF2=d ?
PF1=sqrt((Px-F1x)^2+(Py-F1y)^2)
PF2=sqrt((Px-F2x)^2+(Py-F2y)^2)
PF1-PF2
d

plot(F1x,F1y,'m.','markersize',20)
hold
plot(F2x,F2y,'g.','markersize',20)
plot(Px,Py,'k.','markersize',20)
p1=plot([F1x,Px],[F1y,Py],'m','linewidth',1.5)
p2=plot([F2x,Px],[F2y,Py],'g','linewidth',1.5)
p3=plot(xa,y,'b','linewidth',2)
p4=plot(xb,y,'b')

r=12
plot([x0-r,x0+r],[0,0],'k')
%plot([0,0],[-r,r],'k')
hold off
grid
axis([x0-r,x0+r,-r,r])
axis square
title({['L = ',num2str(L),', d = ',num2str(d),...
      ', a = ',num2str(a),', b = ',num2str(b)]},'fontsize',12)
legend([p1,p2,p3,p4],{'PF_1','PF_2','PF_1 - PF_2 = d','PF_1 - PF_2 = - d'},'fontsize',10)

%%  s.93,  paikannus
clear
%close all
%F1=[0,0]
x1=0;
y1=0;
%F2=[x2,0]
x2=6;
y2=0;
%F3=[x3,y3]
x3=4;
y3=4;

L12=sqrt((x1-x2)^2+(y1-y2)^2); %F1F2
L13=sqrt((x1-x3)^2+(y1-y3)^2); %F1F3

%etaisyyksien erotukset
%joko
d12=4; %PF1-PF2, valilla -L12 ... L12
d13=2; %PF1-PF3, valilla -L13 ... L13
%tai 
Px=4;
Py=-2;
PF1=sqrt((x1-Px)^2+(y1-Py)^2);
PF2=sqrt((x2-Px)^2+(y2-Py)^2);
PF3=sqrt((x3-Px)^2+(y3-Py)^2);
d12=PF1-PF2;
d13=PF1-PF3;


%hyperbelien leikkauspiste
%merkinnat kuten hyperbeli_paikannus.pdf:ssa

%a*x+b*y=e1*R+e2
%c*x+d*y=f1*R+f2
a=2*x2;
b=2*y2;
c=2*x3;
d=2*y3;
e1=2*d12;
e2=x2^2+y2^2-d12^2;
f1=2*d13;
f2=x3^2+y3^2-d13^2;
%x=alfa*R+beta
%y=gamma*R+delta
alfa=(d*e1-b*f1)/(a*d-b*c);
beta=(d*e2-b*f2)/(a*d-b*c);
gamma=(a*f1-c*e1)/(a*d-b*c);
delta=(a*f2-c*e2)/(a*d-b*c);
%A*R^2+B*R+C=0
A=alfa^2+gamma^2-1;
B=2*(alfa*beta+gamma*delta);
C=beta^2+delta^2;
D=B^2-4*A*C;

if D>=0
Rp=(-B+sqrt(D))/(2*A);
Rm=(-B-sqrt(D))/(2*A);
xp=alfa*Rp+beta;
yp=gamma*Rp+delta;
xm=alfa*Rm+beta;
ym=gamma*Rm+delta;
%tarkastus
PpF1=sqrt((x1-xp)^2+(y1-yp)^2);
PpF2=sqrt((x2-xp)^2+(y2-yp)^2);
PpF3=sqrt((x3-xp)^2+(y3-yp)^2);
PpF1-PpF2
d12
PpF1-PpF3
d13
PmF1=sqrt((x1-xm)^2+(y1-ym)^2);
PmF2=sqrt((x2-xm)^2+(y2-ym)^2);
PmF3=sqrt((x3-xm)^2+(y3-ym)^2);
PmF1-PmF2
d12
PmF1-PmF3
d13
end



%kuva
%hyperbelit
tmax=2;
t=-tmax:tmax/100:tmax;
%PF1-PF2=d12
%polttopisteiden F1 ja F2 valinen etaisyys = L12
%keskipiste
x01=L12/2; %y01=0
%puoliakselit
a1=abs(d12/2);
b1=sqrt((L12/2)^2-a1^2);

if d12>0
h1x=x01+a1*cosh(t); %oikea haara
else
h1x=x01-a1*cosh(t); %vasen haara
end
h1y=b1*sinh(t);


%PF1-PF3=d13
%polttopisteiden valinen etaisyys = L13
%hyperbelin pisteiden koordinaatit, kun polttopisteet ovat [0,0] ja [L13,0]
%keskipiste
x02=L13/2; %y02=0
%puoliakselit
a2=abs(d13/2);
b2=sqrt((L13/2)^2-a2^2);
if d13>0
h2x=x02+a2*cosh(t); %oikea haara
else
h2x=x02-a2*cosh(t); %vasen haara
end
h2y=b2*sinh(t);

%napakoordinaatit
r=sqrt(h2x.^2+h2y.^2);
th=atan2(h2y,h2x);
%kierretaan pisteita F1F3:n suuntakulman alfa verran
alfa=atan2(y3,x3);
%kierretty hyperbeli
h2x=r.*cos(th+alfa);
h2y=r.*sin(th+alfa);


plot(x1,y1,'r.','markersize',20)
hold on
plot(x2,y2,'g.','markersize',20)
plot(x3,y3,'m.','markersize',20)
plot(h1x,h1y,'b','linewidth',1.5) %PF1-PF2=d12
plot(h2x,h2y,'c','linewidth',1.5) %PF1-PF3=d13
if D>=0
     %[xm,ym]
     P1F1=sqrt((xm-x1)^2+(ym-y1)^2);
     P1F2=sqrt((xm-x2)^2+(ym-y2)^2);
     P1F3=sqrt((xm-x3)^2+(ym-y3)^2);
     if abs(P1F1-P1F2-d12)<10^-10 & abs(P1F1-P1F3-d13)<10^-10
     plot(xm,ym,'k.','markersize',20)
     end
     %[xp,yp]
     P2F1=sqrt((xp-x1)^2+(yp-y1)^2);
     P2F2=sqrt((xp-x2)^2+(yp-y2)^2);
     P2F3=sqrt((xp-x3)^2+(yp-y3)^2);
     if abs(P2F1-P2F2-d12)<10^-10 & abs(P2F1-P2F3-d13)<10^-10
     plot(xp,yp,'k.','markersize',20)
     end
end
hold off
grid
axis equal
title(['d_{12} = ',num2str(d12),', d_{13} = ',num2str(d13)])
legend('F_1','F_2','F_3','PF_1-PF_2=d_{12}','PF_1-PF_3=d_{13}','P')

