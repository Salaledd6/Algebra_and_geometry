%% s.5 mediaanien leikkauspiste
clear
Ax=1
Ay=1
Bx=6
By=2
Cx=4
Cy=5

%sivujen keskipisteet
Dx=1/2*(Bx+Cx)
Dy=1/2*(By+Cy)
Ex=1/2*(Ax+Cx)
Ey=1/2*(Ay+Cy)
Fx=1/2*(Ax+Bx)
Fy=1/2*(Ay+By)

%mediaanien leikkauspiste, kolmion ABC painopiste
Px=1/3*(Ax+Bx+Cx)
Py=1/3*(Ay+By+Cy)

AP=sqrt((Px-Ax)^2+(Py-Ay)^2)
AD=sqrt((Dx-Ax)^2+(Dy-Ay)^2)
2/3*AD %=AP ?

BP=sqrt((Px-Bx)^2+(Py-By)^2)
BE=sqrt((Ex-Bx)^2+(Ey-By)^2)
2/3*BE %=BP ?

CP=sqrt((Px-Cx)^2+(Py-Cy)^2)
CF=sqrt((Fx-Cx)^2+(Fy-Cy)^2)
2/3*CF %=CP ?


%% kuva

plot([Ax,Bx,Cx,Ax],[Ay,By,Cy,Ay],'linewidth',1.5) %kolmio ABC
hold
plot([Ax,Dx],[Ay,Dy],'r') %viiva AD
plot([Bx,Ex],[By,Ey],'r') %viiva BE
plot([Cx,Fx],[Cy,Fy],'r') %viiva CF
plot(Px,Py,'k.','markersize',20) %piste P
hold off
grid 
axis equal %akseleiden mittakaavat yhtasuuriksi 

%% s.9
clear 
%close all
%napakoordinaatit
r=4 %OP
th=-130 %P:n suuntakulma
%suorakulmaiset koordinaatit
Px=r*cosd(th)
Py=r*sind(th)

%kuva
plot([-r,r],[0,0],'k','linewidth',2) %x-akseli
hold %pito paalle, jotta seuraava plot ei pyyhi edellista
plot([0,0],[-r,r],'k','linewidth',2) %y-akseli
plot([0,Px],[0,Py],'r','linewidth',2) %viiva OP
plot(Px,Py,'r.','markersize',20) %piste P
hold off %pito pois
grid %taustaristikko
%joko
axis([-r,r,-r,r]) %akseleiden rajat, [xmin,xmax,ymin,ymax]
axis square %kuva nelion muotoiseksi (mittakaavat samat, kun xmax-xmin=ymax-ymin)
%tai
axis equal
title(['r = ',num2str(r),', \theta = ',...
      num2str(th),'^o,  Px  = ',num2str(Px),', Py = ',num2str(Py)])
%num2str=number to string, nayttaa muuttujan arvon 

%% 
clear
%suorakulmaiset 
Px=-3
Py=-2
%napa
th=atan2d(Py,Px) %OP:n suuntakulma
r=sqrt(Px^2+Py^2) %OP:n pituus

plot([-r,r],[0,0],'k','linewidth',2) %x-akseli
hold 
plot([0,0],[-r,r],'k','linewidth',2) %y-akseli
plot([0,Px],[0,Py],'r','linewidth',2) %viiva OP
plot(Px,Py,'r.','markersize',20) %piste P
hold off
grid %taustaristikko
axis([-r,r,-r,r]) %akseleiden rajat, [xmin,xmax,ymin,ymax]
axis square %kuva nelion muotoiseksi
title(['Px  = ',num2str(Px),', Py = ',num2str(Py),', r = ',num2str(r),', \theta = ',num2str(th),'^o'])



%% s. 17, pisteen A kautta kulkeva suora, suuntakulma th
clear
Ax=1
Ay=2
th=40
%piste suoralla
r=3
Px=Ax+r*cosd(th)
Py=Ay+r*sind(th)

p1=plot(Ax,Ay,'r.','markersize',20)
hold
plot([Ax-r*cosd(th),Ax+2*r*cosd(th)],...
     [Ay-r*sind(th),Ay+2*r*sind(th)],'r','linewidth',1)
p2=plot(Px,Py,'b.','markersize',20)  
hold off
grid
axis equal %akseleiden mittakaavat yhtasuuret
legend([p1,p2],{'A','P'},'fontsize',12)
title(['\theta = ',num2str(th),'^o, r = ',num2str(r)])


%% s.21, suorien leikkauspiste
clear
Ax=3
Ay=2
theta=-130
Bx=3
By=1
delta=40

a=cosd(theta)
b=-cosd(delta)
c=sind(theta)
d=-sind(delta)
e=Bx-Ax
f=By-Ay
r=(d*e-b*f)/(a*d-b*c)
t=(a*f-c*e)/(a*d-b*c)

Px=Ax+r*cosd(theta)
Py=Ay+r*sind(theta)

%tarkastus
Bx+t*cosd(delta) %=Px ? 
By+t*sind(delta) %=Py ?

%%
%suoraan funktiolla suorien_leikkauspiste.m
%(tiedoston pitaa olla  oletushakemistossa (current directory)
[Px,Py,r,t]=suorien_leikkauspiste(Ax,Ay,theta,Bx,By,delta) %theta ja delta asteina
%[Qx,Qy,u,v]=suorien_leikkauspiste(Mx,My,alfa,Nx,Ny,beta)
%%
plot(Ax,Ay,'r.','markersize',20)
hold
plot(Bx,By,'g.','markersize',20)
plot(Px,Py,'k.','markersize',20)
plot([Ax-r*cosd(theta),Ax+2*r*cosd(theta)],...
     [Ay-r*sind(theta),Ay+2*r*sind(theta)],'r','linewidth',1.5)
plot([Bx-t*cosd(delta),Bx+2*t*cosd(delta)],...
     [By-t*sind(delta),By+2*t*sind(delta)],'g','linewidth',1.5) 
hold off
grid
axis equal
title(['r = ',num2str(r),', t = ',num2str(t),', \theta = ',num2str(theta),', \delta = ',num2str(delta)])
legend('A','B','P')

%% s.25 paikannus
clear
Ax=6
Ay=4
theta=-60
Bx=1
By=5
delta=-150
[Px,Py,r,t]=suorien_leikkauspiste(Ax,Ay,theta,Bx,By,delta)

plot([Px,Ax],[Py,Ay],'r')
hold on
plot([Px,Bx],[Py,By],'b')
plot([Px,Px+1],[Py,Py],'k')
p1=plot(Ax,Ay,'r.','markersize',20)
p2=plot(Bx,By,'b.','markersize',20)
p3=plot(Px,Py,'k.','markersize',20)
hold off
grid
axis equal
legend([p1,p2,p3],'A','B','P')
title(['\theta = ',num2str(theta),', \delta = ',num2str(delta)])

%% s.27 Fermatin piste
clear
%oletus: kiertosuunta A->B->C vastapaivaan
Ax=0
Ay=0
Bx=6
By=2
Cx=1
Cy=5

AB=sqrt((Bx-Ax)^2+(By-Ay)^2)
BC=sqrt((Cx-Bx)^2+(Cy-By)^2)
CA=sqrt((Ax-Cx)^2+(Ay-Cy)^2)

thAB=atan2d(By-Ay,Bx-Ax)
thBC=atan2d(Cy-By,Cx-Bx)
thCA=atan2d(Ay-Cy,Ax-Cx)

AT=AB
thAT=thAB-60

BR=BC
thBR=thBC-60

CS=CA
thCS=thCA-60

Tx=Ax+AT*cosd(thAT)
Ty=Ay+AT*sind(thAT)

Rx=Bx+BR*cosd(thBR)
Ry=By+BR*sind(thBR)

Sx=Cx+CS*cosd(thCS)
Sy=Cy+CS*sind(thCS)

thAR=atan2d(Ry-Ay,Rx-Ax)
thBS=atan2d(Sy-By,Sx-Bx)

[Px,Py,r,t]=suorien_leikkauspiste(Ax,Ay,thAR,Bx,By,thBS)



plot([Ax,Bx,Cx,Ax],[Ay,By,Cy,Ay],'linewidth',2)
hold on
plot([Ax,Tx,Bx],[Ay,Ty,By],'r','linewidth',1.5)
plot([Bx,Rx,Cx],[By,Ry,Cy],'m','linewidth',1.5)
plot([Cx,Sx,Ax],[Cy,Sy,Ay],'g','linewidth',1.5)
plot([Ax,Rx],[Ay,Ry],'k')
plot([Bx,Sx],[By,Sy],'k')
plot([Cx,Tx],[Cy,Ty],'k')
plot([Ax,Px,Bx],[Ay,Py,By],'k','linewidth',1.5)
plot([Cx,Px],[Cy,Py],'k','linewidth',1.5)
plot(Px,Py,'k.','markersize',20)
hold off
grid
axis equal

%% s.29, suoran ja kolmion leikkauspisteet
clear
%suoran pisteet
Px=1
Py=4
Qx=7
Qy=2
%kolmion pisteet
Ax=0
Ay=0
Bx=5
By=1
Cx=3
Cy=6

%suorien PQ ja AB leikkauspiste
%suuntakulmat 
thPQ=atan2d(Qy-Py,Qx-Px)
thAB=atan2d(By-Ay,Bx-Ax)
[S1x,S1y,r1,t1]=suorien_leikkauspiste(Ax,Ay,thAB,Px,Py,thPQ)

%suorien PQ ja BC leikkauspiste
%suuntakulmat 
%thPQ=atan2d(Qy-Py,Qx-Px)
thBC=atan2d(Cy-By,Cx-Bx)
[S2x,S2y,r2,t2]=suorien_leikkauspiste(Bx,By,thBC,Px,Py,thPQ)

%suorien PQ ja CA leikkauspiste
%suuntakulmat 
%thPQ=atan2d(Qy-Py,Qx-Px)
thCA=atan2d(Ay-Cy,Ax-Cx)
[S3x,S3y,r3,t3]=suorien_leikkauspiste(Cx,Cy,thCA,Px,Py,thPQ)

%%kuva
%kolmio ABC
plot([Ax,Bx,Cx,Ax],[Ay,By,Cy,Ay],'linewidth',1.5)
hold 
%suora PQ
PQ=sqrt((Qx-Px)^2+(Qy-Py)^2)
tmin=min([t1,t2,t3,0])-1
tmax=max([t1,t2,t3,PQ])+1
Rx=Px+tmin*cosd(thPQ)
Ry=Py+tmin*sind(thPQ)
Sx=Px+tmax*cosd(thPQ)
Sy=Py+tmax*sind(thPQ)
plot([Rx,Sx],[Ry,Sy],'r','linewidth',1.5)
plot([Px,Qx],[Py,Qy],'r.','markersize',20)

AB=sqrt((Bx-Ax)^2+(By-Ay)^2)
if (r1>=0) && (r1<=AB)
  plot(S1x,S1y,'g.','markersize',20)
end

BC=sqrt((Cx-Bx)^2+(Cy-By)^2)
if (r2>=0) && (r2<=BC)
  plot(S2x,S2y,'g.','markersize',20)
end

CA=sqrt((Ax-Cx)^2+(Ay-Cy)^2)
if (r3>=0) && (r3<=CA)
  plot(S3x,S3y,'g.','markersize',20)
end

hold off
grid
axis equal


%% s.31, ympyra, keskipiste P ja sade r
clear
Px=1
Py=2
r=4
%kiertokulma
th=0:1:360; %th=[0,1,2,...,360], alku:askel:loppu
Qx=Px+r*cosd(th);
Qy=Py+r*sind(th);
%kaari
ka=20%alkukulma
kl=50%loppukulma
thk=ka:(kl-ka)/100:kl;
Qkx=Px+r*cosd(thk);
Qky=Py+r*sind(thk);

plot(Qx,Qy,'b','linewidth',2)
hold 
plot(Px,Py,'b.','markersize',20)
plot(Qkx,Qky,'r','linewidth',3)
hold off
grid
axis equal

%% s.33, iso ja pieni ympyra
clear 
R=1 %ison ympyran sade 
r=(3-sqrt(8))*R %pienen ympyran sade 
%ympyroiden pisteet 
%iso 
%keskipiste 
X0=0
Y0=0
th=0:1:360; %kiertokulma
c=cosd(th);
s=sind(th);
X=X0+R*c;
Y=Y0+R*s;

%pieni
%keskipiste 
x0=X0+(R-r)
y0=Y0-(R-r)

x=x0+r*c;
y=y0+r*s;

plot([-R,R,R,-R,-R],[-R,-R,R,R,-R],'k','linewidth',1.5) %nelio 
hold
plot(X,Y,'b','linewidth',1.5)
plot(X0,Y0,'b.','markersize',20)
plot(x,y,'r','linewidth',1.5)
plot(x0,y0,'r.','markersize',20)
hold off 
grid
axis equal

%% s.35 (pisteiden A,B ja C kautta kulkeva ympyra)
clear 
Ax=1
Ay=1
Bx=5
By=2
Cx=4
Cy=4

%sivun AB keskipiste 
Fx=(Ax+Bx)/2
Fy=(Ay+By)/2
%ja suuntakulma
kAB=atan2d(By-Ay,Bx-Ax)
theta=kAB+90 %keskinormaalin suuntakulma

%sivun AC keskipiste 
Ex=(Ax+Cx)/2
Ey=(Ay+Cy)/2
%ja suuntakulma
kAC=atan2d(Cy-Ay,Cx-Ax)
delta=kAC+90 %keskinormaalin suuntakulma

%ympyran keskipiste eli keskinormaalien D,theta ja E,delta leikkauspiste P

[Px,Py,r,t]=suorien_leikkauspiste(Fx,Fy,theta,Ex,Ey,delta)

%ympyran sade
R=sqrt((Px-Ax)^2+(Py-Ay)^2)
%%
%suoraan funktiolla ympyraABC.m 
%tiedostojen ympyraABC.m ja suorien_leikkauspiste.m
%pitaa olla oletushakemistossa
[Px,Py,R]=ympyraABC(Ax,Ay,Bx,By,Cx,Cy)
%[Qx,Qy,r]=ympyraABC(Kx,Ky,Lx,Ly,Mx,My)
%%
%ympyran pisteet
th=0:1:360;
ympx=Px+R*cosd(th);
ympy=Py+R*sind(th);



plot(ympx,ympy,'b','linewidth',2)
hold
plot(Px,Py,'b.','markersize',20)
plot([Ax,Bx,Cx],[Ay,By,Cy],'r.','markersize',20)
hold off
grid
axis equal

%% s. 39 (tangentit ympyralle)
clear
Ax=5
Ay=4
Px=5
Py=1
r=2

APx=Px-Ax
APy=Py-Ay
AP=sqrt(APx^2+APy^2)
d=asin(r/AP)
AS=sqrt(AP^2-r^2)
AT=AS
th=atan2(APy,APx)
Sx=Ax+AS*cos(th+d)
Sy=Ay+AS*sin(th+d)
Tx=Ax+AT*cos(th-d)
Ty=Ay+AT*sin(th-d)
%%
%%ympyran pisteet
th=0:1:360; 
x=Px+r*cosd(th);
y=Py+r*sind(th);

plot(x,y,'k','linewidth',2)
hold
plot([Ax,Px],[Ay,Py],'c')
plot([Sx,Ax,Tx],[Sy,Ay,Ty],'r.-','linewidth',2,'markersize',20)
plot([Sx,Px,Tx],[Sy,Py,Ty],'g','linewidth',2)
plot(Px,Py,'k.','markersize',20)
hold off
grid
axis equal

%% s.43, suoran ja ympyran leikkauspisteet
clear
close all
Ax=-1
Ay=0
theta=60
Px=3
Py=1 
r=2

[Qx,Qy,u,v]=suorien_leikkauspiste(Ax,Ay,theta,Px,Py,theta+90)

PQ=sqrt((Qx-Px)^2+(Qy-Py)^2)

th=0:1:360;
cx=Px+r*cosd(th);
cy=Py+r*sind(th);

plot(Ax,Ay,'r.','markersize',20)
hold on
plot([Ax,Ax+cosd(theta)],[Ay,Ay+sind(theta)],'r','linewidth',3)
plot(cx,cy,'k','linewidth',1.5)
plot([Px,Qx],[Py,Qy],'c','linewidth',1.5)
plot(Px,Py,'k.','markersize',20)
plot(Qx,Qy,'c.','markersize',20)

if PQ<=r %suora ja ympyra leikkaavat

QS=sqrt(r^2-PQ^2)
QT=QS

Sx=Qx-QS*cosd(theta)
Sy=Qy-QS*sind(theta)

Tx=Qx+QT*cosd(theta)
Ty=Qy+QT*sind(theta)

plot([Ax,Sx,Tx],[Ay,Sy,Ty],'r','linewidth',1.5)  
plot([Sx,Px,Tx],[Sy,Py,Ty],'g','linewidth',1.5)
plot([Sx,Tx],[Sy,Ty],'g.','markersize',20)
else
plot([Ax,Qx],[Ay,Qy],'r','linewidth',1.5)
end 
hold off
grid
axis equal




%% s.47, kasivarsi
%suora kinematiikka
clear
%varsien pituudet
OM=5
MP=3
alfa=130
beta=40

Mx=OM*cosd(alfa)
My=OM*sind(alfa)
Px=Mx+MP*cosd(alfa+beta)
Py=My+MP*sind(alfa+beta)
L=OM+MP


figure(1)
plot([0,Mx,Px],[0,My,Py],'r.-','linewidth',2,'markersize',20)
hold 
plot([-L,L],[0,0],'k')
plot([0,0],[-L,L],'k')
hold off
grid
axis([-L,L,-L,L])
axis square
title(['OM = ',num2str(OM),', MP = ',num2str(MP),...
       ', \alpha = ',num2str(alfa),'^o, \beta = ',num2str(beta),...
       '^o: Px = ',num2str(Px),', Py = ',num2str(Py)])

%% kaanteinen kinematiikka 
clear
%varsien pituudet
OM=5
MP=3
Px=-6
Py=-2

th=atan2d(Py,Px)
OP=sqrt(Px^2+Py^2)
gamma=acosd((OP^2+OM^2-MP^2)/(2*OP*OM))
alfa=th-gamma
beta=180-acosd((MP^2+OM^2-OP^2)/(2*MP*OM))

Mx=OM*cosd(alfa)
My=OM*sind(alfa)
Px=Mx+MP*cosd(alfa+beta)
Py=My+MP*sind(alfa+beta)
L=OM+MP
r=OM-MP
R=OM+MP
k=0:1:360;
x=r*cosd(k);
y=r*sind(k);
X=R*cosd(k);
Y=R*sind(k);

figure(2)
plot([0,Mx,Px],[0,My,Py],'r.-','linewidth',2,'markersize',20)
hold 
plot([-L,L],[0,0],'k')
plot([0,0],[-L,L],'k')
plot(x,y,'b',X,Y,'b')
hold off
grid

axis([-L,L,-L,L])
axis square
title(['OM = ',num2str(OM),', MP = ',num2str(MP),...
       ', Px = ',num2str(Px),', Py = ',num2str(Py),...
       ': \alpha = ',num2str(alfa),'^o, \beta = ',num2str(beta),...
       '^o'])

%% s.51
clear
%close all
%mitat
OA=50
OB=10
AC=70
CD=40
h=30

%kulma
alfa=320

%koordinaatit
Ox=0
Oy=0
Ax=0
Ay=-OA
Bx=OB*cosd(-90+alfa)
By=OB*sind(-90+alfa)
%AB:n ja AC:n suuntakulma
th=atan2d(By-Ay,Bx-Ax)
Cx=Ax+AC*cosd(th)
Cy=Ay+AC*sind(th)
Dy=h
Dx=Cx-sqrt(CD^2-(Dy-Cy)^2)

plot([Ox,Bx],[Oy,By],'k.-','linewidth',2,'markersize',20)
hold on
plot([Ax,Cx],[Ay,Cy],'b.-','linewidth',2,'markersize',20)
plot([Cx,Dx],[Cy,Dy],'r.-','linewidth',2,'markersize',20)
hold off
grid
axis equal

%% animaatio
clear
%mitat
OA=50
OB=10
AC=70
CD=40
h=30
%koordinaatit
Ox=0
Oy=0
Ax=0
Ay=-OA
%kulma
for alfa=0:1:(2*360) %alfa=0,1,2,...,360
Bx=OB*cosd(-90+alfa)
By=OB*sind(-90+alfa)
%AB:n ja AC:n suuntakulma
th=atan2d(By-Ay,Bx-Ax)
Cx=Ax+AC*cosd(th)
Cy=Ay+AC*sind(th)
Dy=h
Dx=Cx-sqrt(CD^2-(Dy-Cy)^2)
plot([Ox,Bx],[Oy,By],'k.-','linewidth',2,'markersize',20)
hold on
plot([Ax,Cx],[Ay,Cy],'b.-','linewidth',2,'markersize',20)
plot([Cx,Dx],[Cy,Dy],'r.-','linewidth',2,'markersize',20)
plot([-70,30],[h,h])
hold off
grid
%axis equal
axis([-70,30,-60,40])
axis square
pause(0.001)
end



%% s.55, four-bar
clear
OA=1.5
AB=3.5
BC=4
OC=5

alfa=0:1:360;
Ox=0
Oy=0
Cx=OC
Cy=0
Ax=OA*cosd(alfa);
Ay=OA*sind(alfa);
AC=sqrt((Ax-Cx).^2+(Ay-Cy).^2);
theta=atan2d(Ay-Cy,Ax-Cx);
phi=acosd((AC.^2+BC^2-AB^2)./(2*AC*BC));
beta=theta-phi;
Bx=Cx+BC*cosd(beta);
By=Cy+BC*sind(beta);
Px=(Ax+Bx)/2;
Py=(Ay+By)/2;

Px=1/3*Ax+2/3*Bx;
Py=1/3*Ay+2/3*By;
%% kuva asennossa numero k
k=72
plot([Ox,Ax(k),Bx(k),Cx],[Oy,Ay(k),By(k),Cy],'k.-','linewidth',1.5,'markersize',20)
hold on
plot(Px(k),Py(k),'r.','markersize',20)
plot(Px,Py,'r','linewidth',1.3)
R=2
plot([-R,OC+R],[0,0],'k')
hold off
grid
axis equal
title(['\alpha = ',num2str(alfa(k)),'^o, \beta = ',num2str(beta(k)),'^o'])

%% animaatio

for k=1:1:length(alfa)
plot([Ox,Ax(k),Bx(k),Cx],[Oy,Ay(k),By(k),Cy],'k.-','linewidth',1.5,'markersize',20)
hold on
plot(Px(k),Py(k),'r.','markersize',20)
plot(Px,Py,'r','linewidth',1.3)
R=2
plot([-R,OC+R],[0,0],'k')
hold off
grid
%axis equal
axis([-2,6,-3,5])
axis square
pause(0.01)
end

%% s.61 , Peaucellier
clear
%varsien pituudet
L=5
r=1.5
a=2

%Q liikkuu pitkin suoraa x=xQ
xQ=(L^2-r^2)/(2*a)

%CP:n maksimikulma
alfamax=180-acosd((2*a^2-(L-r)^2)/(2*a^2))-0.001;
%Qy:n maksimiarvo
hmax=(L+r)*sind(alfamax/2);

alfa=-alfamax:alfamax/200:alfamax;

%koordinaatit 
Cx=a
Cy=0

Px=Cx+a*cosd(alfa);
Py=Cy+a*sind(alfa);
OP=sqrt(Px.^2+Py.^2);
thOP=atan2d(Py,Px);
beta=acosd((L^2+OP.^2-r^2)./(2*L*OP));
thOA=thOP+beta;
Ax=L*cosd(thOA);
Ay=L*sind(thOA);
thOB=thOP-beta;
Bx=L*cosd(thOB);
By=L*sind(thOB);
AB=sqrt((Bx-Ax).^2+(By-Ay).^2);
thAB=atan2d(By-Ay,Bx-Ax);
gamma=acosd(AB/(2*r));
thAQ=thAB+gamma;
Qx=Ax+r*cosd(thAQ);
Qy=Ay+r*sind(thAQ);

for k=1:length(alfa)
%k=241
plot([Ax(k),Px(k),Bx(k),Qx(k),Ax(k)],[Ay(k),Py(k),By(k),Qy(k),Ay(k)],'g','linewidth',1.5)
hold
plot([Ax(k),0,Bx(k)],[Ay(k),0,By(k)],'b.-','linewidth',1.5,'markersize',20)
plot([Cx,Px(k)],[Cy,Py(k)],'r.-','linewidth',1.5,'markersize',20)
plot(Px(k),Py(k),'m.','markersize',25)
plot(Qx(k),Qy(k),'k.','markersize',25)
plot([xQ,xQ],[-hmax,hmax],'k','linewidth',1.5)
plot(Cx+a*cosd(alfa),Cy+a*sind(alfa),'r')
hold off
grid
axis equal
pause(0.01)
end