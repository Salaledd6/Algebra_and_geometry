%% s.7-8
clear %tyhjentaa muuttujat
h=3.5
x=4.5
dx=1.3
dy=sqrt((x+dx)^2+h^2)-sqrt(x^2+h^2)
%sqrt(...) = neliojuuri(...)
%% s.9-10
%MATLAB + Symbolic math toolbox
clear
syms R h L %kirjainmuuttujat
solve((R+h)^2==R^2+L^2,L)
%h^(1/2)*(2*R + h)^(1/2)
%-h^(1/2)*(2*R + h)^(1/2)

%%
clear
R=6400
h=0.1
L=sqrt(2*R*h+h^2)
sqrt(h)*sqrt(h+2*R)
h^(1/2)*(2*R + h)^(1/2)
%% s.13-14
clear
a=12
b=7
c=8
h=2
L=sqrt(h^2+((a-c)/2)^2)

u=a-c
L=sqrt(h^2+(u/2)^2)

A1=1/2*b*L
D=sqrt(h^2+(b/2)^2)
A2=1/2*(a+c)*D
A=2*A1+2*A2

%% s.15-16
%MATLAB + Symbolic math toolbox
clear
syms H v1 v2 L
solve(H/v1==sqrt(L^2+H^2)/v2,H)
simplify(ans)
%%
clear
v1=2
v2=7
L=5
H=v1*L/sqrt(v2^2-v1^2)

%tarkastus
H/v1
sqrt(L^2+H^2)/v2

%% s.17-18
clear
a=5
b=6
x=2

MA=sqrt(x^2+(b/2)^2)
MB=MA
DC=sqrt(a^2-(b/2)^2)
MC=DC-x

s=MA+MB+MC
%%  s:n kuvaaja, kun x = 0 ... DC
clear
a=5
b=6
DC=sqrt(a^2-(b/2)^2)
%lista (vektori) x:n arvoja tasavalein 0...DC
x=0:DC/1000:DC; %alku:askel:loppu
%%
MA=sqrt(x.^2+(b./2).^2); %huom: .^2
MB=MA;
MC=DC-x;
s=MA+MB+MC
%%
%s:n pienin arvo
[smin,indmin]=min(s) %max(s) = suurin arvo
%smin = pienin arvo
%indmin = sen jarjestysnumero
xmin=x(indmin) %vastaava x:n arvo
%x=[x(1),x(2),...,x(n)]
%%
%kuvaaja
plot(x,s,'r','linewidth',2) %plot(vaaka,pysty), r = red = punainen
hold on %pito paalle (ensimmaisen plotin jalkeen)
plot(xmin,smin,'k.','markersize',20) %k. = musta pallero
hold off %pito pois (viimeisen plotin jalkeen)
xlabel('korkeus x') %vaaka-akselin otsikko
ylabel('s = MA+MB+MC') %pysty-akselin otsikko
title('kuvan otsikko')
grid %taustaristikko

%% s.19-20
clear
R=1
r=0.3*R
H=2*sqrt(R^2-r^2)
Vlierio=pi*r^2*H
h=R-H/2
Vsegm=pi*h^2*(R-h/3)
Asegm=2*pi*R*h
Vpois=Vlierio+2*Vsegm
Apois=2*Asegm
Vpallo=4/3*pi*R^3
Apallo=4*pi*R^2
Vsuhde=Vpois/Vpallo
Asuhde=Apois/Apallo
%% Vsuhteen ja Asuhteen kuvaajat, kun r = 0...R
clear
R=1
r=0:R/1000:R;
H=2*sqrt(R^2-r.^2); %huom: .^2
Vlierio=pi*r.^2.*H; %huom: .^2 ja .*
h=R-H/2;
Vsegm=pi*h.^2.*(R-h/3); %huom: .^ja .*
Asegm=2*pi*R*h;
Vpois=Vlierio+2*Vsegm;
Apois=2*Asegm;
Vpallo=4/3*pi*R^3;
Apallo=4*pi*R.^2;
Vsuhde=Vpois/Vpallo;
Asuhde=Apois/Apallo;

plot(r/R,Vsuhde,'r','linewidth',1.5)
hold on
plot(r/R,Asuhde,'b','linewidth',1.5)
hold off
grid
%xticks(0:0.1:1) %yticks
legend('Vsuhde','Asuhde','location','northwest')
xlabel('r/R')
ylim([0,1]) %xlim

%% s.21-22
clear
alfa=0:1:360;
r_kerroin=alfa./360;
h_kerroin=sqrt(1-(alfa./360).^2);

plot(alfa,r_kerroin,'b','linewidth',1.5)
hold on
plot(alfa,h_kerroin,'r','linewidth',1.5)
hold off
grid
xlim([0,360])
xticks(0:30:360)

legend('r/R','h/R')
xlabel('kulma \alpha')

%% s.27-28
clear
alfa=20
sind(alfa) %kulma asteina
sin(alfa*pi/180) %kulma radiaaneina
%% s.37-38
clear
OM=4
MP=3
alfa=30
beta=20
OQ=OM*cosd(alfa)
MN=MP*cosd(alfa+beta)
x=OQ+MN
QM=OM*sind(alfa)
NP=MP*sind(alfa+beta)
y=QM+NP
%% s.39-40
clear
R=2
L=5
alfa=48

OQ=R*cosd(alfa)
QP=R*sind(alfa)
QM=sqrt(L^2-QP^2)
OM=OQ+QM

%% s.41-42
clear
a=1
b=2
alfa=1:0.01:89;
AP=b./sind(alfa);
PB=a./cosd(alfa);
AB=AP+PB;
%%
[ABmin,indmin]=min(AB)
alfamin=alfa(indmin)

plot(alfa,AB,'linewidth',1.5)
hold
plot(alfamin,ABmin,'r.','markersize',20)
hold off
grid
xlabel('kulma \alpha')
ylabel('pituus AB')
ylim([0,20])
%% s.43-44
clear
AB=200
AC=80
BD=40
theta=20
CE=AC*cosd(theta)
AE=AC*sind(theta)
FB=CE-BD
AF=sqrt(AB^2-FB^2)
CD=AF-AE
%% s.45-46
clear
n=3:1:20;
kerroin=(1-sind(180./n))./sind(180./n);

plot(n,kerroin,'r.-','linewidth',1.5,'markersize',20)
grid
xlabel('lukumaara n')
ylabel('suhde R/r')
xticks(n)


%% s.47-48
clear
syms H x L a b
solve(H/x==tan(b),H/(L+x)==tan(a),H,x)
H=ans.H
x=ans.x

%H=-(L*tan(a)*tan(b))/(tan(a) - tan(b))
%x=-(L*tan(a))/(tan(a) - tan(b))

%%
clear
L=2
a=30*pi/180
b=50*pi/180
H=-(L*tan(a)*tan(b))/(tan(a) - tan(b))
x=-(L*tan(a))/(tan(a) - tan(b))

%tarkastus
H/x
tan(b)
H/(L+x)
tan(a)

%% s.55-56
clear
r=2
R=3
H=4
h=0:H/1000:H;
x=h/H*(R-r);
s=r+x;
Vastia=1/3*pi*(R^2+R*r+r^2)*H;
Vneste=1/3*pi*(s.^2+s*r+r^2).*h;
plot(h,Vneste,'linewidth',1.5)
grid
xlabel('korkeus h')
ylabel('Vneste')
title(['r = ',num2str(r),', R = ',num2str(R),', H = ',num2str(H),', Vastia = ',num2str(Vastia)])
%% s.57-58
clear
r=1
h=3
H=5
L=2*pi*r
x=0:L/100:L;%kaaren pituus
alfa=x/r;  %kulma (rad)
w=r-r*cos(alfa);
z=w/(2*r)*(H-h);
y=h+z;

plot(x,y,'linewidth',1.5)
grid
ylim([0,H])
axis equal %akseleiden mittakaavat samoiksi

%% s.59-60
asind(0.35) %kulma asteina
asin(0.35) %kulma radiaaneina

%% s.65-66
clear
H=5
h=2
x=0.01:0.01:10;
alfa=atand(H./x)-atand(h./x);
[alfamax,indmax]=max(alfa)
xmax=x(indmax)

plot(x,alfa,'linewidth',1.5)
hold
plot(xmax,alfamax,'r.','markersize',20)
hold off
grid


%% s.69-70
clear
r=1
h=0:r/100:2*r;
alfa=2*acos((r-h)/r); %rad
Vsuhde=1/(2*pi)*(alfa-sin(alfa));

plot(h/r,Vsuhde,'linewidth',1.5)
grid
xlabel('h/r')
ylabel('Vneste/Vsailio')
%% s.71-72
clear
R=3
r=1
L=5
P1Q1=sqrt(L^2-(R-r)^2)
P2Q2=P1Q1
alfa=acos((R-r)/L)
Q1QQ2=2*alfa*r
P1PP2=(2*pi-2*alfa)*R
hihna=P1Q1+P2Q2+Q1QQ2+P1PP2

%% s.75-76
clear
syms x a b h L
solve(tan(b)==x/h,tan(a+b)==(x+L)/h,x,b)
x=ans.x
b=ans.b

%x =
% -(L*tan(a) - (tan(a)*(tan(a)*L^2 + 4*L*h - 4*tan(a)*h^2))^(1/2))/(2*tan(a))
% -(L*tan(a) + (tan(a)*(tan(a)*L^2 + 4*L*h - 4*tan(a)*h^2))^(1/2))/(2*tan(a))
 
 
%b =
% -atan((L*tan(a) - (tan(a)*(tan(a)*L^2 + 4*L*h - 4*tan(a)*h^2))^(1/2))/(2*h*tan(a)))
% -atan((L*tan(a) + (tan(a)*(tan(a)*L^2 + 4*L*h - 4*tan(a)*h^2))^(1/2))/(2*h*tan(a)))
 


%% 

h=4
a=30*pi/180
L=5

%x =
x1=-(L*tan(a) - (tan(a)*(tan(a)*L^2 + 4*L*h - 4*tan(a)*h^2))^(1/2))/(2*tan(a))
b1=atan(x1/h)
 
x2=-(L*tan(a) + (tan(a)*(tan(a)*L^2 + 4*L*h - 4*tan(a)*h^2))^(1/2))/(2*tan(a))
b2=atan(x2/h) 
 
%tarkastus
tan(a+b1)
(x1+L)/h
 
%% wolfram alphan ratkaisu
x1=1/2*(h*csc(a)*sqrt(sin(a)^2*(4*h*L*cot(a)-4*h^2+L^2)/h^2)-L)
b1=atan(x1/h)
 
x2=1/2*(-h*csc(a)*sqrt(sin(a)^2*(4*h*L*cot(a)-4*h^2+L^2)/h^2)-L)
b2=atan(x/h)
  

%% s.94-95
clear
a=4.3
b=6.2
gamma=25
c=sqrt(a^2+b^2-2*a*b*cosd(gamma))
%% s.95-96
clear
a=2.7
b=3.3
c=1.4
alfa=acosd((b^2+c^2-a^2)/(2*b*c))
beta=acosd((a^2+c^2-b^2)/(2*a*c))
gamma=acosd((a^2+b^2-c^2)/(2*a*b))
%% s.97-98
clear
OM=5
MP=3
x=4
y=4

OP=sqrt(x^2+y^2)
delta=atand(y/x)

gamma=acosd((OP^2+OM^2-MP^2)/(2*OP*OM))
alfa=delta-gamma
beta=180-acosd((OM^2+MP^2-OP^2)/(2*OM*MP))
%% s.99-100
clear
AB=200
AC=80
BD=40
CD=175
BC=sqrt(BD^2+CD^2)
alfa=atand(BD/CD)
beta=acosd((BC^2+AC^2-AB^2)/(2*BC*AC))
theta=alfa+beta-90
%% s.103-104
clear
a=8.6
alfa=54
beta=97
gamma=29

%b/sin(beta)=a/sin(alfa)
b=sind(beta)*a/sind(alfa)

%c/sin(gamma)=a/sin(alfa)
c=sind(gamma)*a/sind(alfa)

%% s.105-106
clear
L=2
alfa=30
beta=50
gamma=beta-alfa
AP=sind(alfa)*L/sind(gamma)
H=AP*sind(beta)
x=AP*cosd(beta)
%% s.113-114
clear
OA=320
OB=80
BC=120
dbeta=15
%alku
OC=sqrt(OB^2+BC^2)
alfa=asind(BC/OA)
beta=atand(BC/OB)
delta=beta
%loppu
gamma=asind(sind(beta+dbeta)/OA*OC)
theta=180-(beta+dbeta)-gamma
dalfa=(180-theta-delta)-alfa
















