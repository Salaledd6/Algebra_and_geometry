%% s.9
clear %tyhjentaa muuttujat
syms x %luodaan kirjainmuuttuja x
simplify((2*x+3)*(4-5*x))
expand((2*x+3)*(4-5*x))
%-10*x^2-7*x+12
%%
x=0.73
(2*x+3)*(4-5*x)
-10*x^2-7*x+12
%% s.12
clear 
syms x
simplify(2*x*y^2+4*x^2*y)
factor(2*x*y^2+4*x^2*y)
%% s.16
clear
syms a
simplify(9*a^4-4)
factor(9*a^4-4)
%% s.19
clear
syms a b
simplify(2*a*b/(a^2+a*b))
%% s.23
clear
syms R H 
simplify(R-R^2/(R+H))

%%  s.27
clear
syms s
simplify((1+s/3)/(1+2*s/(s+3)))

%% s.29-30
clear
syms L v1 v2
simplify(2*L/(L/v1+L/v2))
simplify(2*v1*v2/(v1+v2)-(v1+v2)/2)

%%
v1=80
v2=60
ka=(v1+v2)/2
kn=2*v1*v2/(v1+v2)

%% s.34
clear
syms R s t positive
simplify(sqrt(R^2-(R/2)^2))
simplify(s*t*sqrt(1/t^2-1/s^2))
%% s.39
clear
syms x
solve((2*x-1)/(x+4)==3,x)
%% tarkastus
x=-13
(2*x-1)/(x+4)
%% s.41
syms E e R r
solve(E/e==(R+r)/r,r)
%r=(R*e)/(E - e)
%% tarkastus
E=6.3
e=1.5
R=3.7
r=(R*e)/(E - e)

E/e
(R+r)/r
%% s.43
clear
syms v a1 a2 t1 t2
solve(v==a1*t1+a2*(t2-t1),t1)
%t1=(v - a2*t2)/(a1 - a2)
%% tarkastus

v=8.3
a1=0.4
a2=6.2
t2=9.8
t1=(a2*t2-v)/(a2-a1)

a1*t1+a2*(t2-t1)
v

%% s.46
clear
syms tA tB tAB
solve(1/tA+1/tB==1/tAB,tB)
simplify(ans)
%tB=(tA*tAB)/(tA - tAB)
%% tarkastus
tA=80
tAB=48
tB=tA*tAB/(tA-tAB)

1/tA+1/tB
1/tAB

%% s.48
clear
syms p1 p2 M x y
solve(p1*M/x==p2,x)
solve((p1*M+y)/(M+y)==p2,y)
%x=(M*p1)/p2
%y=(M*p1 - M*p2)/(p2 - 1)
%% tarkastus

p1=0.045
M=6.00
p2=0.11
x=(M*p1)/p2
y=(M*p1 - M*p2)/(p2 - 1)

p1*M/x
p2

(p1*M+y)/(M+y)
p2

%% s.50
clear
syms r1 r2 r12 x
solve(x/r1+(1-x)/r2==1/r12,x)
%x=-(r1*r2 - r1*r12)/(r1*r12 - r2*r12)
%% tarkastus

r1=0.79
r2=1.00
r12=0.875
x=-(r1*r2 - r1*r12)/(r1*r12 - r2*r12)

x/r1+(1-x)/r2
1/r12
%% s.51
clear
a=3
b=2
c=-1
D=b^2-4*a*c
if D>=0
    x1=(-b-sqrt(D))/(2*a)
    x2=(-b+sqrt(D))/(2*a)
end
%% s.55
syms x
solve(3*x^2+2*x-1==0,x)

%% s.58
syms h v t g
solve(h+v*t-1/2*g*t^2==0,t)
%%
h=10
v=2
g=9.81
t=(v+sqrt(v^2+2*g*h))/g
%% s.60
clear
d=0.001
r=0.3
n=100
s=pi*d*n^2+pi*(2*r-d)*n

%% ratkaistaan n yhtalosta a*n^2+b*n+c=0
a=pi*d
b=pi*(2*r-d)
c=-s
n=(-b+sqrt(b^2-4*a*c))/(2*a)
%% 
clear
syms p d r s n
solve(s==p*d*n^2+p*(2*r-d)*n,n)
%% s.61
clear
syms v
solve(sqrt(3*v+4)==1+2*v,v)
%% s.62
clear
syms v p p0 r
solve(v==sqrt(2*(p-p0)/r),p0)
%% s.64
clear
syms x y
solve(2*x+5*y==1,x+3*y==2,x,y)
x=ans.x
y=ans.y

%% s.67
clear
syms x y
yht1=0.20*x+0.70*y+0.20*(2-x-y)==0.40*2;
yht2=0.50*x+0.10*y+0.20*(2-x-y)==0.2725*2;
solve(yht1,yht2,x,y)
x=ans.x
y=ans.y
x=double(x)
y=double(y)
2-x-y
%% s.68
clear
syms PC PZ P1C P1Z P2C P2Z P3C P3Z M x y
solve(P1C*x+P2C*y+P3C*(M-x-y)==PC*M,P1Z*x+P2Z*y+P3Z*(M-x-y)==PZ*M,x,y)
x=ans.x
y=ans.y
%x = (M*P2C*P3Z - M*P3C*P2Z + M*P2Z*PC - M*P3Z*PC - M*P2C*PZ + M*P3C*PZ)/(P1C*P2Z - P2C*P1Z - P1C*P3Z + P3C*P1Z + P2C*P3Z - P3C*P2Z)
%y =-(M*P1C*P3Z - M*P3C*P1Z + M*P1Z*PC - M*P3Z*PC - M*P1C*PZ + M*P3C*PZ)/(P1C*P2Z - P2C*P1Z - P1C*P3Z + P3C*P1Z + P2C*P3Z - P3C*P2Z)
 

%%

P1C=0.20
P1Z=0.50
P2C=0.70
P2Z=0.10
P3C=0.20
P3Z=0.20
PC=0.40
PZ=0.2725
M=2

x =(M*P2C*P3Z - M*P3C*P2Z + M*P2Z*PC - M*P3Z*PC - M*P2C*PZ + M*P3C*PZ)/(P1C*P2Z - P2C*P1Z - P1C*P3Z + P3C*P1Z + P2C*P3Z - P3C*P2Z)  
y = -(M*P1C*P3Z - M*P3C*P1Z + M*P1Z*PC - M*P3Z*PC - M*P1C*PZ + M*P3C*PZ)/(P1C*P2Z - P2C*P1Z - P1C*P3Z + P3C*P1Z + P2C*P3Z - P3C*P2Z) 
M-x-y
%% s.69 
clear
syms x y
yht1=1.00*x+1.26*y==1.11*(x+y)
yht2=1.00*x+1.26*y+0.79*(200-x-y)==201
solve(yht1,yht2,x,y)
x=ans.x
y=ans.y
x=double(x)
y=double(y)
200-x-y
%% s.70
clear
syms rv rg rvg re V M x y
solve(rv*x+rg*y==rvg*(x+y),rv*x+rg*y+re*(V-x-y)==M,x,y)
x=ans.x
y=ans.y
%x = -((rg - rvg)*(M - V*re))/(re*rg - re*rv - rg*rvg + rv*rvg)
%y = ((rv - rvg)*(M - V*re))/(re*rg - re*rv - rg*rvg + rv*rvg)
 
%%

rv=1.00
rg=1.26
re=0.79
rvg=1.11
V=200
M=201
x=-((rg - rvg)*(M - V*re))/(re*rg - re*rv - rg*rvg + rv*rvg)
y=((rv - rvg)*(M - V*re))/(re*rg - re*rv - rg*rvg + rv*rvg)
V-x-y

%% s.72
clear
R1=1
R2=2
R3=3
R4=4
R5=5
R6=6
V=7

syms I1 I2 I3
yht1=V-R1*I1-R2*(I1-I2)-R3*I1==0
yht2=-R4*(I2-I3)-R2*(I2-I1)==0
yht3=-R5*I3-R6*I3-R4*(I3-I2)==0
solve(yht1,yht2,yht3,I1,I2,I3)
I1=ans.I1
I2=ans.I2
I3=ans.I3
%% desimaaliluvuiksi
double(I1)
double(I2)
double(I3)
%%
syms I1 I2 I3 R1 R2 R3 R4 R5 R6 V
yht1=V-R1*I1-R2*(I1-I2)-R3*I1==0
yht2=-R4*(I2-I3)-R2*(I2-I1)==0
yht3=-R5*I3-R6*I3-R4*(I3-I2)==0
solve(yht1,yht2,yht3,I1,I2,I3)
I1=ans.I1
I2=ans.I2
I3=ans.I3


%%
R1=1
R2=2
R3=3
R4=4
R5=5
R6=6
V=7

I1=(V*(R2*R4 + R2*R5 + R2*R6 + R4*R5 + R4*R6))/(R1*R2*R4 + R1*R2*R5 + R1*R2*R6 + R2*R3*R4 + R1*R4*R5 + R2*R3*R5 + R1*R4*R6 + R2*R3*R6 + R2*R4*R5 + R2*R4*R6 + R3*R4*R5 + R3*R4*R6)
 
I2=(R2*V*(R4 + R5 + R6))/(R1*R2*R4 + R1*R2*R5 + R1*R2*R6 + R2*R3*R4 + R1*R4*R5 + R2*R3*R5 + R1*R4*R6 + R2*R3*R6 + R2*R4*R5 + R2*R4*R6 + R3*R4*R5 + R3*R4*R6)
 
I3=(R2*R4*V)/(R1*R2*R4 + R1*R2*R5 + R1*R2*R6 + R2*R3*R4 + R1*R4*R5 + R2*R3*R5 + R1*R4*R6 + R2*R3*R6 + R2*R4*R5 + R2*R4*R6 + R3*R4*R5 + R3*R4*R6)
 
