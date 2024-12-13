function [Px,Py,r]=ympyraABC(Ax,Ay,Bx,By,Cx,Cy)
%P=[Px,Py] on pisteiden A,B ja C kautta kulkevan ympyran kp ja r sen sade

%sivun AB keskipiste 
Fx=(Ax+Bx)/2;
Fy=(Ay+By)/2;
%ja suuntakulma
kAB=atan2d(By-Ay,Bx-Ax);
theta=kAB+90; %keskinormaalin suuntakulma

%sivun AC keskipiste 
Ex=(Ax+Cx)/2;
Ey=(Ay+Cy)/2;
%ja suuntakulma
kAC=atan2d(Cy-Ay,Cx-Ax);
delta=kAC+90; %keskinormaalin suuntakulma

%ympyran keskipiste eli keskinormaalien F,theta ja E,delta leikkauspiste P
[Px,Py,u,v]=suorien_leikkauspiste(Fx,Fy,theta,Ex,Ey,delta);

r=sqrt((Px-Ax)^2+(Py-Ay)^2);