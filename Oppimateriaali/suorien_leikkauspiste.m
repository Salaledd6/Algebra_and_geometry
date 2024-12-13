function [Px,Py,r,t]=suorien_leikkauspiste(Ax,Ay,theta,Bx,By,delta)
%P=[Px,Py] on suorien A,theta ja B,delta leikkauspiste
%Px=A+r*cos(theta)=Bx+t*cos(delta),Py=Ay+r*sin(theta)=By+t*sin(delta)
a=cosd(theta);
b=-cosd(delta);
c=sind(theta);
d=-sind(delta);
e=Bx-Ax;
f=By-Ay;

r=(d*e-b*f)/(a*d-b*c);
t=(a*f-c*e)/(a*d-b*c);

Px=Ax+r*cosd(theta);
Py=Ay+r*sind(theta);
