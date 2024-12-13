clear
v=[4,2,5,3,7]
N=length(v)
%v=[v(1),v(2),...,v(N)]
v(3)
%%
x=0:0.1:1 %alku:askel:loppu
% ; = rivi ei tulostu komentoikkunaan
%x=[0,0.1,0.2,0.3,...,0.9,1]

%%
%laskutoimitukset alkioittain
x+5 % = [x(1)+5,x(2)+5,...,x(N)+5]
x-5 % = [x(1)-5,x(2)-5,...,x(N)-5]
5*x % = [5*x(1),5*x(2),...,5*x(N)]
x/5 % = [x(1)/5,x(2)/5,...,x(N)/5]
sqrt(x) % = [sqrt(x(1)),sqrt(x(2)),...,sqrt(x(N))]
x+sqrt(x) % = [x(1)+sqrt(x(1)),x(2)+sqrt(x(2)),...,x(N)+sqrt(x(N))]
x-sqrt(x) % = [x(1)-sqrt(x(1)),x(2)-sqrt(x(2)),...,x(N)-sqrt(x(N))]


%HUOM !!!!!!!!!!!!!!!
% .^    .*   ./


x.^2  % = [x(1)^2,x(2)^2,...,x(N)^2]
%x^2 % -> virheilmoitus ( = x*x, * on matriisikertolasku)

x.*(1+x) % = [x(1)*(1+x(1)),x(2)*(1+x(2)),...,x(N)*(1+x(N))]
%x*(1+x) % -> virheilmoitus ( * on matriisikertolasku )

x./(1+x) % = [x(1)/(1+x(1)),x(2)/(1+x(2)),...,x(N)/(1+x(N))]
%x/(1+x) % -> tulos on 1 luku

5./(1+x) % = [5/(1+x(1)),5/(1+x(2)),...,5/(1+x(N))]
%5/(1+x) % -> virheilmoitus


