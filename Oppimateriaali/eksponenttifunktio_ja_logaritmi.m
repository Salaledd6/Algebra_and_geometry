%% s.3-4, y=a^x
clear
a=2
xmin=-3
xmax=3
x=xmin:(xmax-xmin)/100:xmax;
y=a.^x;

plot(x,y,'linewidth',1.5)
grid
xlabel('x')
ylabel('y')
title(['y = a^x, a = ',num2str(a)])

%% s.9-10, y=lg(x)
xmin=10^-3
xmax=10
x=xmin:(xmax-xmin)/1000:xmax;
y=log10(x);

plot(x,y,'linewidth',1.5)
grid



%% s.17-18, logaritminen asteikko 
clear
w0=10 %"vahvistuksen rajataajuus"
w=0:w0/10:10*w0;
%"vahvistus"
K=1./sqrt(1+(w/w0).^2);

figure(1)
plot(w,K)
grid
xlabel('w')
ylabel('K   ','rotation',0)

w=logspace(0,4,100); %w=10^0 ... 10^4, 100 arvoa tasavalein logaritmisella asteikolla
%eli lg(w):t tasavalein 0...4
K=1./sqrt(1+(w/w0).^2);
%vahvistus desibeleina
KdB=20*log10(K);


figure(2)
semilogx(w,KdB) %vaaka-akselilla logaritminen asteikko 
grid
xlabel('w')
ylabel('KdB ','rotation',0)

figure(3)
plot(log10(w),KdB)
grid
xlabel('lg(w)')
ylabel('KdB ','rotation',0)



%% s.19, Eulerin luku e=2.718....
e=exp(1) %exp(x)=e^x
n=1:500;
en=(1+1./n).^n;

plot(n,en)
grid
xlabel('n')
ylabel('(1+1/n)^n')
title(['e = ',num2str(e)])
%% s.23-24 jannite
clear
U0=5
tau=0.5
t=0:tau/100:5*tau;
U1=U0*exp(-t/tau);
U2=U0*(1-exp(-t/tau));

plot(t,U1,'b','linewidth',1.5)
hold 
plot(t,U2,'r','linewidth',1.5)
hold off 
grid
xlabel('aika t')
title(['U_0 = ',num2str(U0),', \tau = ',num2str(tau)])
legend({'purku','lataus'},'fontsize',12)

%% s.31-32, lampotila
clear
T0=50
Ty=20
k=1.5
tau=1/k
t=0:tau/100:5*tau;
T=Ty+(T0-Ty)*exp(-k*t);
plot(t,T,'linewidth',1.5)
grid
title(['T_0 = ',num2str(T0),', T_y = ',num2str(Ty),', k = ',num2str(k),', \tau = ',num2str(1/k)])
xlabel('aika t')

%% s.33-34 vaimennettu varahtely
A=3
tau=1.5
w=2*pi
phi=0


t=0:tau/100:5*tau;
y=A*exp(-t/tau).*sin(w*t+phi);
y1=A*exp(-t/tau);

plot(t,y,'b','linewidth',1.5)
hold on
plot(t,y1,'r',t,-y1,'r')
hold off
grid
xlabel('aika t')
%% s.39, piano keyboard frequencies

clear
r=2^(1/12)
f1=440/r^48
k=0:87;
f=f1*r.^k; %taajuudet


varit=ones(1,88); %1 = valkea
varit(2)=2; %2 = musta
for k=1:7
    varit(3+12*(k-1)+[2,4,7,9,11])=2;
end

varit(49)=3; %A = keltainen
varit(40)=4; %C = cyan

colors={'w','k','y','c'}


figure(1)
for k=1:88
bar(k,f(k),colors{varit(k)})
if k==1
    hold on
end
end
hold off
grid
k=-4:4;
yticks(2.^k*440)
set(gca,'yticklabel',{'','55','110','220','440','880','1760','3520'})


notes={'A','A#','B','C','C#','D','D#','E','F','F#','G','G#'}
for k=1:88
    n=mod(k,12); %jakojaannos n/12   
    if n==0
    n=12;
    end
    text(k-0.1,f(k)+20,notes{n},'rotation',90,'fontsize',5) 
end
xlim([-1,90])
xticks(1:12:88)
xlabel('k')
ylabel('f_k','rotation',0)
set(gca,'fontsize',6)
figure(2)
for k=1:88
bar(k,log(f(k))/log(2),colors{varit(k)})
if k==1
    hold on
end
end
hold off
grid

for k=1:88
    n=mod(k,12);    
    if n==0
    n=12;
    end
    text(k-0.1,log(f(k))/log(2)+0.05,notes{n},'rotation',90,'fontsize',5) 
end
xlabel('k')
ylabel('log_2(f_k)')
xlim([-1,90])
xticks(1:12:88)
yticks(0:13)
set(gca,'fontsize',6)