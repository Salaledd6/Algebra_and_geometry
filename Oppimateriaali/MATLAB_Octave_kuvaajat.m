clear
x=[1,2,4,7,9] 
y=[-1,0,4,1,-2]
figure(1)
%plot(x,y,'linewidth',1.5)
%plot(x,y,'r','linewidth',1.5)
%plot(x,y,'r-o','linewidth',1.5,'markersize',10)
plot(x,y,'ks','markersize',10)
grid
xlabel('x-akselin otsikko')
ylabel('y-akselin otsikko')
title('kuvan otsikko')

%% 
clear
x=0:0.01:1; %puolipiste = ei tulostu komentoikkunaan

y1=x./sqrt(1+x.^2); %huomaa pisteet ./ ja .^
y2=x.*(1-x).^2; %huomaa pisteet .* ja .^

figure(2)
plot(x,y1,'r','linewidth',2) %r = punainen viiva
hold on %ensimmaisen plotin jalkeen
plot(x,y2,'k-.','linewidth',1)
hold off %viimeisen plotin jalkeen
grid %taustaristikko
legend('y_1','y_2') %"varikartta"
%legend({'y_1','y_2'},'location','northwest','fontsize',12) %"varikartta"
% xlim([0,1]) %x-akselin rajat, [vasen,oikea]
% ylim([0,0.8]) %y-akselin rajat, [ala,yla] 
% xticks(0:0.1:1) %x-akselin apuviivat
% yticks(0:0.1:0.8) %y-akselin apuviivat
xlabel('x')
ylabel('y','rotation',0)
title('y_1:n ja y_2:n kuvaajat')









