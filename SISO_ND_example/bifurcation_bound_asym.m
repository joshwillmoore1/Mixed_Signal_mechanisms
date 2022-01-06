close all
clc
clear

nb1 = 2;
nl1 = 2;
T = 1.5316;
w2 = 1;
nl2 = linspace(0,2,100);
nb2 = linspace(0,4,100);

dsc_nl2 = [1,2,3];

fs = 28;
[N1,N2] = meshgrid(nl2,nb2);

[C,h] = contourf(N1,N2,asymm_upper_bound(2,N2,2,N1,T,w2),...
    'linewidth',1.5,'color','w');
hcb = colorbar('Location','Northoutside');
set(hcb,'TickLabelInterpreter','latex','ticks',[0,0.25])

xlabel("$n_{L,2}$", 'FontSize',fs)
ylabel("$n_{B,2}$",'FontSize',fs)
hold on
scatter([1,1,1],[1,2,3],300,'filled','x','MarkerEdgeColor','w','linewidth',2)
box off
yticks([0, 2, 4])
xticks([0, 1, 2])

ax = gca; 
ax.FontSize = fs;

function out = asymm_upper_bound(nb1,nb2,nl1,nl2,T,w2)

out = ( (  -(nb1.*nl2 + nl1.*nb2) + sqrt( (nb1.*nl2 + nl1.*nb2).^2 +...
    4*(T^2 - 1).*nb1.*nb2.*nl1.*nl2  ) )./(2.*(T+1).*nb1.*nl1)).*w2;

end