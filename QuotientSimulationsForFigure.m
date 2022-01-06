%% VN quotient simulations
close all
clc
clear

perturbation = 0.01;
endTime = 10000;
MarkerSize = 500;

w1_Delta = 0.02;
w2_Delta = 1;

w1_Ecad = 0.25;
w2_Ecad = 1;

n1 = 2;
n2 = 1;

Ecad_W = Quotient_bilayer_adj(n1,n2,n1,n2,w1_Ecad,w2_Ecad);
Delta_W = Quotient_bilayer_adj(n1,n2,n1,n2,w1_Delta,w2_Delta);

%import ODE parameters
NEDparameters;
paras = [a1,b1,a2,b2,b3,k1,k2,h1,h2,h3];

%FIND THE HSS BEFORE THIS
N_ss = fsolve(@(N) HSS_notch(N,a1,a2,b1,b2,k1,k2,h1,h2,b3,h3), 0.5);
E_ss = (N_ss.^k2)./(a2 + N_ss.^k2);
D_ss = 1./(1 + b3.*(N_ss.^h3));

%initial conditions
lum_IC = [ min(N_ss+ N_ss*perturbation,1),min(E_ss+E_ss*perturbation ,1) ,max(D_ss-D_ss*perturbation,0)] ;
bas_IC = [max(N_ss - N_ss*perturbation,0) , max(E_ss - E_ss*perturbation,0) , min(D_ss + D_ss*perturbation,1)  ];

IC = [lum_IC,  bas_IC  ];

time = [0,endTime];
num_state_variables = 3*length(Ecad_W(:,1));

%solve ODE for stiff solver
[t,y] = ode15s(@(t,y) Ecad_quotient_bilayer(t,y,paras,Ecad_W ,Delta_W,num_state_variables), time,IC);

%plot quotient bilayer
timeStep = 1;
median(divisors(length(t)));
dt = 1:timeStep:length(t);

fig=gcf;
fig.Position(3:4)=[0.2,0.1];
colormap parula
colorsLines = [0,0,0; 0.7529,0,0];
plot([0.997,0.997],[1,1.03], 'linewidth',3,'color', [0,0,0]);
hold on
plot([1.003,1.003],[1,1.03], 'linewidth',3,'color', [0.7529,0,0]);
scatter([1,1],[1,1.03],MarkerSize,y(dt(end),1:3:end-2),'filled','linewidth',1,'MarkerEdgeColor','k')
c = colorbar('location','south');
caxis([0,0.5])
set(c,'TickLabelInterpreter','latex')
set(c,'visible','off')
axis off
xlim([0.9,1.1])
ylim([0.99,1.05])


%% Tri quotient simulations
close all
clc
clear

perturbation = 0.1;
endTime = 10000;
MarkerSize = 500;

w1_Delta = 0.05;
w2_Delta = 1;

w1_Ecad = 0.5;
w2_Ecad = 1;

n1 = 2;
n2 = 2;

Ecad_W = Quotient_bilayer_adj(n1,n2,n1,n2,w1_Ecad,w2_Ecad);
Delta_W = Quotient_bilayer_adj(n1,n2,n1,n2,w1_Delta,w2_Delta);

%import ODE parameters
NEDparameters;
paras = [a1,b1,a2,b2,b3,k1,k2,h1,h2,h3];

%FIND THE HSS BEFORE THIS
N_ss = fsolve(@(N) HSS_notch(N,a1,a2,b1,b2,k1,k2,h1,h2,b3,h3), 0.5);
E_ss = (N_ss.^k2)./(a2 + N_ss.^k2);
D_ss = 1./(1 + b3.*(N_ss.^h3));

%initial conditions
lum_IC = [ min(N_ss+ N_ss*perturbation,1),min(E_ss+E_ss*perturbation ,1) ,max(D_ss-D_ss*perturbation,0)] ;
bas_IC = [max(N_ss - N_ss*perturbation,0) , max(E_ss - E_ss*perturbation,0) , min(D_ss + D_ss*perturbation,1)  ];

IC = [lum_IC,  bas_IC  ];

time = [0,endTime];
num_state_variables = 3*length(Ecad_W(:,1));

%solve ODE for stiff solver
[t,y] = ode15s(@(t,y) Ecad_quotient_bilayer(t,y,paras,Ecad_W ,Delta_W,num_state_variables), time,IC);

%plot quotient bilayer
timeStep = 1;
median(divisors(length(t)));
dt = 1:timeStep:length(t);

fig=gcf;
fig.Position(3:4)=[0.2,0.1];
colormap parula
colorsLines = [0,0,0; 0.7529,0,0];
plot([0.997,0.997],[1,1.03], 'linewidth',3,'color', [0,0,0]);
hold on
plot([1.003,1.003],[1,1.03], 'linewidth',3,'color', [0.7529,0,0]);
scatter([1,1],[1,1.03],MarkerSize,y(dt(end),1:3:end-2),'filled','linewidth',1,'MarkerEdgeColor','k')
c = colorbar('location','south');
caxis([0,0.5])
set(c,'TickLabelInterpreter','latex')
set(c,'visible','off')
axis off
xlim([0.9,1.1])
ylim([0.99,1.05])


%% functions

function out = HSS_notch(N,a1,a2,b1,b2,k1,k2,h1,h2,b3,h3)

U2 = 1./(1 + b3.*(N.^h3));
U1 = (N.^k2)./(a2 + N.^k2);
E = U1;

out = ((U2.^k1)./(a1 + U2.^k1)).*( 1./(1 + b1.*((U1).^h1))).*( 1./(1 + b2.*((E).^h2))) - N;

end