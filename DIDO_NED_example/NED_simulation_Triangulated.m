close all
clc
clear all

% roots of unity approach to equally spaced points on circle

%14
NumOfLuminal = 28;

rotationConstant = -0.12;
LumRadius = 1.8;
BasRadius = 2.5;

projectionDistance = 0.05;
perturbation = 0.2;
endTime = 1000;

w1_Delta = 0.1;
w2_Delta = 1;

w1_Ecad = 1;
w2_Ecad = 1;

AdjMatrixVN = bilayer_Adj_Matrix_periodic_Tri_weighted(NumOfLuminal,w1_Delta,w2_Delta);
AdjMatrixVN_type2 = bilayer_Adj_Matrix_periodic_Tri_weighted(NumOfLuminal,w1_Ecad,w2_Ecad);




NumOfBasal = length(AdjMatrixVN(1,:)) - NumOfLuminal;


% ODE systems for the cells
a1 = 0.01;
k1 = 2;

b1 = 50;
h1 = 1;

a2 = 10;
k2 = 1;

b2 = 100;
h2 = 1;




paras = [a1,b1,a2,b2,k1,k2,h1,h2];

%FIND THE HSS BEFORE THIS



N_ss = fsolve(@(N) HSS_notch(N,a1,a2,b1,b2,k1,k2,h1,h2), 0.5)

E_ss = (N_ss.^k2)./(a2 + N_ss.^k2)

D_ss = 1./(1 + b2.*(N_ss.^h2))





lum_IC = repmat([ min(N_ss+ N_ss*perturbation,1),min(E_ss+E_ss*perturbation,1) ,max(D_ss-D_ss*perturbation,0)],1,NumOfLuminal) ;
bas_IC = repmat([max(N_ss - N_ss*perturbation,0) , max(E_ss - E_ss*perturbation,0) , min(D_ss + D_ss*perturbation,1)  ],1,NumOfBasal);

IC = [lum_IC,  bas_IC  ];


%IC = rand(1,3*length(AdjMatrixVN(1,:)));

time = [0,endTime];

num_state_variables = 3*length(AdjMatrixVN(:,1));

[t,y] = ode15s(@(t,y) Ecad_quotient_bilayer(t,y,paras,AdjMatrixVN_type2 ,AdjMatrixVN,num_state_variables), time,IC);

% plotting the cells

CellPositionsLum = [];
CellPositionsBas = [];


for i = 1:NumOfLuminal
    
    CellPositionsLum = [CellPositionsLum; [cos( ((i - 1)*2*pi)/NumOfLuminal),sin( ((i - 1)*2*pi/NumOfLuminal))]   ];
end

for i = 1:NumOfBasal
    
    CellPositionsBas = [CellPositionsBas; [cos(rotationConstant+ ((i - 1)*2*pi)/NumOfBasal),sin(rotationConstant+ ((i - 1)*2*pi/NumOfBasal))]   ];
end


CellPositionsLumScaled = LumRadius*CellPositionsLum;
CellPositionsBasScaled = BasRadius*CellPositionsBas;


%always luminal cells first
FullCells = [CellPositionsLumScaled;CellPositionsBasScaled];

Colors  = lines(8);

Adj1Xedges1 = [];
Adj1Xedges2 = [];
Adj1Yedges1 = [];
Adj1Yedges2 = [];

Adj2Xedges1 = [];
Adj2Xedges2 = [];
Adj2Yedges1 = [];
Adj2Yedges2 = [];


CheckForOverlap1 = [0,0];
CheckForOverlap2 = [0,0];
for i = 1:length(AdjMatrixVN(:,1))
    for j = 1:length(AdjMatrixVN(:,1))-1
        
        
        if AdjMatrixVN_type2(i,j) ~=0
            
            plotFlag2 = 0;
            for k = 1:length(CheckForOverlap2(:,1))
                c1 = [i,j];
                c2 = [j,i];
                if (sum(c1 == CheckForOverlap2(k,:))==2 || sum(c2 == CheckForOverlap2(k,:))==2 )
                    plotFlag2 = 1;
                end
            end
            
            if plotFlag2 == 0
                
                v12 = [FullCells(i,1);  FullCells(i,2)   ];
                v22 = [FullCells(j,1);  FullCells(j,2)   ];
                n2 = (v22 - v12)./(norm(v22-v12,2));
                v2 = [-n2(2);n2(1)];
                
                if sum([NumOfLuminal-1,NumOfLuminal] == [j,i]) ==2
                    shift12 = v12 + projectionDistance.*v2;
                    shift22 = v22 + projectionDistance.*v2;
                elseif sum([1,NumOfLuminal] == [i,j]) ==2
                    shift12 = v12 + projectionDistance.*v2;
                    shift22 = v22 + projectionDistance.*v2;
                elseif sum([2*NumOfLuminal-1,2*NumOfLuminal] == [j,i]) == 2
                    shift12 = v12 + projectionDistance.*v2;
                    shift22 = v22 + projectionDistance.*v2;
                elseif sum([1,2.*NumOfLuminal] == [i,j]) ==2
                    shift12 = v12 + projectionDistance.*v2;
                    shift22 = v22 + projectionDistance.*v2;
                elseif sum([2.*NumOfLuminal,1] == [i,j]) ==2
                    shift12 = v12 + projectionDistance.*v2;
                    shift22 = v22 + projectionDistance.*v2;
                else
                    shift12 = v12 - projectionDistance.*v2;
                    shift22 = v22 - projectionDistance.*v2;
                    
                end

                Adj2Xedges1 = [Adj2Xedges1,shift12(1)];
                Adj2Xedges2 = [Adj2Xedges2,shift22(1)];
                
                Adj2Yedges1 = [Adj2Yedges1,shift12(2)];
                Adj2Yedges2 = [Adj2Yedges2,shift22(2)];
                
                
                CheckForOverlap2 = [CheckForOverlap2; [i,j]];
                
            end
            
            
        end
        
        if AdjMatrixVN(i,j) ~=0
            
            plotFlag1 = 0;
            for k = 1:length(CheckForOverlap1(:,1))
                c1 = [i,j];
                c2 = [j,i];
                if (sum(c1 == CheckForOverlap1(k,:))==2 || sum(c2 == CheckForOverlap1(k,:))==2 )
                    plotFlag1 = 1;
                end
            end
            
            if plotFlag1 == 0
                
                v11 = [FullCells(i,1);  FullCells(i,2)   ];
                v21 = [FullCells(j,1);  FullCells(j,2)   ];
                n1 = (v21 - v11)./(norm(v21-v11,2));
                v1 = [-n1(2);n1(1)];
                
                if sum([NumOfLuminal-1,NumOfLuminal] == [j,i]) ==2
                    shift11 = v11 - projectionDistance.*v1;
                    shift21 = v21 - projectionDistance.*v1;
                elseif sum([1,NumOfLuminal] == [i,j]) ==2
                    shift11 = v11 - projectionDistance.*v1;
                    shift21 = v21 - projectionDistance.*v1;
                    
                elseif sum([2*NumOfLuminal-1,2*NumOfLuminal] == [j,i]) == 2
                    shift11 = v11 - projectionDistance.*v1;
                    shift21 = v21 - projectionDistance.*v1;
                elseif sum([1,2.*NumOfLuminal] == [i,j]) ==2
                    shift11 = v11 - projectionDistance.*v1;
                    shift21 = v21 - projectionDistance.*v1;
                elseif sum([2.*NumOfLuminal,1] == [i,j]) ==2
                    shift11 = v11 - projectionDistance.*v1;
                    shift21 = v21 - projectionDistance.*v1;
                else
                    shift11 = v11 + projectionDistance.*v1;
                    shift21 = v21 + projectionDistance.*v1;
                    
                end

                Adj1Xedges1 = [Adj1Xedges1,shift11(1)];
                Adj1Xedges2 = [Adj1Xedges2,shift21(1)];
                
                Adj1Yedges1 = [Adj1Yedges1,shift11(2)];
                Adj1Yedges2 = [Adj1Yedges2,shift21(2)];
                
                
                CheckForOverlap1 = [CheckForOverlap1; [i,j]];
                
            end
            
        end
        
    end
end

timeStep = 1;

median(divisors(length(t)));
dt = 1:timeStep:length(t);

MeanNotchActivationLum = mean(y(:,1:3:3*NumOfLuminal-2),2);
MeanNotchActivationBas = mean(y(:,3*(NumOfLuminal+1)-2:3:end-2),2);

MeanEcadActivationLum = mean(y(:,2:3:3*NumOfLuminal-1),2);
MeanEcadActivationBas = mean(y(:,3*(NumOfLuminal+1)-1:3:end-1),2);

MeanDeltaActivationLum = mean(y(:,3:3:3*NumOfLuminal),2);
MeanDeltaActivationBas = mean(y(:,3*(NumOfLuminal+1):3:end),2);

figure()
colormap parula

for k = 1:length(dt)
    clf
    subplot(3,2,[1 3 5])
    plot( [Adj1Xedges1;  Adj1Xedges2 ]   , [Adj1Yedges1;  Adj1Yedges2 ]  , '-.','linewidth',2, 'color',[0,0,0] )
    hold on
    plot( [Adj2Xedges1;  Adj2Xedges2 ]   , [Adj2Yedges1;  Adj2Yedges2  ]  , '-','linewidth',2, 'color',[0,0,0]   )
    
    scatter(FullCells(:,1),FullCells(:,2),450,y(dt(k),1:3:end-2),'filled','linewidth',1.5,'MarkerEdgeColor','k')
    axis equal
    axis off
    title(strcat("$t = $ ",num2str(round(t(dt(k)),2))), 'fontsize',30)
    ylim([-(BasRadius + 0.1*BasRadius),(BasRadius + 0.1*BasRadius)])
    xlim([-(BasRadius + 0.1*BasRadius),(BasRadius + 0.1*BasRadius)])
    c = colorbar('location','south');
    caxis([0,1])
    set(c,'TickLabelInterpreter','latex')
    %c.Label.String = 'Notch';
    c.Position = c.Position +[0,-0.15,0,0];
    
    
    subplot(3,2,2)
    plot(t(1:dt(k)),MeanNotchActivationLum(1:dt(k))  )
        hold on
        plot(t(1:dt(k)),MeanNotchActivationBas(1:dt(k))  )
        legend("Lumial cells",  "Basal cells")
        xlabel("Time")
        ylabel("Notch")
        xlim([0,t(end)])
        ylim([0,1])
        box off
        
        subplot(3,2,4)
        plot(t(1:dt(k)),MeanEcadActivationLum(1:dt(k))  )
        hold on
        plot(t(1:dt(k)),MeanEcadActivationBas(1:dt(k))  )
        xlabel("Time")
        ylabel("Ecad")
        xlim([0,t(end)])
        ylim([0,0.1])
        box off
        
        subplot(3,2,6)
        plot(t(1:dt(k)),MeanDeltaActivationLum(1:dt(k))  )
        hold on
        plot(t(1:dt(k)),MeanDeltaActivationBas(1:dt(k))  )
        xlabel("Time")
        ylabel("Delta")
        xlim([0,t(end)])
        ylim([0,0.3])
        box off
        
         drawnow
end


%% end plot
scatter(FullCells(:,1),FullCells(:,2),800,y(end,1:3:end-2),'filled','linewidth',1.5,'MarkerEdgeColor','k')
axis equal
axis off
title(strcat("$t = $ ",num2str(round(t(end),2))), 'fontsize',30)
ylim([-(BasRadius + 0.1*BasRadius),(BasRadius + 0.1*BasRadius)])
xlim([-(BasRadius + 0.1*BasRadius),(BasRadius + 0.1*BasRadius)])
c = colorbar('location','south');
caxis([0,1])
set(c,'TickLabelInterpreter','latex')
%c.Label.String = 'Notch';
c.Position = c.Position +[0,-0.1,0,0];


%% plot relative difference in  each layer over time


close all
figure()
plot(t, MeanNotchActivationLum);
hold on
plot(t, MeanNotchActivationBas);








function out = HSS_notch(N,a1,a2,b1,b2,k1,k2,h1,h2)

out = Ecad_N(N,(N.^k2)./(a2 + N.^k2),(N.^k2)./(a2 + N.^k2),1./(1 + b2.*(N.^h2)),a1,k1,b1,h1);


end
