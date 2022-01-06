close all
clc
clear

% parameters global
paras = [1,1,2,1,0.01,100];

nb1 = 2;
nl1 = nb1;

nb2 = 3;
nl2 = 1;

w1 = 0.1;
w2 = 1;

time = [0,1000];

% quotient system
num_of_lum_cells = 1;
A_g = Quotient_bilayer_adj(nb1,nb2,nl1,nl2,w1,w2);
num_state_variables = 4*num_of_lum_cells;

lum_IC = [0.6,0.4];
bas_IC = [0.4,0.6];

IC = [bas_IC,lum_IC];
IC_rand = rand(4,1);

[t,y] = ode15s(@(t,y) Collier_quotient_bilayer(t,y,paras,A_g,num_of_lum_cells), time,IC);

% full system
num_of_lum_cells_large = 8;
A_g_large = twoD_zero_curvature_adjacency_mat(num_of_lum_cells_large,w1,w2);
num_state_variables_large = 8*num_of_lum_cells_large;

IC_large = [repmat(lum_IC,1,3*num_of_lum_cells_large),repmat(bas_IC,1,num_of_lum_cells_large) ];
IC_large_rand = rand(8*num_of_lum_cells_large,1);

[t_large,y_large] = ode15s(@(t,y) Collier_quotient_bilayer(t,y,paras,A_g_large,2*num_of_lum_cells_large), time,IC_large);

% plot nodes on a bilayer
total_num_of_cells = 4*num_of_lum_cells_large;

luminal_locations_x = 1:3*num_of_lum_cells_large;
luminal_locations_y = zeros(1,3*num_of_lum_cells_large);

basal_locations_x = 2:3:(3*num_of_lum_cells_large-1);
basal_locations_y = ones(1,num_of_lum_cells_large);

cell_loc_x = [luminal_locations_x,basal_locations_x];
cell_loc_y = [luminal_locations_y,basal_locations_y];

width = 1.5;
figure()
hold on
for i = 1:total_num_of_cells
    for j = 1:total_num_of_cells
        
        if A_g_large(i,j) ~= 0
            plot([cell_loc_x(i),cell_loc_x(j)],[cell_loc_y(i),cell_loc_y(j)],'Linewidth',2.*width,'color','k'  )
        end
   
    end
end

set(gcf,'GraphicsSmoothing','on')
scatter(cell_loc_x, cell_loc_y,500, y_large(end,1:2:end-1),'filled','MarkerEdgeColor','k','Linewidth',width)
ylim([-0.1,1.1])
box off
axis off

hcb = colorbar;
set(hcb,'TickLabelInterpreter','latex','ticks',[0,1])
caxis([0,1])
colorbar off

% Quotient graph figure
theta = linspace(0,2*pi,100);
r = 0.1;
delta_y = 0.1;
top_loop_x = r.*cos(theta);
top_loop_y = r.*sin(theta) + 0.8 + delta_y ;

bot_loop_x = r.*cos(theta);
bot_loop_y = r.*sin(theta) + 0.2 - delta_y;

figure()
plot([0,0],[0.8,0.2],'Linewidth',2.*width,'color','k')
hold on
plot(top_loop_x,top_loop_y,'Linewidth',2.*width,'color','k')
plot(bot_loop_x,bot_loop_y,'Linewidth',2.*width,'color','k')
scatter([0,0],[0.8,0.2],1000,[y(end,1),y(end,3)],'filled','MarkerEdgeColor','k','Linewidth',width)
hcb = colorbar;
set(hcb,'TickLabelInterpreter','latex','ticks',[0,1])
caxis([0,1])
colorbar off
box off
axis off
axis equal
ylim([-0.1,1.1])
xlim([-(r+0.1),(r+0.1)])