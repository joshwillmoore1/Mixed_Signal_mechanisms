close all
clc
clear

% solve for HSS
NEDparameters;
NSS = fsolve(@(N) HSS_notch(N,a1,a2,b1,b2,k1,k2,h1,h2,b3,h3), 0.5);
ESS = f_i(NSS,a2,k2);
DSS = g_i(NSS,b3,h3);


% compute the derivative of the transfer functions
DT = Trans_derivative(NSS,ESS,DSS,a1,a2,b1,b2,b3,k1,k2,h1,h2,h3);

%sweep space
w1_Delta = linspace(0,0.1,600);
w1_Ecad = linspace(0,1.5,600);

InstabilityID_vn = [];% zeros(length(w1_Delta),length(w2_Ecad));
InstabilityID_tri = [];

%connectivity parameters
n1_d = 2;
n2_d_vn = 1;

n1_e = 2;
n2_e_vn  = 1;


n2_d_tri = 2;
n2_e_tri = 2;

%check HSS instability condition
for i = 1:length(w1_Ecad)
    
   for j = 1:length(w1_Delta)
       
       
       a1_vn = (n1_e.*w1_Ecad(i))./(n1_e.*w1_Ecad(i) + n2_e_vn);
       b1_vn = (n1_e.*w1_Ecad(i))./(n1_e.*w1_Ecad(i) + n2_e_vn);
       
       a2_vn = (n1_d.*w1_Delta(j))./(n1_d.*w1_Delta(j) + n2_d_vn);
       b2_vn = (n1_d.*w1_Delta(j))./(n1_d.*w1_Delta(j) + n2_d_vn);
       
       
       a1_tri = (n1_e.*w1_Ecad(i))./(n1_e.*w1_Ecad(i) + n2_e_tri);
       b1_tri = (n1_e.*w1_Ecad(i))./(n1_e.*w1_Ecad(i) + n2_e_tri);
       
       a2_tri = (n1_d.*w1_Delta(j))./(n1_d.*w1_Delta(j) + n2_d_tri);
       b2_tri = (n1_d.*w1_Delta(j))./(n1_d.*w1_Delta(j) + n2_d_tri);
              
       
       LamMat_vn = diag([a1_vn+b1_vn-1,a2_vn+b2_vn-1]);
       LamMat_tri = diag([a1_tri+b1_tri-1,a2_tri+b2_tri-1]);
       
       instab_cond_vn =  trace(LamMat_vn*DT) - det(LamMat_vn*DT);
       instab_cond_tri =  trace(LamMat_tri*DT) - det(LamMat_tri*DT);
       
       if instab_cond_vn>1
           InstabilityID_vn(i,j) = 1;
           
       else
           InstabilityID_vn(i,j) = 0;
       end
       
       if instab_cond_tri>1
           InstabilityID_tri(i,j) = 1;
           
       else
           InstabilityID_tri(i,j) = 0;
       end
       
       
   end
   
end

%point for example sim
w1Ecad = 0.25;
w1Delta = 0.02;

%plot HSS region
[L1,L2] = meshgrid(w1_Ecad,w1_Delta);

fig=gcf;
fig.Position(3:4)=[0.2,0.3];
colors = [1,1,1;0.6,0.6,0.6];

co = [1 0 0;0 1 0;0 0 1];
set(gca,'nextplot','replacechildren')
set(gca,'colororder',co)

fs = 30;
colormap(colors)
hold on
contourf((L1),(L2),(InstabilityID_vn)','-', 'linewidth',2,'edgecolor',[0,0,0])
scatter(w1Ecad,w1Delta,150,'x','MarkerFaceColor',[1,1,1],'markeredgecolor',[1,1,1],'linewidth',2)
xticks([0,1.5])
yticks([0,0.1])
grid off
view(0,90)
shading interp
box off
ylabel("$w_{1}^{ \left[ N \right] }$",'fontsize',fs)
xlabel("$w_{1}^{ \left[ E \right] }$",'fontsize',fs)
axis square

%% analytical functions
function out = f_i(x,a,k)

out = (x.^k)./(a + x.^k);

end

function out = g_i(x,b,h)

out = 1./(1 + b.*x.^h);

end

function out = f_i_prime(x,a,k)

out = (a.*k.*x.^(k-1))./((a+ x.^k).^2);

end

function out = g_i_prime(x,b,h)

out = -(b.*h.*x.^(h-1))./((1+ b.*x.^h).^2);

end

function out = HSS_notch(N,a1,a2,b1,b2,k1,k2,h1,h2,b3,h3)

out = g_i(f_i(N,a2,k2),b1,h1).*g_i(f_i(N,a2,k2),b2,h2).*f_i(g_i(N,b3,h3),a1,k1) - N;

end

function out = Trans_derivative(NSS,ESS,DSS,a1,a2,b1,b2,b3,k1,k2,h1,h2,h3)

f1 = f_i(DSS,a1,k1);
f2 = f_i(NSS,a2,k2);

g1 = g_i(ESS,b1,h1);
g2 = g_i(ESS,b2,h2);
g3 = g_i(NSS,b3,h3);

f1_prime = f_i_prime(DSS,a1,k1);
f2_prime = f_i_prime(NSS,a2,k2);

g1_prime = g_i_prime(ESS,b1,h1);
g2_prime = g_i_prime(ESS,b2,h2);
g3_prime = g_i_prime(NSS,b3,h3);

DetA = f1.*g1.*g2_prime.*f2_prime - 1;


T11 = f1.*g2.*g1_prime.*f2_prime;
T12 = g1.*g2.*f1_prime.*f2_prime;
T21 = f1.*g2.*g1_prime.*g3_prime;
T22 = g1.*g2.*f1_prime.*g3_prime;

out = (-1./(DetA)).*[T11,T12;T21,T22];
end 


