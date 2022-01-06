function dydt = Ecad_quotient_bilayer(t,y,paras,A_g_1,A_g_2,total_num_of_state_vars)

dydt = zeros(total_num_of_state_vars,1);

a1 = paras(1);
b1 = paras(2);
a2 = paras(3);
b2 = paras(4);
b3 = paras(5);
k1 = paras(6);
k2 = paras(7);
h1 = paras(8);
h2 = paras(9);
h3 = paras(10);

N = y(1:3:end-2);
E = y(2:3:end-1);
D = y(3:3:end);

%connect cells via adjacency matrices
U1 = A_g_1*E;
U2 = A_g_2*D;

dydt(1:3:end-2) = Ecad_N(N,E,U1,U2,a1,k1,b1,h1,b2,h2);
dydt(2:3:end-1) = Ecad_E(N,E,a2,k2);
dydt(3:3:end) = Ecad_D(N,D,b3,h3);

end