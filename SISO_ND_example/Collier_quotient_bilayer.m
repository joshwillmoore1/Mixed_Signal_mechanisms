function dydt = Collier_quotient_bilayer(t,y,paras,A_g,row_num)

dydt = zeros(4*row_num,1);

mu1 = paras(1);
mu2 = paras(2);
k = paras(3);
h = paras(4);
a = paras(5);
b = paras(6);

N = y(1:2:end-1);
D = y(2:2:end);

U = A_g*D;

dydt(1:2:end-1) = Collier_ND_f(N,D,U,mu1,k,a);
dydt(2:2:end) = Collier_ND_g(N,D,U,mu2,h,b);

end