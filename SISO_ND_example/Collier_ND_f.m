function f = Collier_ND_f(N,D,U,mu1,k,a)
    f = -mu1.*N + (U.^k)./(a + U.^k);
end