function g = Collier_ND_g(N,D,U,mu2,h,b)
    g = -mu2.*D + 1./(1 + b.*N.^h);
end