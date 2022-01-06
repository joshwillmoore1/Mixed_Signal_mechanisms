function f = Ecad_N(N,E,U1,U2,a1,k1,b1,h1,b2,h2)

    f = ((U2.^k1)./(a1 + U2.^k1)).*( 1./(1 + b1.*((U1).^h1))).*( 1./(1 + b2.*((E).^h2))) - N;
    
end