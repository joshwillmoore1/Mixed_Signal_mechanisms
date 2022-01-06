function f = Ecad_E(N,E,a2,k2)

    f = (N.^k2)./(a2 + N.^k2) - E;
    
end