function f = Ecad_D(N,D,b3,h3)

    f = 1./(1 + b3.*(N.^h3)) - D;

end