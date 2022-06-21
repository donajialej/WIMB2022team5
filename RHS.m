function dy = RHS(t,y,par,thetas)
beta  =  betacomp(t,10,thetas,0,60); 
N = par(1);
p = par(2);
alpha = par(3);
gammaSR = par(4);
gammaSD = par(5);
gammaVR = par(6);
gammaVD = par(7);
dSdt = - beta * y(1)* ((y(3) +y(4)) / (N - y(5))) - p*y(1);
dVdt = p*y(1) - (1-alpha) *beta * y(2) *( (y(3) +y(4))/ (N - y(5)) ) ;
dIsdt = beta * y(1)* ((y(3) +y(4)) / (N - y(5))) - (gammaSR + gammaSD) * y(3);
dIvdt = (1-alpha) *beta * y(2) *( (y(3) +y(4))/ (N - y(5)) ) - (gammaVR + gammaVD) * y(4);
dDdt = gammaVD*y(4) + gammaSD*y(3) ;

dy = [dSdt; dVdt; dIsdt; dIvdt; dDdt];
end


