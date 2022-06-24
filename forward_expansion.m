function dy = forward_expansion(t,u,thetas,par, n_u)
N = par(1);
p = par(2);
alpha = par(3);
gammaSR = par(4);
gammaSD = par(5);
gammaVR = par(6);
gammaVD = par(7);
end_time = par(8);

x = t/end_time;
i=1:1:max(n_u);
coss= cos(2*pi*i'*x);
sinn= sin(2*pi*i'*x);

A = u(1,:);
B = u(2,:);


n_u_cum = cumsum(n_u);
[S,dS] = fourier(coss(1:n_u(1)),sinn(1:n_u(1)),A(1:n_u(1)),B(1:n_u(1)));
[V,dV] = fourier(coss(1:n_u(2)),sinn(1:n_u(2)),A(n_u(1)+1:n_u_cum(2)),B(n_u+1:n_u_cum(2)));
[Is,dIs] = fourier(coss(1:n_u(3)),sinn(1:n_u(3)),A(n_u_cum(2)+1:n_u_cum(3)),B(n_u_cum(2)+1:n_u_cum(3)));
[Iv,dIv] = fourier(coss(1:n_u(4)),sinn(1:n_u(4)),A(n_u_cum(3)+1:n_u_cum(4)),B(n_u_cum(3)+1:n_u_cum(4)));
[D,dD] = fourier(coss(1:n_u(5)),sinn(1:n_u(5)),A(n_u_cum(4)+1:n_u_cum(5)),B(n_u_cum(4)+1:n_u_cum(5)));


beta  =  betacomp(t,10,thetas,0,60); 

dSdt = - beta * S* ((Is +Iv) / (N - D)) - p*S - dS;
dVdt = p*S - (1-alpha) *beta * V *( (Is +Iv)/ (N - D) ) -dV;
dIsdt = beta * S* ((Is +Iv) / (N - D)) - (gammaSR + gammaSD) * Is - dIs;
dIvdt = (1-alpha) *beta * V *( (Is +Iv)/ (N - D) )- (gammaVR + gammaVD) * Iv - dIv;
dDdt = gammaVD*Iv + gammaSD*Is -dD;

dy = [dSdt; dVdt; dIsdt; dIvdt; dDdt];
end