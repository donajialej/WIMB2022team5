function [coss, sinn] = fourier_functions(x,n)

i=1:1:n;
coss= cos(2*pi*i'*x);
sinn= sin(2*pi*i'*x);


end