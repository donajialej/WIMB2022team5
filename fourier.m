function [fun,dfun] = fourier(coss,sinn,A,B)
i=1:1:length(B);
fun = A(1) + A(2:end)*coss(1:end-1) + B*sinn;
dfun = 2*pi*((i.*B)*coss - (i(2:end).*A(2:end))*sinn(1:end-1));
end