function y = betacomp(t,m,theta,a,b)     
       y = 0;
   for j = 1:m
       y = y + theta(j).*leg(j-1,t,a,b);
   end
   
end