
%%% Calculates the value of theta
function [th1, th2] = theta(x,x0,t,r)
    th1 = exp(-(( (x-x0) -t)./r).^2);
    th2 = -exp(-(( (x-x0) +t)./r).^2);
end