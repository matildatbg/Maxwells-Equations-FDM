
%%% Runge-Kutta 4 for timestepping the solvers
%%% to Maxwell's equations with an Interface.
%%% Uncomment stuff in order to plot solution

function v = RK4(A, v0, x, T, k, m)
x_l = x(1:m);
x_r = x(m+1:end);
temp = v0;

%figure(1);
%plot(x(1:m),temp(1:m),x(m:2*m-1),temp(2*m+1:3*m));

k = T/ceil(T/k);    %Step all the way to T
A = k.*(sparse(A));
t=0;
i = 1;
while t<T
    w1 = A*temp;
    w2 = A*(temp+0.5*w1);
    w3 = A*(temp+0.5*w2);
    w4 = A*(temp+w3);
    w = (w1+2*w2+2*w3+w4)/6;
    temp = (temp+w);     
    %plot(x(1:m),temp(1:m), x(m:2*m-1),temp(2*m+1:3*m));
    %plot([temp(1:m);
    %      temp(2*m+1:3*m)]);
    
    %if (mod(i,15) == 0)
    %    plot(x_l, temp(1:m), 'b', x_l, temp(m+1:2*m), 'r', ...
    %        x_r, temp(2*m+1:3*m), 'b', x_r, temp(3*m+1:end), 'r')
    %    ylim([-2 2])
    %    drawnow;
    %end
    i = i+1;
    t=t+k;
end
v = temp;
end


