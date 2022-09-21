
%%% SBP-SAT and SBP-Projection solvers for simplification
%%% of Maxwell's equations in 1D, imposing two different
%%% BC:s: characterisitic and Dirichlet. 
%%% Calls the function RK4.m

clear all; close all;

% From project description
eps = 1;
mu = 1;
r = 0.1;
g_l = 0;
g_r = 0;
A = [0 1; 1 0];
C = [eps 0; 0 mu];
T = 1.8;

%Domain properties
%m = 201;
m = [61 101 161 201 401 601];
x_l = -1;    x_r = 1;
len = x_r-x_l;
h = (x_r-x_l)./(m-1);

n = input('1 for Characteristic \n2 for Dirichlet:          ');
disp(" ");
nn = input('1 for SAT \n2 for Projection:          ');
disp(" ");

for i=1:length(m)
    
    x = x_l:h(i):x_r;
    x_0 = 0;            %Initial wave starts at x = 0
    [theta_1, theta_end2] = theta(x,x0, 0,r);
    u_initial = [(theta_end2-theta_1) (theta_1+theta_end2)];
    
    [H, HI, D1, e_1, e_m] = SPB4_BV3(m(i),h(i));
    
    % Chooses component from vector
    % Forming vector as
    % u^(1,l)_1... u^(1,l)_m u^(2,l)_1... u^(2,l)_m
    e_first = [1 0];    %u^(1)
    e_second = [0 1];   %u^(2)
    
    switch n
        % Characteristic BC
        case 1
            BC = "Characteristic";
            tau_l = [-1/2; 1/2];
            tau_r = [-1/2; -1/2];
            L_l = kron(e_first, e_1') - kron(e_second, e_1')-g_l;
            L_r = kron(e_first, e_m') + kron(e_second, e_m')-g_r;
            
        % Dirichlet BC
        case 2
            BC = "Dirichlet";
            tau_l = [0; 1];
            tau_r = [-1; 0];
            L_l = kron(e_first, e_1');
            L_r = kron(e_second, e_m');
    end
    
    switch nn
        % Set SBP-SAT
        case 1
            method = "SAT";
            SAT_l = kron(tau_l, HI*e_1*L_l);
            SAT_r = kron(tau_r, HI*e_m*L_r);
            SBPx = kron(A, D1);
            B = SBPx + SAT_l + SAT_r;
            
        % Set SBP-Projection
        case 2
            method = "Projection";
            L_P = [L_l; L_r]';
            HI_P = kron(eye(2), HI);
            temp = inv(L_P'*HI_P*L_P);
            P = eye(2*m(i))-HI_P*L_P*temp*L_P';
            B = P*kron(A, D1)*P;
    end
    
    % Time stepping
    eigB = eig(B);
    CFL = 2.8/max(abs(eigB));
    k = (CFL/5)*h(i);
    
    %k = 0.023*h(i);
    v0 = u_initial';
    v = RK4(B,v0,x,T,k);
    
    % Calculating analytical solution
    [theta_end1, theta_end2]= theta(x,len-T,r);
    u_exact = [-(theta_end1+theta_end2) -(theta_end1-theta_end2)]';
    
    % Calculating error
    err = u_exact-v;
    L2E(i)= sqrt(err'*kron(eye(2),H)*err);
    
    % Plotting
    figure(i)
    subplot(2,1,1)
    plot(x,u_exact(1:m(i)), x, u_exact(m(i)+1:2*m(i)));
    title(["Analytical, T = " + T, BC + " BC"]);
    xlabel("Domain");
    ylabel("Amplitude");
    
    subplot(2,1,2)
    plot(x,v(1:m(i)), x, v(m(i)+1:2*m(i)));
    title( ["SBP-" + method + ", T = " + T, BC + " BC", "m = " m(i)]);
    xlabel("Domain, x");
    ylabel("Amplitude, [H, E]");
    
    figure(i+length(m))
    plot(eigB, '*');
    title(["Eigenvalues " + method, BC + "-BC, m = " + m(i)]);
    xlabel("Re(\lambda)");
    ylabel("Im(\lambda)");
end


%Convergence stuff
q = zeros(1,length(m));
for j=1:length(m)-1
    q(j+1) = (log(L2E(j+1))-log(L2E(j))) / (log(h(j+1))-log(h(j)));
end

slope = polyfit(log(h), log(L2E),1);

figure(length(m)+i+1)
loglog(h, L2E, 'b', h, h.^(slope(1)),'r');
title(["l^2-Error, SBP-" + method, BC + " BC"]);
xlabel("Spatial step size, log(h)");
ylabel("log(||u-v||_{l^2})");
legend("Error", "h^P, P =" + slope(1),'Location', 'southeast');




