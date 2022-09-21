
%%% SBP-Projection solvers for simplification
%%% of Maxwell's equations in 1D with an interface
%%% Calls the function RK4.m

clear all; close all;

% From project description
eps_l = 2;
eps_r = 1;
mu = 1;
r = 0.1;
g_l = 0;
g_r = 0;
A = [0 1; 1 0];
C_l = [eps_l 0; 0 mu];
C_r = [eps_r 0; 0 mu];
T_end = 1.1; % T = 1.1 for reflecton without interference

% Setting the domain
x_l = -1;  x_I = 0;  x_r = 1;
x_0 = -0.5;
len = x_I-x_l;
gridpnts = ([61 101 161 201 501 601]-1)./2+1;
h = len./(gridpnts-1);

% Chooses component from vector
% Forming vector as
% u^(1,l)_1... u^(1,l)_m u^(2,l)_1... u^(2,l)_m u^(1,r)_1... u^(1,r)_m etc
e_1_left = [1 0 0 0];       %u^(1) left domain
e_2_left = [0 1 0 0];       %u^(2) left domain
e_1_right = [0 0 1 0];      %u^(1) right domain
e_2_right = [0 0 0 1];      %u^(2) right domain

n = input('1 for Characteristic \n2 for Dirichlet:          ');
disp(" ");
for i=1:length(gridpnts)
    
    m = gridpnts(i);
    
    xl = x_l:h(i):x_I;
    xr = x_I:h(i):x_r;
    x = [xl xr];
    
    [H, HI, D1, e_1, e_m] = SPB4_BV3(m,h(i));
    [theta_1, theta_end2] = theta(xl,x_0,0,r);
    u_initial = [theta_end2-theta_1 theta_1+theta_end2 zeros(1,2*m)];
    
    % Constructs bonudary-operators dep. on BC.
    switch n
        case 1
            BC = "Characteristic";
            L_left = kron(e_1_left, e_1') - kron(e_2_left, e_1');
            L_int_1 = kron(e_1_left, e_m') - kron(e_1_right, e_1');
            L_int_2 = kron(e_2_right, e_1')-kron(e_2_left, e_m');
            L_right = kron(e_2_right, e_m') + kron(e_2_right, e_m');
            
        case 2
            BC = "Dirichlet";
            L_left = kron(e_1_left, e_1');
            L_int_1 = kron(e_1_left, e_m') - kron(e_1_right, e_1');
            L_int_2 = kron(e_2_left, e_m') - kron(e_2_right, e_1');
            L_right = kron(e_2_right, e_m');
    end
    
    L = [L_left; L_int_1; L_int_2; L_right]';
    
    % Constructs C-matrix
    C = [kron(C_l, eye(m)) zeros(size(kron(C_l, eye(m))));
        zeros(size(kron(C_r, eye(m)))) kron(C_r, eye(m))];
    
    % Projection method
    HI_P = kron(eye(4), HI);
    temp = (L'*HI_P*L)\eye(4);
    P = eye(4*m)-HI_P*L*temp*L';
    D_bar = C\(kron(eye(2), kron(A, D1)));
    Av = P*D_bar*P;
    
    % Time stepping
    eigAv = eig(Av);
    CFL = 2.8/max(abs(eigAv));
    k = (CFL/5)*h(i);
    v0 = u_initial';
    v = RK4(Av, v0, x, T_end, k, m);
    
    % Analytical values
    eta_l = sqrt(eps_l);
    eta_r = sqrt(eps_r);
    T = 2*eta_l/(eta_l+eta_r);       %transmission
    R = (eta_l-eta_r)/(eta_l+eta_r); %reflection
    E_I = 1;                         %height incident wave
    E_T = T*E_I;                     %       transmitted
    E_R = R*E_I;                     %       reflected
    
    
    % Numerical values
    % at t =
    E_R_num = abs(min(v(1:m)));
    E_T_num = abs(min(v(2*m+1:3*m)));
    T_num = E_T_num/E_I;
    R_num = E_R_num/E_I;
    
    err_T(i) = abs(T - T_num);
    err_R(i) = abs(R - R_num);
    err(i) = mean(err_T+err_R);
    
    % Plotting stuff
    figure(i)
    %plot(xl, v(1:m), 'b', xr, v(2*m+1:3*m), 'b');
    hold on;
    plot(xl, v(1:m), 'b', xl, v(m+1:2*m), 'r', ...
        xr, v(2*m+1:3*m), 'b', xr, v(3*m+1:end), 'r');
    plot([x_I x_I], [-2 2],'k');
    title(["SBP-Projection with Interface", "T = " + T_end + ", m = "+ (2*m-1)]);
    legend("E", "H");
    xlabel("x");
    ylabel("Amplitude");
end

%Convergence stuff
q = zeros(1,length(m));
for j=1:length(gridpnts)-1
    q(j+1) = (log(err_T(j+1))-log(err_T(j))) / (log(h(j+1))-log(h(j)));
end

slope = polyfit(log(h), log(err_T),1);

figure(length(gridpnts)+i+1)
loglog(h, err_T, '-b', h, h.^(slope(1)),'r');
title(["Error, SBP-Projection"]);
xlabel("Spatial step size, log(h)");
ylabel("log(e_T)");
legend("Error", "h^P, P =" + slope(1),'Location', 'southeast');


