%% Initial Setup
clear; close all;
% Assuming dimensions given in m
A = 0;
B = 2;
C = 6;
D = 2;
E = 3;
a1 = 0.35;
a2 = 0.75;

zero_celsius = 273.15; %K
k = 0.026; %W/m-K
gamma = 1.4;
R = 287.05; %J/kg-K
M_inf = 0.2;
P_amb = 101300; %Pa
T_amb = 300; %K
mu = 1.85e-5; %Pa-s
v_inf = M_inf*sqrt(gamma*R*T_amb);

%% Grid Generation
% Plot the C-D nozzle separately to make sure nozzle grid looks right
N = 100;
xrange_1 = linspace(A, B, N/2);
xrange_2 = linspace(B, C, N/2);
xrange = [xrange_1 xrange_2];
curve1_1 = (D - a1) - a1*cos(pi*(xrange_1 - A)/(B-A));
curve1_2 = (D - a2) + a2*cos(pi*(xrange_2 - B)/(C-B));
curve1 = [curve1_1 curve1_2];
curve2_1 = (E + a1) + a1*cos(pi*(xrange_1 - A)/(B-A));
curve2_2 = (E + a2) - a2*cos(pi*(xrange_2 - B)/(C-B));
curve2 = [curve2_1 curve2_2];
% Generate the grid
IL = 100;
JL = 40;
[X, Y] = two_boundary_nozzle_grid(IL, JL, 0.08, A, B, C, D, E, a1, a2);

%% Calculate Metric Coefficients

deltaXi = 1 / (IL - 1);
deltaEta = 1 / (JL - 1);

J = zeros(IL, JL);
Xi_x = zeros(IL, JL);
Xi_y = zeros(IL, JL);
Eta_x = zeros(IL, JL);
Eta_y = zeros(IL, JL);
Y_xi = zeros(IL, JL);
X_xi = zeros(IL, JL);
Y_eta = zeros(IL, JL);
X_eta = zeros(IL, JL);
for i = 1:IL
    for j = 1:JL
        if i == 1
            Y_xi(i,j) = (-3*Y(i,j) + 4*Y(i+1,j) - Y(i+2,j)) / (2 * deltaXi);
            X_xi(i,j) = (-3*X(i,j) + 4*X(i+1,j) - X(i+2,j)) / (2 * deltaXi);
        elseif i == IL
            Y_xi(i,j) = (3*Y(i,j) - 4*Y(i-1,j) + Y(i-2,j)) / (2 * deltaXi);
            X_xi(i,j) = (3*X(i,j) - 4*X(i-1,j) + X(i-2,j)) / (2 * deltaXi);
        else
            Y_xi(i,j) = (Y(i+1,j) - Y(i-1,j)) / (2 * deltaXi);
            X_xi(i,j) = (X(i+1,j) - X(i-1,j)) / (2 * deltaXi);
        end
        if j == 1
            Y_eta(i,j) = (-3*Y(i,j) + 4*Y(i,j+1) - Y(i,j+2)) / (2 * deltaEta);
            X_eta(i,j) = (-3*X(i,j) + 4*X(i,j+1) - X(i,j+2)) / (2 * deltaEta);
        elseif j == JL
            Y_eta(i,j) = (3*Y(i,j) - 4*Y(i,j-1) + Y(i,j-2)) / (2 * deltaEta);
            X_eta(i,j) = (3*X(i,j) - 4*X(i,j-1) + X(i,j-2)) / (2 * deltaEta);
        else
            Y_eta(i,j) = (Y(i,j+1) - Y(i,j-1)) / (2 * deltaEta);
            X_eta(i,j) = (X(i,j+1) - X(i,j-1)) / (2 * deltaEta);
        end
        J(i,j) = X_xi(i,j)*Y_eta(i,j) - Y_xi(i,j)*X_eta(i,j);
        Xi_x(i,j) = Y_eta(i,j) / J(i,j);
        Xi_y(i,j) = -X_eta(i,j) / J(i,j);
        Eta_x(i,j) = -Y_xi(i,j) / J(i,j);
        Eta_y(i,j) = X_xi(i,j) / J(i,j);
    end
end

%% Initialize the Problem
delta_t = 1/6000;
t_final = 10;
c_amb = sqrt(gamma*R*T_amb);
RHO = 1;
M = 2;
N = 3;
E = 4;

rho_amb = P_amb / (R*T_amb);

U = zeros(4, IL, JL);
P = ones(IL, JL) .* P_amb;
T = ones(IL, JL) .* T_amb;
U(1,:,:) = rho_amb;
U(2,:,:) = rho_amb*M_inf*c_amb;
U(3,:,:) = 0;
U(4,:,:) = P_amb / (gamma - 1) + 0.5 * rho_amb * (M_inf * c_amb)^2;

figure('Position', [100, 100, 1600, 1200], 'visible', 'off');
t = 0;
counter = 0;
while t < t_final
    [U,P,T] = recalculate_parameters_mex(U,P,T,gamma,R,T_amb,P_amb,M_inf);
    % Output Current Timestep
    fname = sprintf('out/ns_%05i.png', counter);
    tiledlayout(2,1);
    nexttile;
    s = pcolor(X, Y, P, EdgeColor='none');
    c = colorbar();
    c.Label.String = 'Pressure (Pa)';
    xlim([0 6]);
    ylim([0 5]);
    clim([0 3e5]);
    set(gca, 'color', [0.6 0.6 0.6]);
    colormap('turbo');
    nexttile;
    vel = sqrt(squeeze((U(M,:,:).^2 + U(N,:,:).^2) ./ U(RHO,:,:).^2));
    Mach = vel ./ sqrt(gamma*R*T);
    s = pcolor(X, Y, Mach, EdgeColor='none');
    c = colorbar();
    c.Label.String = 'Mach';
    xlim([0 6]);
    ylim([0 5]);
    clim([0 1.5]);
    set(gca, 'color', [0.6 0.6 0.6]);
    colormap('turbo');
    saveas(gcf, fname);
    U = advance_time_mex(U,P,T,gamma,R,mu,k,T_amb,M_inf,deltaXi,deltaEta,delta_t,Xi_x,Xi_y,Eta_x,Eta_y,J);
    t = t + delta_t;
    counter = counter + 1;
end