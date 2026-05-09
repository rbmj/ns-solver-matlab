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
U(3,:,:) = rho_amb*Y_xi*M_inf*c_amb;
U(4,:,:) = P_amb / (gamma - 1) + 0.5 * rho_amb * (M_inf * c_amb)^2;
U_xi = zeros(4, IL, JL);
U_eta = zeros(4, IL, JL);
T_xi = zeros(IL, JL);
T_eta = zeros(IL, JL);
Jtau_xx = zeros(IL, JL);
Jtau_yy = zeros(IL, JL);
Jtau_xy = zeros(IL, JL);
A = zeros(4,4,IL,JL);
B = zeros(4,4,IL,JL);
specrad_A = zeros(IL,JL);
specrad_B = zeros(IL,JL);
Av = zeros(4,4,IL,JL);
Bv = zeros(4,4,IL,JL);
Av_xi = zeros(4,4,IL,JL);
Bv_xi = zeros(4,4,IL,JL);
Av_eta = zeros(4,4,IL,JL);
Bv_eta = zeros(4,4,IL,JL);
RHS = zeros(4,IL,JL);
DeltaJU = zeros(4,IL,JL);

dPdU = @(Ui) [
    (gamma-1)*(Ui(M)^2+Ui(N)^2)/(2*Ui(RHO)^2)
    -Ui(M)*(gamma-1)/Ui(RHO)
    -Ui(N)*(gamma-1)/Ui(RHO)
    gamma-1
    ];
dTdU = @(Ui) [
    (gamma-1)*(Ui(M)^2+Ui(N)^2-Ui(E)*Ui(RHO))/(R*Ui(RHO)^3)
    -Ui(M)*(gamma-1)/(R*Ui(RHO)^2)
    -Ui(N)*(gamma-1)/(R*Ui(RHO)^2)
    (gamma-1)/(R*Ui(RHO))
    ];
dudU = @(Ui) [-Ui(M)/Ui(RHO)^2; 1/Ui(RHO); 0; 0];
dvdU = @(Ui) [-Ui(N)/Ui(RHO)^2; 0; 1/Ui(RHO); 0];

figure('Position', [100, 100, 800, 600], 'visible', 'off');
t = 0;
counter = 0;
while t < t_final
    % Update Pressure/Temperature
    for i=1:IL
        for j=1:JL
            q = (U(2,i,j)^2 + U(3,i,j)^2) / (2*U(1,i,j));
            P(i,j) = (gamma - 1)*(U(4,i,j) - q);
            T(i,j) = P(i,j) / (U(1,i,j)*R);
        end
    end
    % Enforce Boundary Conditions
    % Inflow
    for j=1:JL
        P(1,j) = 2*P(2,j) - P(3,j);
        T(1,j) = T_amb;
        U(RHO,1,j) = P(1,j) / (R*T_amb);
        U(M,1,j) = M_inf * c_amb * U(RHO,1,j);
        U(N,1,j) = 0;
        U(E,1,j) = P(1,j)/(gamma - 1) + U(M,1,j)^2 / (2*U(RHO,1,j));
    end
    if any(U(RHO,:,:) <= 0, 'all')
        [ii,jj] = find(squeeze(U(RHO,:,:)) <= 0, 1);
        error("Negative density at t=%g, i=%d, j=%d, P=%g, rho=%g, e=%g", ...
            t, ii, jj, P(ii,jj), U(RHO,ii,jj), U(E,ii,jj));
    end
    
    if any(P <= 0, "all")
        [ii,jj] = find(P <= 0, 1);
        error("Negative pressure at t=%g, i=%d, j=%d, P=%g, rho=%g, e=%g", ...
            t, ii, jj, P(ii,jj), U(RHO,ii,jj), U(E,ii,jj));
    end
    % Outflow
    for j=1:JL
        T(end,j) = 2*T(end-1,j) - T(end-2,j);
        P(end,j) = P_amb;
        U(RHO,end,j) = P(end,j) / (R*T(end,j));
        U(M,end,j) = U(RHO,end,j) * ( ...
            2*U(M,end-1,j)/U(RHO,end-1,j) - ...
            U(M,end-2,j)/U(RHO,end-2,j));
        U(N,end,j) = U(RHO,end,j) * ( ...
            2*U(N,end-1,j)/U(RHO,end-1,j) - ...
            U(N,end-2,j)/U(RHO,end-2,j));
        U(E,end,j) = P(end,j)/(gamma-1) + ...
            (U(M,end,j)^2 + U(N,end,j)^2) / (2*U(RHO,end,j));
    end
    % Walls
    for i=1:IL
        T(i,1) = T_amb;
        T(i,end) = T_amb;
        P(i,1) = P(i,2);
        P(i,end) = P(i,end-1);
        U(RHO,i,1) = P(i,1) / (R*T(i,1));
        U(RHO,i,end) = P(i,end) / (R*T(i,end));
        U(M,i,1) = 0;
        U(M,i,end) = 0;
        U(N,i,1) = 0;
        U(N,i,end) = 0;
        U(E,i,1) = P(i,1)/(gamma-1);
        U(E,i,end) = P(i,end)/(gamma-1);
    end
    % Output Current Timestep
    fname = sprintf('out/ns_%05i.png', counter);
    s = pcolor(X, Y, P, EdgeColor='none');
    c = colorbar();
    c.Label.String = 'Pressure (Pa)';
    xlim([0 6]);
    ylim([0 6]);
    clim([0 1e7]);
    set(gca, 'color', [0.6 0.6 0.6]);
    saveas(gcf, fname);
    % Calculate Gradients for Jacobians + Shear Stresses
    [T_xi, T_eta] = grad2d(@(i,j) T(i,j), IL, JL, deltaXi, deltaEta);
    [P_xi, P_eta] = grad2d(@(i,j) P(i,j), IL, JL, deltaXi, deltaEta);
    [U_xi, U_eta] = grad2d(@(i,j) U(:,i,j), IL, JL, deltaXi, deltaEta, 4);
    
    [u_xi, u_eta] = grad2d(@(i,j) U(M,i,j)/U(RHO,i,j), IL, JL, deltaXi, deltaEta);
    [v_xi, v_eta] = grad2d(@(i,j) U(N,i,j)/U(RHO,i,j), IL, JL, deltaXi, deltaEta);
    u_x = Xi_x .* u_xi + Eta_x .* u_eta;
    u_y = Xi_y .* u_xi + Eta_y .* u_eta;
    v_x = Xi_x .* v_xi + Eta_x .* v_eta;
    v_y = Xi_y .* v_xi + Eta_y .* v_eta;
    Jtau_xx = J.*mu.*(4*u_x - 2*v_y)/3;
    Jtau_yy = J.*mu.*(4*v_y - 2*u_x)/3;
    Jtau_xy = J.*mu.*(u_y+v_x);

    T_x = @(i,j) T_xi(i,j)*Xi_x(i,j) + T_eta(i,j)*Eta_x(i,j);
    T_y = @(i,j) T_xi(i,j)*Xi_y(i,j) + T_eta(i,j)*Eta_y(i,j);
   
    % Calculate the A/B Matrices
    for i=1:IL
        for j=1:JL
            A(:,:,i,j) = calcA(P(i,j), U(E,i,j), gamma, U(M,i,j), U(N,i,j), U(RHO,i,j));
            B(:,:,i,j) = calcB(P(i,j), U(E,i,j), gamma, U(M,i,j), U(N,i,j), U(RHO,i,j));
            specrad_A(i,j) = abs(U(M,i,j)/U(RHO,i,j)) + ...
                sqrt(gamma*P(i,j)/U(RHO,i,j));
            specrad_B(i,j) = abs(U(N,i,j)/U(RHO,i,j)) + ...
                sqrt(gamma*P(i,j)/U(RHO,i,j));
            Av(:,:,i,j) = dFv_dU(R, U(E,i,j), U_eta(E,i,j), U_xi(E,i,j), ...
                Eta_x(i,j), Eta_y(i,j), gamma, k, U(M,i,j), ...
                U_eta(M,i,j), U_xi(M,i,j), mu, U(N,i,j), U_eta(N,i,j), ...
                U_xi(N,i,j), U(RHO,i,j), U_eta(RHO,i,j), U_xi(RHO,i,j), ...
                Xi_x(i,j), Xi_y(i,j));
            Bv(:,:,i,j) = dGv_dU(R, U(E,i,j), U_eta(E,i,j), U_xi(E,i,j), ...
                Eta_x(i,j), Eta_y(i,j), gamma, k, U(M,i,j), ...
                U_eta(M,i,j), U_xi(M,i,j), mu, U(N,i,j), U_eta(N,i,j), ...
                U_xi(N,i,j), U(RHO,i,j), U_eta(RHO,i,j), U_xi(RHO,i,j), ...
                Xi_x(i,j), Xi_y(i,j));
            Av_xi(:,:,i,j) = dFv_dU_xi(R, U(E,i,j), gamma, k, U(M,i,j), ...
                mu, U(N,i,j), U(RHO,i,j), Xi_x(i,j), Xi_y(i,j));
            Bv_xi(:,:,i,j) = dGv_dU_xi(R, U(E,i,j), gamma, k, U(M,i,j), ...
                mu, U(N,i,j), U(RHO,i,j), Xi_x(i,j), Xi_y(i,j));
            Av_eta(:,:,i,j) = dFv_dU_eta(R, U(E,i,j), Eta_x(i,j), Eta_y(i,j), ...
                gamma, k, U(M,i,j), mu, U(N,i,j), U(RHO,i,j));
            Bv_eta(:,:,i,j) = dGv_dU_eta(R, U(E,i,j), Eta_x(i,j), Eta_y(i,j), ...
                gamma, k, U(M,i,j), mu, U(N,i,j), U(RHO,i,j));
        end
    end
    % Explicit portion
    energy_term_f = @(i,j) ...
        (U(M,i,j)*Jtau_xx(i,j) + U(N,i,j)*Jtau_xy(i,j))/U(RHO,i,j) + ...
        J(i,j)*k*T_x(i,j);
    energy_term_g = @(i,j) ...
        (U(M,i,j)*Jtau_xy(i,j) + U(N,i,j)*Jtau_yy(i,j))/U(RHO,i,j) + ...
        J(i,j)*k*T_y(i,j);
    [Fhat_v_explicit, ~] = grad2d(@(i,j) ...
        Xi_x(i,j)*[0; Jtau_xx(i,j); Jtau_xy(i,j); energy_term_f(i,j)] + ...
        Xi_y(i,j)*[0; Jtau_xy(i,j); Jtau_yy(i,j); energy_term_g(i,j)], ...
        IL, JL, deltaXi, deltaEta, 4);
    [~, Ghat_v_explicit] = grad2d(@(i,j) ...
        Eta_x(i,j)*[0; Jtau_xx(i,j); Jtau_xy(i,j); energy_term_f(i,j)] + ...
        Eta_y(i,j)*[0; Jtau_xy(i,j); Jtau_yy(i,j); energy_term_g(i,j)], ...
        IL, JL, deltaXi, deltaEta, 4);
    for i=2:IL-1
        for j=2:JL-1
            RHS(:,i,j) = [0;0;0;0];
            % Xi direction
            Aplus_upwind = (A(:,:,i-1,j) + eye(4)*specrad_A(i-1,j))/2;
            Aplus = (A(:,:,i,j) + eye(4)*specrad_A(i,j))/2;
            Aminus = (A(:,:,i,j) - eye(4)*specrad_A(i,j))/2;
            Aminus_downwind = (A(:,:,i+1,j) - eye(4)*specrad_A(i+1,j))/2;

            Bplus_upwind = (B(:,:,i-1,j) + eye(4)*specrad_B(i-1,j))/2;
            Bplus = (B(:,:,i,j) + eye(4)*specrad_B(i,j))/2;
            Bminus = (B(:,:,i,j) - eye(4)*specrad_B(i,j))/2;
            Bminus_downwind = (B(:,:,i+1,j) - eye(4)*specrad_B(i+1,j))/2;

            RHS(:,i,j) = RHS(:,i,j) - Fhat_v_explicit(:,i,j) + (...
                (Xi_x(i,j)*J(i,j)*Aplus*U(:,i,j) ...
                    - Xi_x(i-1,j)*J(i-1,j)*Aplus_upwind*U(:,i-1,j)) ...
                + (Xi_x(i+1,j)*J(i+1,j)*Aminus_downwind*U(:,i+1,j) ...
                    - Xi_x(i,j)*J(i,j)*Aminus*U(:,i,j)) ...
                + (Xi_y(i,j)*J(i,j)*Bplus*U(:,i,j) ...
                    - Xi_y(i-1,j)*J(i-1,j)*Bplus_upwind*U(:,i-1,j)) ...
                + (Xi_y(i+1,j)*J(i+1,j)*Bminus_downwind*U(:,i+1,j) ...
                    - Xi_y(i,j)*J(i,j)*Bminus*U(:,i,j))...
            )/deltaXi;
            % Eta direction
            Aplus_upwind = (A(:,:,i,j-1) + eye(4)*specrad_A(i,j-1))/2;
            Aminus_downwind = (A(:,:,i,j+1) - eye(4)*specrad_A(i,j+1))/2;
            Bplus_upwind = (B(:,:,i,j-1) + eye(4)*specrad_B(i,j-1))/2;
            Bminus_downwind = (B(:,:,i,j+1) - eye(4)*specrad_B(i,j+1))/2;

            RHS(:,i,j) = RHS(:,i,j) - Ghat_v_explicit(:,i,j) + ( ...
                (Eta_x(i,j)*J(i,j)*Aplus*U(:,i,j) ...
                    - Eta_x(i,j-1)*J(i,j-1)*Aplus_upwind*U(:,i,j-1)) ...
                + (Eta_x(i,j+1)*J(i,j+1)*Aminus_downwind*U(:,i,j+1) ...
                    - Eta_x(i,j)*J(i,j)*Aminus*U(:,i,j)) ...
                + (Eta_y(i,j)*J(i,j)*Bplus*U(:,i,j) ...
                    - Eta_y(i,j-1)*J(i,j-1)*Bplus_upwind*U(:,i,j-1)) ...
                + (Eta_y(i,j+1)*J(i,j+1)*Bminus_downwind*U(:,i,j+1) ...
                    - Eta_y(i,j)*J(i,j)*Bminus*U(:,i,j)) ...
            )/deltaEta;
            % Cross derivatives
            [DeltaJU_xi, DeltaJU_eta] = grad2d(@(i,j) DeltaJU(:,i,j)/J(i,j), ...
                IL, JL, deltaXi, deltaEta, 4);
            RHS(:,i,j) = RHS(:,i,j) - ( ...
                Xi_x(i+1,j)*J(i+1,j)*Av_eta(:,:,i+1,j)*DeltaJU_eta(:,i+1,j) ...
                - Xi_x(i-1,j)*J(i-1,j)*Av_eta(:,:,i-1,j)*DeltaJU_eta(:,i-1,j) ...
                + Xi_y(i+1,j)*J(i+1,j)*Bv_eta(:,:,i+1,j)*DeltaJU_eta(:,i+1,j) ...
                - Xi_y(i-1,j)*J(i-1,j)*Bv_eta(:,:,i-1,j)*DeltaJU_eta(:,i-1,j) ...
            )/(4*deltaXi);
            RHS(:,i,j) = RHS(:,i,j) - ( ...
                Eta_x(i,j+1)*J(i,j+1)*Av_xi(:,:,i,j+1)*DeltaJU_xi(:,i,j+1) ...
                - Eta_x(i,j-1)*J(i,j-1)*Av_xi(:,:,i,j-1)*DeltaJU_xi(:,i,j-1) ...
                + Eta_y(i,j+1)*J(i,j+1)*Bv_xi(:,:,i,j+1)*DeltaJU_xi(:,i,j+1) ...
                - Eta_y(i,j-1)*J(i,j-1)*Bv_xi(:,:,i,j-1)*DeltaJU_xi(:,i,j-1) ...
            )/(4*deltaEta);
            % Finally... multiply by -delta t (negative to move to RHS)
            RHS(:,i,j) = RHS(:,i,j) * -delta_t;
        end
    end
    UP = zeros(4,4,JL-1);
    MD = zeros(4,4,JL);
    DN = zeros(4,4,JL-1);
    YY = zeros(4,JL);
    % Solve across the eta direction
    for i=2:IL-1
        for j=2:JL-1
            YY(:,j) = RHS(:,i,j)*(2/delta_t);
            UP(:,:,j-1) = ( ...
                Eta_x(i,j-1)*(A(:,:,i,j-1) + eye(4)*specrad_A(i,j-1)) ...
                + Eta_y(i,j-1)*(B(:,:,i,j-1) + eye(4)*specrad_B(i,j-1)) ...
            )/(-2*deltaEta);
            MD(:,:,j) = ( ...
                Eta_x(i,j)*(A(:,:,i,j) + eye(4)*specrad_A(i,j)) ...
                + Eta_y(i,j)*(B(:,:,i,j) + eye(4)*specrad_B(i,j)) ...
                - Eta_x(i,j)*(A(:,:,i,j) - eye(4)*specrad_A(i,j)) ...
                - Eta_y(i,j)*(B(:,:,i,j) - eye(4)*specrad_B(i,j)) ...
            )/(2*deltaEta);
            DN(:,:,j) = ( ...
                Eta_x(i,j+1)*(A(:,:,i,j+1) - eye(4)*specrad_A(i,j+1))...
                + Eta_y(i,j+1)*(B(:,:,i,j+1) - eye(4)*specrad_B(i,j+1)) ...
            )/(2*deltaEta);
            UP(:,:,j-1) = UP(:,:,j-1) - ( ...
                Eta_x(i,j-1)*Av_eta(:,:,i,j) + ...
                Eta_y(i,j-1)*Bv_eta(:,:,i,j) ...
            )/(deltaEta^2);
            MD(:,:,j) = MD(:,:,j) + 2*(...
                Eta_x(i,j)*Av_eta(:,:,i,j) + ...
                Eta_y(i,j)*Bv_eta(:,:,i,j) ...
            )/(deltaEta^2);
            DN(:,:,j) = DN(:,:,j) - ( ...
                Eta_x(i,j+1)*Av_eta(:,:,i,j) + ...
                Eta_y(i,j+1)*Bv_eta(:,:,i,j) ...
            )/(deltaEta^2);
            dAv = (Av_eta(:,:,i,j+1) - Av_eta(:,:,i,j-1))/(2*deltaEta);
            dBv = (Bv_eta(:,:,i,j+1) - Bv_eta(:,:,i,j-1))/(2*deltaEta);
            UP(:,:,j-1) = UP(:,:,j-1) + (...
                dAv*Eta_x(i,j-1) + dBv*Eta_y(i,j-1) ...
                + Av(:,:,i,j-1)*Eta_x(i,j-1) ...
                + Bv(:,:,i,j-1)*Eta_y(i,j-1) ...
            )/(2*deltaEta);
            DN(:,:,j) = DN(:,:,j) - ( ...
                dAv*Eta_x(i,j+1) + dBv*Eta_y(i,j+1) ...
                + Av(:,:,i,j+1)*Eta_x(i,j+1) ...
                + Bv(:,:,i,j+1)*Eta_y(i,j+1) ...
            )/(2*deltaEta);
            MD(:,:,j) = MD(:,:,j) + eye(4)*2/delta_t;
        end
        % Wall Boundary Conditions
        MD(:,:,1) = [
            0 1 0 0
            0 0 1 0
            dPdU(U(:,i,1)).'
            R*T_amb 0 0 (1-gamma)
        ] ./ J(i,1);
        DN(:,:,1) = [
            0 0 0 0
            0 0 0 0
            -dPdU(U(:,i,2)).'
            0 0 0 0
        ] ./ J(i,2);
        YY(:,1) = [0; 0; 0; 0];
        MD(:,:,end) = [
            0 1 0 0
            0 0 1 0
            dPdU(U(:,i,end)).'
            R*T_amb 0 0 (1-gamma)
        ] ./ J(i,end);
        UP(:,:,end-1) = [
            0 0 0 0
            0 0 0 0
            -dPdU(U(:,i,end-1)).'
            0 0 0 0
        ] ./ J(i,end-1);
        YY(:,end) = [0; 0; 0; 0];
        YY = block_tridiag(UP, MD, DN, YY);
        for j=1:JL
            DeltaJU(:,i,j) = YY(:,j);
        end
    end
    % Now in the Xi direction (DeltaJU is DeltaJU* right now)
    UP = zeros(4,4,IL-1);
    MD = zeros(4,4,IL);
    DN = zeros(4,4,IL-1);
    YY = zeros(4,IL);
    for j=2:JL-1
        for i=2:IL-1
            YY(:,i) = DeltaJU(:,i,j)*(2/delta_t);
            UP(:,:,i-1) = ( ...
                Xi_x(i-1,j)*(A(:,:,i-1,j) + eye(4)*specrad_A(i-1,j)) ...
                + Xi_y(i-1,j)*(B(:,:,i-1,j) + eye(4)*specrad_B(i-1,j)) ...
            )/(-2*deltaXi);
            MD(:,:,i) = ( ...
                Xi_x(i,j)*(A(:,:,i,j) + eye(4)*specrad_A(i,j)) ...
                + Xi_y(i,j)*(B(:,:,i,j) + eye(4)*specrad_B(i,j)) ...
                - Xi_x(i,j)*(A(:,:,i,j) - eye(4)*specrad_A(i,j)) ...
                - Xi_y(i,j)*(B(:,:,i,j) - eye(4)*specrad_B(i,j)) ...
            )/(2*deltaXi);
            DN(:,:,i) = ( ...
                Xi_x(i+1,j)*(A(:,:,i+1,j) - eye(4)*specrad_A(i+1,j))...
                + Xi_y(i+1,j)*(B(:,:,i+1,j) - eye(4)*specrad_B(i+1,j)) ...
            )/(2*deltaXi);
            UP(:,:,i-1) = UP(:,:,i-1) - ( ...
                Xi_x(i-1,j)*Av_xi(:,:,i,j) + ...
                Xi_y(i-1,j)*Bv_xi(:,:,i,j) ...
            )/(deltaXi^2);
            MD(:,:,i) = MD(:,:,i) + 2*(...
                Xi_x(i,j)*Av_xi(:,:,i,j) + ...
                Xi_y(i,j)*Bv_xi(:,:,i,j) ...
            )/(deltaXi^2);
            DN(:,:,i) = DN(:,:,i) - ( ...
                Xi_x(i+1,j)*Av_xi(:,:,i,j) + ...
                Xi_y(i+1,j)*Bv_xi(:,:,i,j) ...
            )/(deltaXi^2);
            dAv = (Av_xi(:,:,i+1,j) - Av_xi(:,:,i-1,j))/(2*deltaXi);
            dBv = (Bv_xi(:,:,i+1,j) - Bv_xi(:,:,i-1,j))/(2*deltaXi);
            UP(:,:,i-1) = UP(:,:,i-1) + (...
                dAv*Xi_x(i-1,j) + dBv*Xi_y(i-1,j) ...
                + Av(:,:,i-1,j)*Xi_x(i-1,j) ...
                + Bv(:,:,i-1,j)*Xi_y(i-1,j) ...
            )/(2*deltaXi);
            DN(:,:,i) = DN(:,:,i) - ( ...
                dAv*Xi_x(i+1,j) + dBv*Xi_y(i+1,j) ...
                + Av(:,:,i+1,j)*Xi_x(i+1,j) ...
                + Bv(:,:,i+1,j)*Xi_y(i+1,j) ...
            )/(2*deltaXi);
            MD(:,:,i) = MD(:,:,i) + eye(4)*2/delta_t;
        end
        % Inflow Boundary Conditions
        MD(:,:,1) = [
            -M_inf*c_amb 1 0 0
            0 0 1 0
            dPdU(U(:,1,j)).'
            dTdU(U(:,1,j)).'
        ] ./ J(1,j);
        DN(:,:,1) = [
            0 0 0 0
            0 0 0 0
            -dPdU(U(:,2,j)).'*2
            0 0 0 0
        ] ./ J(2,j);
        P3 = [
            0 0 0 0
            0 0 0 0
            dPdU(U(:,3,j)).'
            0 0 0 0
        ] ./ J(3,j);
        YY(:,1) = [0; 0; 0; 0];
        % Outflow Boundary Conditions
        MD(:,:,end) = [
            dTdU(U(:,end,j)).'
            dudU(U(:,end,j)).'
            dvdU(U(:,end,j)).'
            dPdU(U(:,end,j)).'
        ] ./ J(end,j);
        UP(:,:,end-1) = [
            -2*dTdU(U(:,end-1,j)).'
            -2*dudU(U(:,end-1,j)).'
            -2*dvdU(U(:,end-1,j)).'
            0 0 0 0
        ] ./ J(end-1,j);
        Q3 = [
            dTdU(U(:,end-2,j)).'
            dudU(U(:,end-2,j)).'
            dvdU(U(:,end-2,j)).'
            0 0 0 0
        ] ./ J(end-2,j);
        YY(:,end) = [0; 0; 0; 0];
        YY = block_tridiag(UP, MD, DN, YY, P3, Q3);
        for i=1:IL
            DeltaJU(:,i,j) = YY(:,i);
        end
    end
    % Update with calculated delta
    for i=2:IL-1
        for j=2:JL-1
            U(:,i,j) = U(:,i,j) + DeltaJU(:,i,j) / J(i,j);
        end
    end
    t = t + delta_t;
    counter = counter + 1;
end