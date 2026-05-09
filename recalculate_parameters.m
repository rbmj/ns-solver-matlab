function [U,P,T] = recalculate_parameters(U,P,T,gamma,R,T_amb,P_amb,M_inf)
arguments (Input)
    U (4,:,:) double
    P (:,:) double
    T (:,:) double
    gamma (1,1) double
    R (1,1) double
    T_amb (1,1) double
    P_amb (1,1) double
    M_inf (1,1) double
end
arguments (Output)
    U (4,:,:) double
    P (:,:) double
    T (:,:) double
end
    coder.gpu.kernelfun;
    RHO = 1;
    M = 2;
    N = 3;
    E = 4;

    c_amb = sqrt(gamma*R*T_amb);
    s = size(U);
    IL = s(2);
    JL = s(3);
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
end