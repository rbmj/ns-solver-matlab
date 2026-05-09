% MATLAB translation of Fortran code from Prof Shih
function [X, Y] = two_boundary_nozzle_grid(IL, JL, K, A, B, C, D, E, a1, a2)

    % Basic input checks
    if IL < 3 || JL < 3
        error('IL and JL must both be at least 3.');
    end
    if ~(A < B && B < C)
        error('Require A < B < C.');
    end

    deltaXi = 1.0 / (IL - 1);
    Xi = linspace(0, 1, IL);
    Eta = linspace(0, 1, JL);

    % Throat location in Xi-space
    W = (B - A) / (C - A);

    % Lower (1)/Upper (2) boundary points as a function of Xi
    X1 = (C - A) * Xi + A;
    X2 = X1;
    Y1  = zeros(1, IL);
    Y2  = zeros(1, IL);
    for i = 1:IL
        if Xi(i) < W
            Y1(i) = (D - a1) - a1 * cos(pi * (X1(i) - A) / (B - A));
            Y2(i) = (E + a1) + a1 * cos(pi * (X2(i) - A) / (B - A));
        else
            Y1(i) = (D - a2) + a2 * cos(pi * (X1(i) - B) / (C - B));
            Y2(i) = (E + a2) - a2 * cos(pi * (X2(i) - B) / (C - B));
        end
    end
    
    % Hermite basis as a function of eta
    H1 =  2 * Eta.^3 - 3 * Eta.^2 + 1;
    H2 = -2 * Eta.^3 + 3 * Eta.^2;
    H3 =      Eta.^3 - 2 * Eta.^2 + Eta;
    H4 =      Eta.^3 -     Eta.^2;

    dX1 = zeros(1, IL);
    dY1 = zeros(1, IL);
    dX2 = zeros(1, IL);
    dY2 = zeros(1, IL);

    for i = 2:(IL-1)
        dX1(i) = -K * 0.5 * (Y2(i) - Y1(i)) * (Y1(i+1) - Y1(i-1)) / deltaXi;
        dX2(i) = -K * 0.5 * (Y2(i) - Y1(i)) * (Y2(i+1) - Y2(i-1)) / deltaXi;
    end

    for i = 1:IL
        dY1(i) = K * (Y2(i) - Y1(i)) * (C - A);
        dY2(i) = dY1(i);
    end

    % Construct the grid
    X   = zeros(IL, JL);
    Y   = zeros(IL, JL);
    for j = 1:JL
        for i = 1:IL
            X(i,j) = X1(i) * H1(j) + X2(i) * H2(j) + dX1(i) * H3(j) + dX2(i) * H4(j);
            Y(i,j) = Y1(i) * H1(j) + Y2(i) * H2(j) + dY1(i) * H3(j) + dY2(i) * H4(j);
        end
    end
end