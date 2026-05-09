function Y = block_tridiag(L, M, N, Y, Pa, Qa)
    arguments
        L (:,:,:)
        M (:,:,:)
        N (:,:,:)
        Y (:,:)
        Pa (:,:) = 0
        Qa (:,:) = 0
    end

    n = size(M,3);
    m = size(M,1);
    if nargin < 6
        P = zeros(m);
        Q = zeros(m);
    else
        P = Pa;
        Q = Qa;
    end
    if m ~= size(M,2)
        error('Elements of M must be square');
    end
    if any(size(Y) ~= [m n])
        error('Y must be m by n');
    end
    if any(size(L) ~= [m m n-1])
        error('L must be m by m by n-1');
    end
    if any(size(N) ~= [m m n-1])
        error('N must be m by m by n-1');
    end
    if any(size(P) ~= [m m])
        error('P must be m by m');
    end
    if any(size(Q) ~= [m m])
        error('Q must be m by m');
    end

    % First row:
    % M1*x1 + N1*x2 + P*x3 = Y1
    X = M(:,:,1) \ [N(:,:,1), P, Y(:,1)];
    N(:,:,1) = X(:,1:m);
    P = X(:,m+1:2*m);
    Y(:,1) = X(:,end);

    % Eliminate x1 from row 2.
    % This modifies both M2 and N2 because row 1 contains x3.
    M(:,:,2) = M(:,:,2) - L(:,:,1) * N(:,:,1);
    N(:,:,2) = N(:,:,2) - L(:,:,1) * P;
    Y(:,2) = Y(:,2) - L(:,:,1) * Y(:,1);

    % Standard forward elimination up to row n-2
    for i = 2:(n-2)
        X = M(:,:,i) \ [N(:,:,i), Y(:,i)];
        N(:,:,i) = X(:,1:m);
        Y(:,i) = X(:,end);

        M(:,:,i+1) = M(:,:,i+1) - L(:,:,i) * N(:,:,i);
        Y(:,i+1) = Y(:,i+1) - L(:,:,i) * Y(:,i);
    end

    % Last row has:
    % Q*x_{n-2} + L_{n-1}*x_{n-1} + M_n*x_n = Y_n
    %
    % Row n-2 is already normalized, so eliminate x_{n-2}
    % from the last row.
    L(:,:,n-1) = L(:,:,n-1) - Q * N(:,:,n-2);
    Y(:,n) = Y(:,n) - Q * Y(:,n-2);

    % Eliminate x_{n-1} from last row
    X = M(:,:,n-1) \ [N(:,:,n-1), Y(:,n-1)];
    N(:,:,n-1) = X(:,1:m);
    Y(:,n-1) = X(:,end);

    M(:,:,n) = M(:,:,n) - L(:,:,n-1) * N(:,:,n-1);
    Y(:,n) = Y(:,n) - L(:,:,n-1) * Y(:,n-1);

    % Solve last block
    Y(:,n) = M(:,:,n) \ Y(:,n);

    % Standard back substitution down to row 2
    for i = (n-1):-1:2
        Y(:,i) = Y(:,i) - N(:,:,i) * Y(:,i+1);
    end

    % First row still has x3 coupling
    Y(:,1) = Y(:,1) - N(:,:,1) * Y(:,2) - P * Y(:,3);
end