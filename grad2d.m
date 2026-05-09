function [grad_x, grad_y] = grad2d(fn, IL, JL, dx, dy, depth)
    if nargin < 6
        depth = 1;
    end
    grad_x = zeros(depth, IL, JL);
    grad_y = zeros(depth, IL, JL);
    for i = 1:IL
        for j = 1:JL
            if i == 1
                grad_x(:,i,j) = (-3*fn(i,j) + 4*fn(i+1,j) - fn(i+2,j)) ...
                    / (2 * dx);
            elseif i == IL
                grad_x(:,i,j) = (3*fn(i,j) - 4*fn(i-1,j) + fn(i-2,j)) ...
                    / (2 * dx);
            else
                grad_x(:,i,j) = (fn(i+1,j) - fn(i-1,j)) / (2 * dx);
            end
        end
    end
    for i = 1:IL
        for j = 1:JL
            if j == 1
                grad_y(:,i,j) = (-3*fn(i,j) + 4*fn(i,j+1) - fn(i,j+2)) ...
                    / (2 * dy);
            elseif j == JL
                grad_y(:,i,j) = (3*fn(i,j) - 4*fn(i,j-1) + fn(i,j-2)) ...
                    / (2 * dy);
            else
                grad_y(:,i,j) = (fn(i,j+1) - fn(i,j-1)) / (2 * dy);
            end
        end
    end
    grad_x = squeeze(grad_x);
    grad_y = squeeze(grad_y);
end