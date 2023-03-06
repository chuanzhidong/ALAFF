% Define a function to calculate the squared error, (||A*U-f||_2)^2
% Here's a short script that computes it, taking U, F and h as arguments
% Credits to Jeffrey Cochran
function squared_error = squared_error(U,F,h)
    squared_error = 0;

    N = size(U,1) - 2;

    for i=2:N+1
        for j=2:N+1
            squared_error = squared_error + ( ...
                4*U( i,j ) - ( U( i, j-1 ) + U( i-1, j ) + U( i+1, j ) + U( i, j+1 ) + h^2 * F( i, j ) ) ...
                )^2;
        end
    end
end
