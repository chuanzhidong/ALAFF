% This matlab script is adopted to demonstrate the Poisson SOR solutions.
% Code credits are partially from the authors of Advanced Linear
% Algebra:Foundations to Frontiers, Drs Robert van de Geijn, Margaret Myers
% and the TA of the course of Spring 2023, Jeffrey Cochran

close all;
clear all;
clc;

% Set number of iterations to be performed
nk = 100

% Set parameters alpha and beta
alpha = 2;
beta  = 3;

% Set the number of meshpoints so the interior has N x N such points
N = 50;

% Compute the distance between mesh points, in each direction
h = 1/(N+1);

% We will have arrays that capture the boundary as well as the interior
% meshpoints.  As a result, we need those arrays to be of size (N+2) x
% (N+2) so that those indexed 2:N+1, 2:N+1 represent the interior.  

% Compute the x-values at each point i,j, including the boundary
x = h * [ 0:N+1 ];   % Notice this creates a row vector

% Compute the y-values at each point i,j, including the boundary
y = h * [ 0:N+1 ];   % Notice this creates a row vector

% Create an array that captures the load at each point i,j
for i=1:N+2
    for j=1:N+2
        F( i,j ) = ...
            ( alpha^2 + beta^2 ) * pi^2 * sin( alpha * pi * x( i ) ) * sin( beta * pi * y( j ) );
    end
end

% define convergance factor
epsilon = 1e-4;
F_Fro_norm = norm(F, "fro");
epsilon_sq_itr_SOR = zeros(nk, 1);
epsilon_sq_itr_GS = zeros(nk, 1);
% converg_basis = (epsilon*F_Fro_norm)^2;

% Set the initial values at the mesh points.  Use SOR to indicate
% the array for the SOR iteration.  Use GS to indicate the array
% is for the Gauss-Seidel iteration.
SOR = zeros( N+2, N+2 );
GS = zeros( N+2, N+2 );

% define the omiga in SOR solution
omega = 1.8; % we make omiga as a number between 1 and 2: here 1.2, 1.4, 1.6, 1.8
% Perform nk iterations
for k = 1:nk
    k           % print current iteration index 
    
    % update all the interior points (Jacobi iteration)
    for i=2:N+1
        for j=2:N+1
            SOR( i,j ) = omega*( 4*(1/omega-1)*SOR(i, j) + SOR( i, j-1 ) + SOR( i-1, j ) + SOR( i+1, j ) + SOR( i, j+1 ) + h^2 * F( i, j ) ) / 4;
        end
    end 
    
%     subplot( 3, 1, 1 );  % plot in top graph
    subplot( 3, 1, 1 );  % plot in top graph
    mesh( x, y, SOR );
%     axis( [ 0 1 0 1 -1.5 1.5 ]);

    % update all the interior points (Gauss-Seidel iteration)
    for i=2:N+1
        for j=2:N+1
            GS( i,j ) = ( GS( i, j-1 ) + GS( i-1, j ) + GS( i+1, j ) + GS( i, j+1 ) + h^2 * F( i, j ) ) / 4;
        end
    end 
    
%     subplot( 3, 1, 2 );  % plot in bottom graph
    subplot( 2, 1, 2 );  % plot in bottom graph
    mesh( x, y, GS );
    axis( [ 0 1 0 1 -1.5 1.5 ]);
%     
%     subplot( 3, 1, 3);
%     mesh( x, y, GS - SOR );
%     axis( [ 0 1 0 1 -0.05 0.05 ])
    
    % wait to continue to the next iteration
%     next = input( 'press RETURN to continue' );
    
    % compare with the predefined convergence basis
    epsilon_sq_running_SOR = squared_error(SOR,F,h) / F_Fro_norm^2;
    epsilon_sq_itr_SOR(k) = squared_error(SOR,F,h) / F_Fro_norm^2;
    epsilon_sq_itr_GS(k) = squared_error(GS,F,h) / F_Fro_norm^2;
    % record the converged iterations
    if epsilon_sq_running_SOR < epsilon^2
        fprintf('converge at %d', k)
%         break
        converge_itr(k, 1) = k;
    end
end

figure(2)
plot(sqrt(epsilon_sq_itr_SOR), 'r', LineWidth=2)
hold on
plot(sqrt(epsilon_sq_itr_GS), 'b', LineWidth=2)
hold on
yline(epsilon, '--g', LineWidth=2)
legend('SOR', 'GS', '1e-4')
title(['Convergences with omega = ', num2str(omega)])
xlabel('Iterations')
ylabel('Convergence value')
% xlim([0 70])

figure(3)
plot(log(sqrt(epsilon_sq_itr_SOR)), 'r', LineWidth=2)
hold on
plot(log(sqrt(epsilon_sq_itr_GS)), 'b', LineWidth=2)
hold on
yline(log(epsilon), '--g', LineWidth=2)
legend('SOR', 'GS', '1e-4')
title(['Convergences with omega = ', num2str(omega), ', log plot'])
xlabel('Iterations')
ylabel('Convergence value: log')
% xlim([0 70])

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



