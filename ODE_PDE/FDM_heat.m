% This MATLAB file implements different finite difference methods for solving partial 
% differential equations (PDEs), specifically focusing on the 1D heat equation. It includes 
% three numerical schemes: Crank Nicolson, Forward Time Central Space (FTCS), Leapfrog. 
% The script also a crank nicolson solver for backwards heat equation, which is not 
% well-posed. 


function u = crankNicolson(u, h, k, x_int, t_int)
    % Crank-Nicolson Scheme
    % Solves the heat equation using the Crank-Nicolson method, which is a finite difference method.
    % Inputs:
    % u - Initial condition matrix
    % h - Spatial step size
    % k - Temporal step size
    % x_int - Interval for the spatial variable x
    % t_int - Interval for the temporal variable t
    % Outputs:
    % u - Matrix containing the solution of the heat equation
    % Reference: https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method

    nh = (x_int(2) - x_int(1)) / h + 1; 
    nt = (t_int(2) - t_int(1)) / k + 1;
    a = -2;
    b = 1;
    A = diag(a * ones(1, nh)) + diag(b * ones(1, nh - 1), 1) + diag(b * ones(1, nh - 1), -1);
    A = 1/2 * k / h^2 * A;
    for n = 2:nt
        u(:, n) = linsolve((eye(nh) - A), (eye(nh) + A) * u(:, n - 1));
    end
end

function u = crankNicolsonBackwards(u, h, k, x_int, t_int)
    % Crank-Nicolson Scheme for Backwards Heat Equation
    % Solves the backwards heat equation using a modified Crank-Nicolson method.
    % Inputs and Outputs are similar to the crankNicolson function.
    % Note: The backwards heat equation is not well-posed, and this method is expected to be divergent.
    % Reference: [No direct reference, but related to the Crank-Nicolson method]

    nh = (x_int(2) - x_int(1)) / h + 1; 
    nt = (t_int(2) - t_int(1)) / k + 1;
    a = -2;
    b = 1;
    A = diag(a * ones(1, nh)) + diag(b * ones(1, nh - 1), 1) + diag(b * ones(1, nh - 1), -1);
    A = 1/2 * k / h^2 * A;
    for n = 2:nt
        u(:, n) = linsolve((eye(nh) + A), (eye(nh) - A) * u(:, n - 1));
    end
end

function u = leapfrog(u, h, k, x_int, t_int)
    % Leapfrog Scheme
    % Solves the heat equation using the Leapfrog method, a finite difference method.
    % Inputs and Outputs are similar to the oneStepForwardTimeCentralSpace function.
    % The Leapfrog method is known for its instability in certain conditions.
    % Reference: https://en.wikipedia.org/wiki/Leapfrog_integration

    nh = (x_int(2) - x_int(1)) / h + 1; 
    nt = (t_int(2) - t_int(1)) / k + 1;
    a = -2;
    b = 1;
    A = diag(a * ones(1, nh)) + diag(b * ones(1, nh - 1), 1) + diag(b * ones(1, nh - 1), -1); 
    A = 2 * k / h^2 * A;
    for n = 3:nt
        u(:, n) = u(:, n - 2) + A * u(:, n - 1);
    end
end

function u = oneStepForwardTimeCentralSpace(u, h, k, x_int, t_int)
    % One Step Forward Time Central Space (FTCS) Scheme
    % Solves the heat equation using the FTCS method, a finite difference method.
    % Inputs:
    % u - Initial condition matrix
    % h - Spatial step size
    % k - Temporal step size
    % x_int - Interval for the spatial variable x
    % t_int - Interval for the temporal variable t
    % Outputs:
    % u - Matrix containing the solution of the heat equation
    % Reference: https://en.wikipedia.org/wiki/FTCS_scheme

    nh = (x_int(2) - x_int(1)) / h + 1; 
    nt = (t_int(2) - t_int(1)) / k + 1;
    a = -2;
    b = 1;
    A = diag(a * ones(1, nh)) + diag(b * ones(1, nh - 1), 1) + diag(b * ones(1, nh - 1), -1); 
    A = eye(nh) + k / h^2 * A;
    for n = 2:nt
        u(:, n) = A * u(:, n - 1);
    end
end

