% This MATLAB file demonstrates the implementation of the Alternating Direction Implicit 
% (ADI) method, a numerical technique used to solve partial differential equations (PDEs), 
% specifically for 2D heat equations. The ADI method is particularly effective for 
% problems where standard explicit methods would require very small time steps to maintain 
% stability.


function u = adiMethod(m, g, u_0)
    % Alternating Direction Implicit (ADI) Method
    % Solves 2D heat equations using the ADI method, a finite difference method.
    % Inputs:
    % m - Grid size for spatial discretization
    % g - Function handle representing the boundary conditions
    % u_0 - Initial condition function handle
    % Outputs:
    % u - Solution matrix of the 2D heat equation
    % Reference: https://en.wikipedia.org/wiki/Alternating_direction_implicit_method

    h = 1/m;  
    k = h;
    x = 0:h:1;
    mu = k/h^2;
    A_f = diag(-2*ones(1,m+1)) + diag(ones(1,m),1) + diag(ones(1,m),-1);
    A_s = diag(-2*ones(1,m-1)) + diag(ones(1,m-2),1) + diag(ones(1,m-2),-1);  
    B_r = sparse(kron(eye(m-1), A_f)); 
    B_l = sparse(kron(eye(m-1), A_s));
    A_l = sparse(kron(A_s, eye(m-1))); 
    A_r = sparse(kron(A_f, eye(m-1))); 
    [x_grid, y_grid ] = meshgrid(x,x); 
    u_in = u_0(x_grid, y_grid); 
    u = zeros((m+1)^2, m+1); 
    u(:, 1) = reshape(u_in, [], 1);

    % [ADI Method Implementation...]

    u = reshape(u, [m+1, m+1, m+1]); 
end
