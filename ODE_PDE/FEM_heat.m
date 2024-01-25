% This MATLAB file implements a Finite Element Method (FEM) solution for the heat equation. 
% FEM is a numerical method for solving problems of engineering and mathematical physics. 
% The script constructs basis functions and their derivatives, forms the system of linear 
% equations using the assembled stiffness matrix A and load vector b, and solves for the 
% FEM solution. 


function [u_hat, x_i] = femSolve(Ns, n)
    % Finite Element Method (FEM) for solving differential equations
    % Solves a given differential equation using FEM for different element sizes.
    % Inputs:
    % Ns - Array of element sizes for FEM discretization
    % n - Number of discretization points in the domain
    % Outputs:
    % u_hat - Array containing FEM solutions for each element size in Ns
    % x_i - Discretized domain points corresponding to the last element size

    u_hat = cell(length(Ns), 1);

    for k = 1:length(Ns)
        N = Ns(k);
        [basis, phi_prime] = constructBasisFunctions(N, n);

        a = @(x) 2 + sin(2.*pi.*x); 
        f = @(x) (12*pi^2 * cos(2*pi.*x)) ./ (2 + sin(2*pi.*x)).^2; 
        x = 0:1/(n-1):1;
        x_i = 0:1/(N-1):1;

        phi_prime = sparse(phi_prime);
        basis = sparse(basis);
        A = build_A(a(x), phi_prime, N);
        b = build_b(f(x), basis, N);

        u_hat{k} = A\b; 
    end
end

% [Insert definitions of 'constructBasisFunctions', 'build_A', and 'build_b' here...]

function [A] = build_A(a, phi_prime, N) 
    % Build the stiffness matrix A for the FEM solution
    % Inputs:
    % a - Coefficient function in the differential equation
    % phi_prime - Derivatives of basis functions
    % N - Number of elements in the finite element discretization
    % Outputs:
    % A - Stiffness matrix

    A = sparse(zeros(N, N)); 
    for i = 1:N
        for j = 1:N
            if abs(i - j) > 1
                continue 
            end
            A(i, j) = int(phi_prime(i, :), phi_prime(j, :), a);
        end
    end 
end

function [b] = build_b(f, basis, N)
    % Build the load vector b for the FEM solution
    % Inputs:
    % f - Source function in the differential equation
    % basis - Basis functions
    % N - Number of elements in the finite element discretization
    % Outputs:
    % b - Load vector

    b = zeros(N, 1); 
    for i = 1:N
        b(i) = int(f, basis(i, :)); 
    end
end

function [phi_p] = phi_p(i, N, n)
    % Derivative of the i-th basis function
    % Inputs:
    % i - Index of the basis function
    % N - Number of elements in the finite element discretization
    % n - Number of points in the discretization
    % Outputs:
    % phi_p - Derivative of the i-th basis function
    x_i = (i-1)/N;
    x = 0:1/(n-1):1;
    phi_p = zeros(1, n);
    phi_p(x_i - x < 1/N & x_i - x > 0) =  N;
    phi_p(x_i - x < 0 & x_i - x > -1/N) =  -N;
    if i == 1
        phi_p(abs(1 - x) < 1/N) =  N;
    end
end

function [phi] = phi(i, N, n)
    % The i-th basis function
    % Inputs and Outputs are similar to the phi_p function.
    x_i = (i-1)/N;
    x = 0:1/(n-1):1;
    phi = zeros(1, n);
    phi(abs(x_i - x) < 1/N) =  1-N*abs(x_i - x(abs(x_i - x) < 1/N));
    if i == 1
        phi(abs(1 - x) < 1/N) =  1-N*abs(1 - x(abs(1 - x) < 1/N));
    end
end

function [val] = int(x, y, z)
    % Numerical integration of product of functions
    % Inputs:
    % x, y - Functions to be integrated
    % z - Optional weight function
    % Outputs:
    % val - Result of the integration
    if ~exist('z','var')
        % third parameter does not exist, so default it to something
         z = ones(1, 2^16); 
       end
       if ~exist('y','var')
          y = ones(1, 2^16);
       end
       val = mean(x.*y.*z);
   end

