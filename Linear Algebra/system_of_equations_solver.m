% This MATLAB script is designed to solve linear systems of equations using various fixed 
% point iteration methods. These methods are particularly useful for large systems where 
% direct methods like Gaussian elimination are computationally expensive. The script 
% includes implementations of the Jacobi method, the Gauss-Seidel method, and the 
% Successive Over-Relaxation (SOR) method.

function [x, residual, closeness] = fp_iteration(A, b, R, g, x_0, tol, true)
    % Fixed Point Iteration for Solving Linear Systems
    % Solves Ax = b using fixed point iteration methods (Jacobi, Gauss-Seidel, SOR).
    % Inputs:
    % A - Coefficient matrix
    % b - Right-hand side vector
    % R - Iteration matrix (specific to the method used)
    % g - Constant vector (specific to the method used)
    % x_0 - Initial guess for the solution
    % tol - Tolerance for the stopping criterion
    % true - True solution vector for error calculation
    % Outputs:
    % x - Computed solution vector
    % residual - Residual at each iteration
    % closeness - Closeness to the true solution at each iteration
    % Reference: https://en.wikipedia.org/wiki/Fixed-point_iteration
    residual = norm(A * x_0 - b);
    cur = x_0; 
    i = 1; 
    while residual > tol
       cur(:, i+1) = R * cur(:,i)+ g; 
       residual(i) = norm(A * cur(:, i) - b); 
       closeness(i) = norm(true - cur(:, i)); 
       i = i + 1;
       if i > 20000
           break
       end
    end
    x = cur(:, end); 
end

function [R, g] = jacobi(A, b)
    % Jacobi Method
    % Solves linear systems using the Jacobi iterative method.
    % Inputs:
    % A - Coefficient matrix
    % b - Right-hand side vector
    % Outputs:
    % R - Iteration matrix for Jacobi method
    % g - Constant vector for Jacobi method
    % Reference: https://en.wikipedia.org/wiki/Jacobi_method
    D = diag(diag(A));
    L = -1 * tril(A, -1); 
    U = -1 * triu(A, 1); 
    R = inv(D) * (L + U); 
    g = inv(D) * b; 
end

function [R, g] = gs(A, b)
    % Gauss-Seidel Method
    % Solves linear systems using the Gauss-Seidel iterative method.
    % Inputs:
    % A - Coefficient matrix
    % b - Right-hand side vector
    % Outputs:
    % R - Iteration matrix for Gauss-Seidel method
    % g - Constant vector for Gauss-Seidel method
    % Reference: https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
    D = diag(diag(A));
    L = -1 * tril(A, -1); 
    U = -1 * triu(A, 1); 
    R = inv(D - L) * U; 
    g = inv(D - L) * b; 
end

function [R, g] = sor(A, b, w)
    % Successive Over-Relaxation (SOR) Method
    % Solves linear systems using the SOR iterative method.
    % Inputs:
    % A - Coefficient matrix
    % b - Right-hand side vector
    % w - Relaxation factor
    % Outputs:
    % R - Iteration matrix for SOR method
    % g - Constant vector for SOR method
    % Reference: https://en.wikipedia.org/wiki/Successive_over-relaxation
    D = diag(diag(A));
    L = -1 * tril(A, -1); 
    U = -1 * triu(A, 1); 
    R = inv(D - w * L) * ((1-w) * D + w * U); 
    g = w * inv(D - w * L) * b;
end
