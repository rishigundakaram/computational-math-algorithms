% QR Algorithm for Eigenvalue Computation
% This script implements the QR algorithm for approximating eigenvalues of a matrix.


function [approx_eigenvalues, errors] = qrEigenvalueComputation(C, M)
    % qrEigenvalueComputation
    % This function computes the eigenvalues of a matrix using the QR algorithm.
    % Inputs:
    %   C - The input matrix for which eigenvalues are to be computed.
    %   M - The number of iterations for the QR algorithm.
    % Outputs:
    %   approx_eigenvalues - Approximated eigenvalues at the final iteration.
    %   errors - L2 norm of the difference between true eigenvalues and approximations per iteration.

    
    % QR Algorithm for Eigenvalue Computation
    for con = 1:M
        [Q, R] = qr(C, 0); % QR decomposition of C
        C = R * Q; % Update C to be R*Q for the next iteration
        EIG1 = sort(diag(C)); % Extract and sort the diagonal elements (approximated eigenvalues)
    end

    % Output the approximated eigenvalues from the final iteration
    approx_eigenvalues = diag(C);
end
