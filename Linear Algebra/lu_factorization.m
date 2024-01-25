function y = LowerSolver(L, x)
    % Lower Triangular Matrix Solver
    % Solves Ly = x for y, where L is a lower triangular matrix.
    % Inputs:
    % L - Lower triangular matrix (m x n)
    % x - Right-hand side vector (m x 1)
    % Output:
    % y - Solution vector (m x 1)

    [m, n] = size(L);
    y = zeros(m, 1);
    y(1) = x(1) / L(1, 1);
    for i = 2:m
        y(i) = (x(i) - dot(L(i, 1:i-1), y(1:i-1))) / L(i, i);
    end
end

function y = UpperSolver(U, x)
    % Upper Triangular Matrix Solver
    % Solves Uy = x for y, where U is an upper triangular matrix.
    % Inputs:
    % U - Upper triangular matrix (m x n)
    % x - Right-hand side vector (m x 1)
    % Output:
    % y - Solution vector (m x 1)

    [m, n] = size(U);
    y = zeros(m, 1);
    y(m) = x(m) / U(m, m);
    for i = m-1:-1:1
        y(i) = (x(i) - dot(U(i, i+1:end), y(i+1:end))) / U(i, i);
    end
end



function [L, U, P] = LU_Factorization(A)
    % LU Factorization with Partial Pivoting
    % Performs LU factorization of a matrix A using partial pivoting.
    % Inputs:
    % A - Matrix to be factorized (m x n)
    % Outputs:
    % L - Lower triangular matrix (n x n)
    % U - Upper triangular matrix (n x n)
    % P - Permutation matrix (n x n)

    [m, n] = size(A);
    p = eye(n);
    for k = 1:n-1
        [A, p] = pivot(A, k, p); 
        for i = k+1:n
           A(i,k) = A(i,k) / A(k,k); 
           A(i, k+1:n) = A(i, k+1:n) - A(i,k) * A(k, k+1:n);
        end
    end
    P = p; 
    L = tril(A, -1) + eye(n);
    U = triu(A);
end

function [A, p] = pivot(A, k, p)
    % Pivot Function for LU Factorization
    % Swaps rows in matrix A and permutation matrix p for LU factorization.
    % Inputs:
    % A - Matrix to be factorized
    % k - Current pivot index
    % p - Current permutation matrix
    % Outputs:
    % A - Matrix after row swapping
    % p - Updated permutation matrix

    [m, n] = size(A);
    big = k; 
    for i = k+1:n
        if abs(A(i, k)) > abs(A(big, k))
            big = i; 
        end
    end
    
    temp = A(k, :); 
    A(k, :) = A(big, :); 
    A(big, :) = temp; 
    p = switch_permutation(p, k, big);
end

function p = switch_permutation(p, i, j)
    % Switch Permutation
    % Swaps columns i and j in a permutation matrix p.
    % Inputs:
    % p - Permutation matrix.
    % i, j - Indices of the columns to be swapped.
    % Output:
    % p - Modified permutation matrix after swapping columns.

    temp = p(:, i);
    p(:, i) = p(:, j); 
    p(:, j) = temp;
end

