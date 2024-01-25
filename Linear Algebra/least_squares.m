% This MATLAB file contains implementations of various linear algebra algorithms for solving least squares problems. 
% It includes methods such as the pseudoinverse, Classical Gram-Schmidt, Modified Gram-Schmidt with pivoting, 
% Householder transformations with pivoting, and Singular Value Decomposition (SVD). 
% We start with the most basic methods, and then move on to more advanced methods that 
% are more numerically stable and scale better with larger matrices.

function c = leastSquaresQR(A, f)
    % Solves least squares problem using QR Decomposition
    % This method decomposes matrix A using QR decomposition and solves for coefficients c.
    % Inputs:
    % A - Design matrix
    % f - Vector of function values at interpolation points
    % Outputs:
    % c - Coefficients of the least squares fit
    % Reference: https://en.wikipedia.org/wiki/QR_decomposition

    [Q, R] = qr(A, 0); 
    c = R \ (Q' * f);
end

function c = leastSquaresSVD(A, f)
    % Solves least squares problem using Singular Value Decomposition (SVD)
    % This method decomposes matrix A using SVD and solves for coefficients c.
    % Inputs:
    % A - Design matrix
    % f - Vector of function values at interpolation points
    % Outputs:
    % c - Coefficients of the least squares fit
    % Reference: https://en.wikipedia.org/wiki/Singular_value_decomposition

    [U, S, V] = svd(A, 'econ'); 
    c = V * (S \ (U' * f));
end

function c = leastSquaresLU(A, f)
    % Solves least squares problem using LU Decomposition
    % This method decomposes matrix A using LU decomposition and solves for coefficients c.
    % Inputs:
    % A - Design matrix
    % f - Vector of function values at interpolation points
    % Outputs:
    % c - Coefficients of the least squares fit
    % Reference: https://en.wikipedia.org/wiki/LU_decomposition

    [L, U, P] = lu(A);
    opts.LT = true; opts.UT = false;
    b = linsolve(L, P*f, opts);
    opts.LT = false; opts.UT = true;
    c = linsolve(U, b, opts);
end



% Classical Gram-Schmidt process
function [Q, R] = clgs(A)
    % Classical Gram-Schmidt
    % Performs QR factorization of a matrix A using the classical Gram-Schmidt process.
    % Inputs:
    % A - Matrix to be factorized (m x n)
    % Outputs:
    % Q - Orthogonal matrix (m x n)
    % R - Upper triangular matrix (n x n)

    [m, n] = size(A);
    Q = zeros(m, n); 
    R = zeros(n, n); 
    % Initialization of matrix
    R(1, 1) = norm(A(:, 1)); 
    Q(:, 1) = A(:, 1) / R(1, 1);
    for j = 2:n
        q = A(:, j); 
        for i = 1:j-1
            R(i,j) = dot(A(:, j), Q(:, i));  
            q = q - R(i, j) * Q(:, i); 
        end
        R(j, j) = norm(q);
        if R(j, j) == 0
            break 
        else
            Q(:, j) = q / R(j, j);
        end     
    end
end

function [Q, R] = mgs(A)
    % Modified Gram-Schmidt (MGS) Process
    % Performs QR factorization of a matrix A using the Modified Gram-Schmidt process.
    % This method is more numerically stable compared to the classical Gram-Schmidt process.
    % Inputs:
    % A - Matrix to be factorized (m x n)
    % Outputs:
    % Q - Orthogonal matrix (m x n)
    % R - Upper triangular matrix (n x n)

    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    
    for i = 1:n
        R(i, i) = norm(A(:, i));
        if R(i, i) == 0
            break
        else
            Q(:, i) = A(:, i) / R(i, i);
        end

        for j = i+1:n
            R(i, j) = dot(A(:, j), Q(:, i));
            A(:, j) = A(:, j) - R(i, j) * Q(:, i);
        end
    end
end


% Modified Gram-Schmidt with Pivoting
function [Q, R, permutation] = mgscp(A)
    % Modified Gram-Schmidt with Column Pivoting
    % Performs QR factorization of a matrix A using the Modified Gram-Schmidt process with column pivoting.
    % Inputs:
    % A - Matrix to be factorized (m x n)
    % Outputs:
    % Q - Orthogonal matrix (m x n)
    % R - Upper triangular matrix (n x n)
    % permutation - Permutation matrix representing the pivoting of columns

    [m, n] = size(A);
    permutation = eye(n);
    Q = zeros(m, n); 
    R = zeros(n, n); 
    for i = 1:n
        [A, j] = switch_pair(A, i);
        permutation = switch_permutation(permutation, i, j);
        R = switch_row(R, i, j); 
        R(i, i) = norm(A(:, i)); 
        if R(i, i) == 0
           break 
        else
           Q(:, i) = A(:, i) / R(i, i);
        end 
        for j = i+1:n
            R(i, j) = dot(A(:, j), Q(:, i));
            A(:, j) = A(:,j) - R(i, j) * Q(:, i); 
        end
    end
end

% Householder Transformation with Pivoting
function [Q, R, P] = house_qr_pivot(A)
    % Householder Transformations with Pivoting
    % Performs QR factorization of a matrix A using Householder transformations with column pivoting.
    % Inputs:
    % A - Matrix to be factorized (m x n)
    % Outputs:
    % Q - Orthogonal matrix (m x m)
    % R - Upper triangular matrix (m x n)
    % P - Permutation matrix representing the pivoting of columns

    [m, n] = size(A); 
    P = eye(n);
    Q = eye(m);
    R = A; 
    for i = 1:n
        [R, j] = switch_sub(R, i);
        P = switch_permutation(P, i, j);
        if norm(R(i:end, i)) == 0
           break
        else
            u = house(R(i:end, i));
            c = -sign(u(1));
            R(i:end, i:end) = c * (R(i:end, i:end) - 2 * u * u' * R(i:end, i:end));
            Pt = c * (eye(m - i + 1) - 2 * u * u'); 
            if i > 1
                Pt = [zeros(m-i+1, i-1) , Pt];
                Pt = [zeros(i-1, m); Pt];
                for k = 1:i-1
                   Pt(k, k) = 1;
                end
            end
            Q = Q * Pt; 
        end
    end
end

% Helper Functions
function [A, j] = switch_pair(A, i)
    % Switch Pair
    % This function swaps the ith column of matrix A with another column that has the maximum Euclidean norm among remaining columns.
    % Inputs:
    % A - The matrix whose columns are to be swapped.
    % i - The index of the column to swap.
    % Outputs:
    % A - The modified matrix after swapping columns.
    % j - The index of the column that was swapped with the ith column.

    [~, j] = max(vecnorm(A(:, i:end)));
    j = i + j - 1;
    temp = A(:, i); 
    A(:, i) = A(:, j); 
    A(:, j) = temp; 
end

function R = switch_row(R, i, j)
    % Switch Row
    % Swaps rows i and j in the upper triangular part of matrix R.
    % Inputs:
    % R - Upper triangular matrix.
    % i, j - Indices of the rows to be swapped.
    % Output:
    % R - Modified matrix after swapping rows.

    temp = R(1:i-1, i); 
    R(1:i-1, i) = R(1:i-1, j); 
    R(1:i-1, j) = temp;
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

function u = house(a)
    % Householder Vector
    % Computes the Householder vector for a given vector a.
    % Input:
    % a - Input vector.
    % Output:
    % u - Householder vector that can be used to create a Householder transformation.

    r1 = norm(a); 
    r = zeros(size(a)); 
    r(1) = r1;
    u = (a + sign(a(1))*r) / norm((a + sign(a(1))*r)); 
end

function [A, j] = switch_sub(A, i)
    % Switch Sub-Matrix
    % Swaps the ith column of matrix A with another column that has the maximum Euclidean norm among remaining columns in the sub-matrix starting from row i.
    % Inputs:
    % A - The matrix whose columns are to be swapped.
    % i - The starting row index for the sub-matrix.
    % Outputs:
    % A - The modified matrix after swapping columns.
    % j - The index of the column that was swapped with the ith column.

    [~, j] = max(vecnorm(A(i:end, i:end)));
    j = i + j - 1;
    temp = A(:, i); 
    A(:, i) = A(:, j); 
    A(:, j) = temp; 
end


% (Include documentation for each helper function like 'switch_pair', 'switch_row', 'switch_permutation', 'house')

% [Add documentation for 'switch_pair', 'switch_row', 'switch_permutation', and 'house' functions here...]

