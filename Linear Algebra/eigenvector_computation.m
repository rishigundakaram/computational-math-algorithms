% This MATLAB script is dedicated to exploring various numerical techniques for computing 
% eigenvalues and eigenvectors of matrices, focusing on QR decomposition methods. It 
% encompasses both the QR algorithm without shift and the QR algorithm with Wilkinson 
% shift. Additionally, it includes methods for transforming matrices into Hessenberg form 
% and the inverse power method for eigenvector calculation.


function [E] = QR_method(C, M)
    % QR_method
    % This function computes the eigenvalues of a matrix using the QR algorithm.
    % Inputs:
    %   C - The input matrix for which eigenvalues are to be computed.
    %   M - The number of iterations for the QR algorithm.
    % Outputs:
    %   approx_eigenvalues - Approximated eigenvalues at the final iteration.

    
    % QR Algorithm for Eigenvalue Computation
    for con = 1:M
        [Q, R] = qr(C, 0); % QR decomposition of C
        C = R * Q; % Update C to be R*Q for the next iteration
        EIG1 = sort(diag(C)); % Extract and sort the diagonal elements (approximated eigenvalues)
    end

    % Output the approximated eigenvalues from the final iteration
    E = diag(C);
end

function [E] = QR_no_shift(A)
    % QR Decomposition without Shift
    % This function computes eigenvalues of a matrix A using the QR algorithm without shift.
    % Inputs:
    % A - Square matrix for which eigenvalues are to be computed
    % Outputs:
    % E - Eigenvalues of matrix A
    % Reference: https://en.wikipedia.org/wiki/QR_algorithm
    eps = 10^-6;
    H = hessenberg(A);
    i = 1;
    while i < 100
        [Q, R] = QR_decomposition_givens(H);
        H = R * Q;
        lower = tril(H);
        if norm(lower, inf) <= eps
            break end
        i = i + 1;
    end
    E = diag(H);
end

function [Q, R] = QR_decomposition_givens(A)
    % QR Decomposition using Givens Rotations
    % Computes QR factorization of a matrix A using Givens rotations.
    % Inputs:
    % A - Matrix to be factorized
    % Outputs:
    % Q - Orthogonal matrix
    % R - Upper triangular matrix
    % Reference: https://en.wikipedia.org/wiki/Givens_rotation
    [~,n] = size(A);
    Q = eye(n);
    for i=1:n-1
        G = givens(A, i);
        Q = G * Q;
        A = G * A;
    end
    Q = Q.';
    R = A;
end

function [E] = QR_wilkinson_shift(A)
    % QR Algorithm with Wilkinson Shift
    % Computes eigenvalues of a matrix A using the QR algorithm with Wilkinson shift.
    % Inputs:
    % A - Square matrix for which eigenvalues are to be computed
    % Outputs:
    % E - Eigenvalues of matrix A
    % Reference: https://en.wikipedia.org/wiki/QR_algorithm#Wilkinson_shift
    eps = 10^-6;
    H = hessenberg(A);
    i = 1;
    [~,n] = size(A);
    while i < 1000
        m = choose_index(H);
        if m ==  0
            break 
        end
        sigma = wilkinson_shift(H(m-1:m, m-1:m));
        [Q, R] = QR_decomposition_givens(H - sigma * eye(n));
        H = R * Q + sigma * eye(n);
        lower = tril(H,-1);
        if max(abs(lower)) <= eps
            break 
        end
        i = i + 1;
    end
    E = diag(H);
end

function sigma = wilkinson_shift(A)
    % Wilkinson Shift
    % Calculates the Wilkinson shift for a 2x2 matrix, which is used in the QR algorithm with Wilkinson shift.
    % This shift strategy helps to accelerate the convergence of the QR algorithm.
    % Inputs:
    % A - A 2x2 submatrix from the bottom right corner of a larger matrix.
    % Outputs:
    % sigma - The computed Wilkinson shift.
    % Reference: https://en.wikipedia.org/wiki/QR_algorithm#Wilkinson_shift

    delta = (A(1,1) - A(2,2)) / 2;
    if delta == 0
        sig = -1;
    else
        sig = sign(delta);
    end
    denom = abs(delta) + sqrt(delta^2 + A(2,1)^2);
    sigma = A(2,2) - sig * A(2,1)^2 / denom;
end


function m = choose_index(A)
    % Choose Index for Partitioning
    % Determines the partitioning index for a matrix in the QR algorithm with Wilkinson shift.
    % This function is used to find the submatrix on which to apply the Wilkinson shift.
    % Inputs:
    % A - Matrix being factorized in the QR algorithm.
    % Outputs:
    % m - The partitioning index where the submatrix starts.
    % Reference: [No direct reference, but related to the concept of matrix partitioning in QR algorithms]

    [~,n] = size(A);
    for i = n:-1:2
       if abs(A(i, i-1)) > 10^-6
           m = i;
           return
       end
    end
end



function eigenvectors = inverse_power_method(A, mu)
    % Inverse Power Method
    % Computes eigenvectors of a matrix A for given eigenvalues mu.
    % Inputs:
    % A - Matrix for which eigenvectors are to be computed
    % mu - Eigenvalues corresponding to the eigenvectors
    % Outputs:
    % eigenvectors - Computed eigenvectors of matrix A
    % Reference: https://en.wikipedia.org/wiki/Inverse_iteration
    [m, ~] = size(A);
    [total, ~] = size(mu);
    eigenvectors = zeros(m, total);
    for i=1:total
       x = 2 * rand(m, 1) - 1;
       x = x / norm(x);
       r = x.' * A * x;
       for j=1:3
          y =  (mu(i) * eye(m) - A)\x;
          x = y / norm(y);
          r = x.'* A * x;
        end
       eigenvectors(:, i) = x;
    end

end

function G = givens(A, i)
    % Givens Rotation
    % Creates a Givens rotation matrix for the ith row of matrix A.
    % Inputs:
    % A - Matrix that requires rotation
    % i - Row index for creating Givens rotation
    % Outputs:
    % G - Givens rotation matrix
    % Reference: https://en.wikipedia.org/wiki/Givens_rotation
    c = 0;
    s = -sign(b);
    if abs(b) > abs(a)
       tau = -a / b;
       s = 1 / sqrt(1 + tau^2);
       c = s * tau;
    else
       tau = -b / a;
       c = 1 / sqrt(1 + tau^2);
       s = c * tau;
    end
    a = A(i, i);
    b = A(i, i+1);
    if b == 0
    c = 1;
        s = 0;
    elseif a == 0
        G = eye(size(A));
        c = A(i, i) / sqrt(A(i, i)^2 + A(i+1, i)^2);
        s = - A(i+1, i) / sqrt(A(i, i)^2 + A(i+1, i)^2);
        G(i:i+1,i:i+1) = [c -1 * s; s c];
    end


function [H] = hessenberg(A)
    % Hessenberg Form
    % Converts a given square matrix A to upper Hessenberg form.
    % Inputs:
    % A - Square matrix to be converted
    % Outputs:
    % H - Upper Hessenberg matrix
    % Reference: https://en.wikipedia.org/wiki/Hessenberg_matrix
    [~,n] = size(A);
    for i=1:n-2
        u = house(A(i+1:n,i));
        A(i+1:n, i:n) = A(i+1:n, i:n)  - 2 * u * u.' * A(i+1:n, i:n);
        A(:, i+1:n) = A(:, i+1:n) - A(:, i+1:n) * 2 * u * u.';
    end
    H = A;
end

function u = house(a)
    % Householder Transformation
    % Generates a Householder vector for a given vector a.
    % Inputs:
    % a - Input vector
    % Outputs:
    % u - Householder vector
    % Reference: https://en.wikipedia.org/wiki/Householder_transformation
    r1 = norm(a);
    r = zeros(size(a));
    r(1) = r1;
    u = (a + sign(a(1))*r) / norm((a + sign(a(1))*r));
end
