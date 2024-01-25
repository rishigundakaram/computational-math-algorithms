function L = lagrange_interpolation(x_i, f, x)
    % Lagrange Interpolation
    % Computes the Lagrange interpolating polynomial for a given set of points and a function.
    % Inputs:
    % x_i - Vector of x-coordinates for interpolation points
    % f - Function handle for the function to be interpolated
    % x - Symbolic variable used for constructing the polynomial
    % Outputs:
    % L - Lagrange interpolating polynomial
    % Reference: https://en.wikipedia.org/wiki/Lagrange_polynomial

    m = length(x_i) - 1;
    total = 0; 
    for k = 1:m+1
        sub = 1;
        for j = 1:m+1
            if j == k
                continue
            end
            sub = sub .* (x - x_i(j)) / (x_i(k) - x_i(j));
        end
        total = total + f(x_i(k)) * sub;
    end
    L = total;
end

function N = newton_interpolation(x_i, f, x)
    % Newton Interpolation Method
    % Computes the Newton interpolating polynomial using divided differences.
    % Inputs:
    % x_i - Vector of x-coordinates for interpolation points
    % f - Function handle for the function to be interpolated
    % x - Symbolic variable used for constructing the polynomial
    % Outputs:
    % N - Newton interpolating polynomial
    % Reference: https://en.wikipedia.org/wiki/Newton_polynomial

    m = length(x_i);
    t = zeros(1, m);
    a = zeros(1, m);
    total = 0;

    for i = 1:m
        t(i) = f(x_i(i));
        for j = i-1:-1:1
            t(j) = (t(j+1) - t(j)) / (x_i(i) - x_i(j));
        end
        a(i) = t(1);
        poly = 1;
        for k = 1:i-1
            poly = poly * (x - x_i(k));
        end
        total = total + a(i) * poly;
    end
    N = total;
end

