
function [a,b,c,d] = nkcsi(k, h, f)
    % Natural Cubic Spline Interpolation
    % This function computes the natural cubic spline coefficients for a given set of data points.
    % Reference: https://en.wikipedia.org/wiki/Spline_interpolation#Algorithm_to_find_the_interpolating_cubic_spline

    % Inputs:
    %   k - The input matrix for which eigenvalues are to be computed.
    %   h - The number of iterations for the QR algorithm.
    %   f - The number of iterations for the QR algorithm.
    % Outputs:
    %   [a, b, c, d] - coefficients of the natural cubic spline.

    [~, n] = size(k);
    a = arrayfun(f, k);
    A = 4 * h * eye(n) + triu(tril(h*ones(n),1),1) + tril(triu(h*ones(n),-1),-1);
    A(1,1) = h; 
    A(1,2) = 2*h;
    A(n,n) = h;
    A(n,n-1) = 2*h;
    sol = zeros(n,1);
    for i=1:n
       if i == 1
          sol(i) = 1/(2) * (5*(a(2)-a(1)) - (a(3)-a(2))); 
       elseif i == n
           sol(i) = 1/(2) * (a(i-1) - a(i-2) +5*(a(i)-a(i-1)));
       else
          sol(i) = 3/h * (a(i+1)-2 * a(i) + a(i-1)); 
       end
    end
    c = A\sol; 
    b = zeros(n-1,1); 
    d = zeros(n-1,1);
    for i=1:n-1
        b(i) = 1/h * (a(i+1)-a(i)) - h/3 * (2 * c(i) + c(i+1)); 
        d(i) = 1/(3*h)*(c(i+1)-c(i));
    end
    a = a(1:n-1)';
    c = c(1:n-1);
end

function [a,b,c,d] = ncsi(k, h, f)
    % Natural Cubic Spline Interpolation
    % This function computes the natural cubic spline coefficients for a given set of data points,
    % enforcing the condition that the second derivative at the endpoints is zero.
    % Reference: https://en.wikipedia.org/wiki/Spline_interpolation#Algorithm_to_find_the_interpolating_cubic_spline
    % Inputs:
    %   k - The input matrix for which eigenvalues are to be computed.
    %   h - The number of iterations for the QR algorithm.
    %   f - The number of iterations for the QR algorithm.
    % Outputs:
    %   [a, b, c, d] - coefficients of the natural cubic spline.
    [~, n] = size(k);
    a = arrayfun(f, k);
    A = 4 * h * eye(n) + triu(tril(h*ones(n),1),1) + tril(triu(h*ones(n),-1),-1);
    A(1,1) = 1;
    A(n,n) = 1;
    A(1,2) = 0;
    A(n,n-1) = 0;
    sol = zeros(n,1);
    for i=1:n
       if i == 1 || i == n
          sol(i) = 0; 
       else
          sol(i) = 3/h * (a(i+1)-2 * a(i) + a(i-1)); 
       end
    end
    c = A\sol; 
    b = zeros(n-1,1); 
    d = zeros(n-1,1);
    for i=1:n-1
        b(i) = 1/h * (a(i+1)-a(i)) - h/3 * (2 * c(i) + c(i+1)); 
        d(i) = 1/(3*h)*(c(i+1)-c(i));
    end
    a = a(1:n-1)';
    c = c(1:n-1);
end

function y = eval_cubic_spline(a,b,c,d,k,x)
   [~,n] = size(x);
   y = zeros(1,n);
   cur = 2; 
   s = @(x) a(1) + b(1)*(x-k(1)) + c(1) * (x-k(1))^2 + d(1)*(x-k(1))^3;
   for i=1:n
       if x(i) > k(cur)
          s = @(x) a(cur) + b(cur)*(x-k(cur)) + c(cur) * (x-k(cur))^2 + d(cur)*(x-k(cur))^3;
          cur = cur + 1; 
       end
       y(i) = s(x(i));
   end
end

function [a,b,c,d] = ccsi(k, s_0, s_n, h, f)
   % Clamped Cubic Spline Interpolation (ccsi)
   %
   % This MATLAB function calculates the coefficients for clamped cubic spline interpolation.
   %
   % Inputs:
   % 1. k: Vector (1 x n) - A vector of x-coordinates of the data points.
   % 2. s_0: Scalar - The specified first derivative value at the start point.
   % 3. s_n: Scalar - The specified first derivative value at the end point.
   % 4. h: Scalar - The step size, representing the uniform spacing between points in k.
   % 5. f: Function handle - A function to calculate the y-values for the given x-coordinates.
   
   % Outputs:
   % [a, b, c, d] - coefficients of the caluclated cubic spline.
   %
   % The function implements clamped cubic spline interpolation, useful for fitting a smooth curve through a series of data points with specified slope conditions at the endpoints.

   [~, n] = size(k);
   a = arrayfun(f, k);
   A = 4 * h * eye(n) + triu(tril(h*ones(n),1),1) + tril(triu(h*ones(n),-1),-1);
   A(1,1) = A(1,1) - 2*h; 
   A(n,n) = A(n,n) - 2*h;
   sol = zeros(n,1);
   for i=1:n
      if i == 1
         sol(i) = 3 * (1/h *(a(2)-a(1)) - s_0); 
      elseif i == n
          sol(i) = 3 * (s_n - 1/h *(a(i)-a(i-1)));
      else
         sol(i) = 3/h * (a(i+1)-2 * a(i) + a(i-1)); 
      end
   end
   c = A\sol; 
   b = zeros(n-1,1); 
   d = zeros(n-1,1);
   for i=1:n-1
       b(i) = 1/h * (a(i+1)-a(i)) - h/3 * (2 * c(i) + c(i+1)); 
       d(i) = 1/(3*h)*(c(i+1)-c(i));
   end
   a = a(1:n-1)';
   c = c(1:n-1);
end
