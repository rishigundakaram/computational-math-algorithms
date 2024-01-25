% This MATLAB file is dedicated to exploring numerical solutions of ordinary differential 
% equations (ODEs) using xplicit methods. It includes implementations 
% of third-order, Nystrom, fourth-order ERK schemes, adn the adams bashforth methdo to 
% solve a given ODE

function y_n = thirdOrderERK(f, t, y_0, h)
    % Third-Order Explicit Runge-Kutta (ERK) Scheme
    % Solves an ODE using a third-order ERK method.
    % Inputs:
    % f - Function handle representing the ODE (dy/dt = f(t,y))
    % t - Time interval as a vector [t_start, t_end]
    % y_0 - Initial condition y(t_start)
    % h - Step size for the ERK method
    % Outputs:
    % y_n - Numerical solution of the ODE over the interval
    A = h*[1/2 0; -1 2]; 
    b = [1/6 2/3 1/6]; 
    c = [0 1/2 1]; 
    t_0 = t(1);
    length = t(2) - t(1);
    stages = 3;
    n = length * 1/h + 1;
    y_n = zeros(1, n); 
    y_n(1) = y_0; 
    for k=2:n
        t_cur = t_0 + h * (k-2); 
        k_stages = zeros(stages,1);
        k_stages(1) = f(t_cur, y_n(k-1));
        for j=2:stages
           sum = dot(A(j-1,1:j-1),k_stages(1:j-1));
           k_stages(j) = f(t_cur + h * c(j), y_n(k-1) + sum);
        end
        y_n(k) = y_n(k-1) + h * b * k_stages;
    end
end


function y_n = nystromSchemeERK(f, t, y_0, h)
    % Nystrom ERK Scheme
    % Solves an ODE using the Nystrom ERK method.
    % Inputs and Outputs are similar to thirdOrderERK function.
    A = h*[2/3 0; 0 2/3]; 
    b = [1/4 3/8 3/8]; 
    c = [0 2/3 2/3]; 
    t_0 = t(1);
    length = t(2) - t(1);
    stages = 3;
    n = length * 1/h + 1;
    y_n = zeros(1, n); 
    y_n(1) = y_0; 
    for k=2:n
        t_cur = t_0 + h * (k-2); 
        k_stages = zeros(stages,1);
        k_stages(1) = f(t_cur, y_n(k-1));
        for j=2:stages
           sum = dot(A(j-1,1:j-1),k_stages(1:j-1));
           k_stages(j) = f(t_cur + h * c(j), y_n(k-1) + sum);
        end
        y_n(k) = y_n(k-1) + h * b * k_stages;
    end
end


function y_n = fourthOrderERK(f, t, y_0, h)
    % Fourth-Order Classic ERK Scheme
    % Solves an ODE using a fourth-order classic ERK method.
    % Inputs and Outputs are similar to thirdOrderERK function.
    A = h*[1/2 0 0; 0 1/2 0; 0 0 1]; 
    b = [1/6 1/3 1/3 1/6]; 
    c = [0 1/2 1/2 1]; 
    t_0 = t(1);
    length = t(2) - t(1);
    stages = 4;
    n = length * 1/h + 1;
    y_n = zeros(1, n); 
    y_n(1) = y_0; 
    for k=2:n
        t_cur = t_0 + h * (k-2); 
        k_stages = zeros(stages,1);
        k_stages(1) = f(t_cur, y_n(k-1));
        for j=2:stages
           sum = dot(A(j-1,1:j-1),k_stages(1:j-1));
           k_stages(j) = f(t_cur + h * c(j), y_n(k-1) + sum);
        end
        y_n(k) = y_n(k-1) + h * b * k_stages;
    end
end

function y_n = adamsBashforth(f, t, h, y_n)
    % Adams-Bashforth Method
    % Solves an ODE using the Adams-Bashforth method, specifically a 3-step, 4th-order method.
    % Inputs:
    % f - Function handle representing the ODE (dy/dt = f(t,y))
    % t - Time interval as a vector [t_start, t_end]
    % h - Step size for the method
    % y_n - Initial condition and precomputed values for the first few steps
    % Outputs:
    % y_n - Numerical solution of the ODE over the interval
    % Reference: https://en.wikipedia.org/wiki/Adams%E2%80%93Bashforth_method

    n = 1/h + 1;
    t_0 = t(1);
    t = [1:h:2];
    for k = 4:n
        p1 = f(t(k-1), y_n(k-1));
        p2 = f(t(k-2), y_n(k-2)); 
        p3 = f(t(k-3), y_n(k-3));
        y_n(k) = y_n(k-1) + h * (23/12*p1 - 4/3*p2 + 5/12*p3);
    end
end

