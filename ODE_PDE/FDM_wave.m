%This MATLAB file explores several finite difference schemes for solving the advection 
% equation, a type of partial differential equation (PDE). The file includes 
% implementations of the Forward Time Backward Space (FTBS), Lax-Friedrichs, Lax-Wendroff, 
% and Forward Time Forward Space (FTFS) schemes. It evaluates these methods' accuracy and 
% stability by comparing the numerical solutions with an exact solution and examining the 
% behavior under various grid resolutions.

function u = FTBS(h)
    % Forward Time Backward Space (FTBS) Scheme
    % Solves the advection equation using the FTBS method.
    % Inputs:
    % h - Spatial step size
    % Outputs:
    % u - Solution matrix of the advection equation
    l = .8; 
    k = l*h; 
    x = [-1:h:3];
    t = [0:k:2.4];
    u_0 = (abs(x) <= .5) .* cos(pi * x).^2; 
    x_length = length(x);
    t_length = length(t);
    u = zeros(x_length, t_length); 
    u(:, 1) = u_0;
    t_length = (t(end) - t(1))/k + 1;
    for m=2:t_length+1
        u(2:end, m) = u(2:end, m-1) - l * (u(2:end, m-1) - u(1:end-1, m-1));
    end
end

function u = laxFriedrichs(h)
    % Lax-Friedrichs Scheme
    % Solves the advection equation using the Lax-Friedrichs method.
    % Inputs and Outputs are similar to the FTBS function.
    % [Lax-Friedrichs implementation...]
    l = .8;
    k = l*h; 
    x = [-1:h:3];
    t = [0:k:2.4];
    u_0 = (abs(x) <= .5) .* cos(pi * x).^2; 
    x_length = length(x);
    t_length = length(t);
    u = zeros(x_length, t_length); 
    u(:, 1) = u_0;
    t_length = (t(end) - t(1))/k + 1;
    for m=2:t_length + 1
        u(end, m) = u(end, m-1);
        u(2:end-1, m) = 1/2*(u(1:end-2, m-1) + u(3:end, m-1)) -1/2*k/h * (u(3:end, m-1)-u(1:end-2, m-1));
    end
end

function u = laxWendroff(h)
    % Lax-Wendroff Scheme
    % Solves the advection equation using the Lax-Wendroff method.
    % Inputs and Outputs are similar to the FTBS function.
    l = .8;
    k = l*h; 
    x = [-1:h:3];
    t = [0:k:2.4];
    u_0 = (abs(x) <= .5) .* cos(pi * x).^2; 
    x_length = length(x);
    t_length = length(t);
    u = zeros(x_length, t_length); 
    u(:, 1) = u_0;
    t_length = (t(end) - t(1))/k + 1;
    for m=2:t_length + 1
        u(end, m) = u(end, m-1);
        u(2:end-1, m) = u(2:end-1,m-1)-1/2*(k/h)*(u(3:end, m-1)-u(1:end-2, m-1))+1/2*(k/h)^2*(u(3:end, m-1)-2*u(2:end-1, m-1)+u(1:end-2, m-1));
    end
end

function u = FTFS(h)
    % Forward Time Forward Space (FTFS) Scheme
    % Solves the advection equation using the FTFS method.
    % Inputs and Outputs are similar to the FTBS function.
    l = .8;
    k = l*h; 
    x = [-1:h:3];
    t = [0:k:2.4];
    u_0 = (abs(x) <= .5) .* cos(pi * x).^2; 
    x_length = length(x);
    t_length = length(t);
    u = zeros(x_length, t_length); 
    u(:, 1) = u_0;
    t_length = (t(end) - t(1))/k + 1;
    for m=2:t_length + 1
        u(2:end-1, m) = u(2:end-1, m-1) - l * (u(3:end, m-1)- u(2:end-1, m-1));
    end
end

