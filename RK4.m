function [y] = RK4(odefun, tspan, y0)
% ODEFUN contains the ode functions of the system
% TSPAN  is a 1D vector of equally spaced t values
% Y0     contains the intial conditions for the system variables

    % Initialise step-size variables
    t = tspan(:); % ensure column vector = (0:h:1)';
    h = t(2)-t(1);% define h from t
    N = length(t);

    % Initialise y vector, with a column for each equation in odefun
    y = zeros(N, numel(y0));

    % Starting conditions
    y(1, :) = y0(:)';  % Set intial conditions using row vector of y0

    k = zeros(4, numel(y0));              % Initialise K vectors
    b = repmat([1 2 2 1]', 1, numel(y0)); % RK4 coefficients

    % Iterate, computing each K value in turn, then the i+1 step values
    for i = 1:(N-1)
        
        k(1, :) = odefun(t(i), y(i,:));        
        k(2, :) = odefun(t(i) + (h/2), y(i,:) + (h/2)*k(1,:));        
        k(3, :) = odefun(t(i) + (h/2), y(i,:) + (h/2)*k(2,:));        
        k(4, :) = odefun(t(i) + h, y(i,:) + h*k(3,:));

        y(i+1, :) = y(i, :) + (h/6)*sum(b.*k);    
        
        
    end    
end