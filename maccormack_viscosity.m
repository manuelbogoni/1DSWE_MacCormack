% MACCORMACK SCHEME - ARTIFICIAL VISCOSITY
function [Dv, Uv] = maccormack_viscosity(n, D, U)

% allocate arrays
Dv = D;
Uv = U;

% weight coefficient
alpha = 0.5;

% viscosity coefficient
nu = nan(n,1);
for j = 2:n-1
    nu(j) = abs(D(j+1) - 2* D(j) + D(j-1)) / ...
            ( abs(D(j+1)) + 2 * abs(D(j)) + abs(D(j-1)) );
end

% first point
nu(1) = abs(D(2) - D(1)) / ( abs(D(2)) + abs(D(1)) );

% last point
nu(end) = abs(D(end) - D(end-1)) / ( abs(D(end)) + abs(D(end-1)) );

% redistribution coefficient
beta = nan(length(D)-1,1);
for j = 1:n-1
    beta(j) = alpha * max(nu(j+1), nu(j));
end

% update variables
for j = 2:n-1
    Dv(j) = D(j) + beta(j) * (D(j+1) - D(j)) - ...
                   beta(j-1) * (D(j) - D(j-1));
    Uv(j) = U(j) + beta(j) * (U(j+1) - U(j)) - ...
                   beta(j-1) * (U(j) - U(j-1));
end

% end of the function
return