% MACCORMACK SCHEME - PREDICTION STEP
function [Dp, Up] = maccormack_prediction(g, dt, dx, S0, n, ...
                                    Ut, Dt, Uf, Df, Se, Dus, Uus, Dds, Uds)

% allocate arrays
Dp = nan(n,1);
Up = nan(n,1);

% set boundary conditions
Dp(1)   = Dus;
Up(1)   = Uus;
Dp(end) = Dds;
Up(end) = Uds;

% forward loop over nodes
for j = 2:n-1
    Dp(j) = Dt(j) - dt/dx * ( Uf(j)*Df(j) - Uf(j-1)*Df(j-1) );
    Up(j) = Ut(j) - dt/dx * ( (Uf(j)^2 - Uf(j-1)^2)/2 + ...
            g*(Df(j)-Df(j-1))) - g * dt * (Se(j)-S0(j));
end

% end of the function
return