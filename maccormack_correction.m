% MACCORMACK SCHEME - CORRECTION STEP
function [Dc, Uc] = maccormack_correction(g, dt, dx, S0, n, ...
                                    Ut, Dt, Uf, Df, Se, Dus, Uus, Dds, Uds)

% allocate arrays
Dc = nan(n,1);
Uc = nan(n,1);

% set boundary conditions
Dc(1)   = Dus;
Uc(1)   = Uus;
Dc(end) = Dds;
Uc(end) = Uds;

% backward loop over nodes
for j = n-1:-1:2
    Dc(j) = Dt(j) - dt/dx * ( Uf(j+1)*Df(j+1) - Uf(j)*Df(j) );
    Uc(j) = Ut(j) - dt/dx * ( (Uf(j+1)^2 - Uf(j)^2)/2 + ...
            g*(Df(j+1)-Df(j))) - g * dt * (Se(j)-S0(j));
end

% end of the function
return