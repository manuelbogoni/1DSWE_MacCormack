% COURANT - FRIEDRICHS - LEWY CONDITION
function dt = cfl(g, dx, U, D)

% compute time step
dt = 0.95 * (dx / max( U + sqrt(g*D)));

% end of the function
return