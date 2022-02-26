% FLOW BENEATH A SLUICE GATE
function q = sluice_gate(g, Cc, a, H0)

% discharge coefficient of a sluice gate
Cq  = Cc * sqrt(1 / (1 + Cc * a/ H0));

% unit discharge of a flow beneath a sluice gate
q   = Cq * a * sqrt(2 * g * H0);

% end of the function
return
