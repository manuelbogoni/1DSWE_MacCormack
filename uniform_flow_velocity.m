% UNIFORM FLOW
function U = uniform_flow_velocity(D, Ks, S0)

% unit discharge
U = Ks * D^(2/3) * sqrt(S0);

% end of the function
return
