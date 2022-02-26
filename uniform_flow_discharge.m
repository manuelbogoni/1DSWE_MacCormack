% UNIFORM FLOW
function q = uniform_flow_discharge(D, Ks, S0)

% unit discharge
q = Ks .* D.^(5/3) * sqrt(S0);

% end of the function
return
