% UNIFORM FLOW
function D = uniform_flow_depth(q, Ks, S0, flag)

% normal depth from velocity
if flag == 'q'
    D = ( q ./ (Ks .* sqrt(S0)) ) .^ (3/5);

% normal depth from discharge
elseif flag == 'v' || flag == 'U'
    U = q;
    D = ( U ./ (Ks .* sqrt(S0)) ) .^ (3/2);
end

% end of the function
return
