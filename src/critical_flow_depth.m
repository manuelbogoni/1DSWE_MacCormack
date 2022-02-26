% CRITICAL FLOW
function D = critical_flow_depth(g, q)

% critical depth from discharge
D = ( q.^2 / g).^(1/3);

% end of the function
return
