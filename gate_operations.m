% GATE OPENING
function [gate_a, jgate] = gate_operations(t, gate_times, gate_opening, jgate)


if t >= gate_times(jgate)
    jgate = min(length(gate_times), jgate + 1);
end
if t >= gate_times(end)
    gate_a = gate_opening(end);
else
    gate_a = gate_opening(jgate-1) + ...
            ( gate_opening(jgate) - gate_opening(jgate-1) ) / ...
            ( gate_times(jgate) - gate_times(jgate-1) ) * ...
            max(0, t - gate_times(jgate-1));
end

% end of the function
return