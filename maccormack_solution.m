%MACCORMACK SCHEME - SOLUTION STEP
function [D, U] = maccormack_solution(Dp, Up, Dc, Uc)

D = 0.5 * (Dp + Dc);
U = 0.5 * (Up + Uc);

%end of the function
return