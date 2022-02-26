% FRICTION SLOPE FROM GAUCLKER-STRICKLER FORMULA
function Se = friction_slope(Ks, U, D)

% energy loss per unit length
Se = ( U ./ (Ks .* D.^(2/3)) ).^2;

% end of the function
return