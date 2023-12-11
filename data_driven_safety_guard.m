% Author:       Mehran Attar
% Written:      10-December-2023
% Last update:  --------------
% Last revision: 10-December-2023
% This function checks the safety of the plant using the received control
% signal 
%------------- BEGIN CODE --------------
function [x_1,is_safe] = data_driven_safety_guard(u,x,U,X,AB,W)

x = zonotope(x,0*diag(ones(2,1)));
u = zonotope(u,0*diag(ones(3,1)));

x_1 = AB * (cartProd(x,u))+ W;  % one-step evolution of the system, \hat{\mathcal{S}}

if U.contains(u) && X.contains(x_1) == 1
    is_safe = 0;
else
    is_safe = 1;
end 

end
%------------- END OF CODE --------------

