% Author:       Mehran Attar
% Written:      10-December-2023
% Last update:  --------------
% Last revision: 10-December-2023
% This function computes a zonotopic inner approximation of a polytope 
%------------- BEGIN CODE --------------

function [out_approx,alpha_out] = poly_approx(P,num_gen,generator_matrix)

dim = P.Dim;
alpha = sdpvar(num_gen, 1);
center = sdpvar(dim, 1);
weight = zeros(num_gen, 1);

for i = 1:num_gen
    weight(i, :) = norm(generator_matrix(:, i), 2);
end

constraints = [P.A * center + abs(P.A * generator_matrix) * alpha <= P.b];

opt = sdpsettings('verbose',0);
diagnostics = optimize(constraints,sum(-weight .* log(alpha)),opt);
center = value(center);
alpha = value(alpha);
alpha_out = norm(value(alpha));

for i = 1:num_gen
    generator(:, i) = generator_matrix(:, i) * value(alpha(i));
end

out_approx = zonotope(center,generator);
end

%------------- END CODE --------------