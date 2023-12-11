function alarm = detector_data_driven(x,x_pre)

% x = zonotope(x,0*diag(ones(2,1)));
% u = zonotope(u,0*diag(ones(3,1)));

% x_1 = AB * (cartProd(x,u))+ 1*W;
if x_pre.contains(x) == 1
    alarm = 0;
else
    alarm = 1;
end

end

