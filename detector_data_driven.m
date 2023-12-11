% Author:       Mehran Attar
% Written:      10-December-2023
% Last update:  --------------
% Last revision: 10-December-2023
% This function simulates the data-driven anomaly detector local to the tracking controller, 
% which is in charge of detecting anomalies caused by FDI attacks
% signal 
%------------- BEGIN CODE --------------

function alarm = detector_data_driven(x,x_pre)

if x_pre.contains(x) == 1
    alarm = 0;
else
    alarm = 1;
end

end

