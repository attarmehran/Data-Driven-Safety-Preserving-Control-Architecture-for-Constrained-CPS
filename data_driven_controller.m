% Author:       Mehran Attar
% Written:      10-December-2023
% Last update:  --------------
% Last revision: 10-December-2023
% This function computes a data-driven tracking controller 
%------------- BEGIN CODE --------------

function K = data_driven_controller(A,B,W,X0,U)

dim_x = size(A,1);
dim_u = size(B,2);

%Number of trajectories
initpoints =1;
%Number of time steps
steps = 1;
initpoints =1;
%Number of time steps
steps = 5;
totalsamples = initpoints*steps;
%% initial set and input

%Construct matrix zonotpe \mathcal{M}_w
index=1;
for i=1:size(W.generators,2)
    vec=W.Z(:,i+1);
    GW{index}= [ vec,zeros(dim_x,totalsamples-1)];
    for j=1:totalsamples-1
        GW{j+index}= [GW{index+j-1}(:,2:end) GW{index+j-1}(:,1)];
    end
    index = j+index+1;
end
Wmatzono= matZonotope(zeros(dim_x,totalsamples),GW);


% randomly choose constant inputs for each step / sampling time
for i=1:totalsamples
    u(:,i) = randPoint(U);
end

%simulate the system to get the data
x0 = X0.center;
x(:,1) = x0;
index=1;
for j=1:dim_x:initpoints*dim_x
    x(j:j+dim_x-1,1) = randPoint(X0);
    for i=1:steps
        utraj(j,i) = u(index);
        noise(:,index) = randPoint(W);
        x(j:j+dim_x-1,i+1) = A*x(j:j+dim_x-1,i) + B*u(:,index) + randPoint(W);
        index=index+1;
    end
end


% concatenate the data trajectories
index_0 =1;
index_1 =1;
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1
        x_meas_vec_1(:,index_1) = x(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        u_mean_vec_0(:,index_0) = utraj(j,i);
        x_meas_vec_0(:,index_0) = x(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end

% X_+ is X_1T
% X_- is X_0T
U_full = u_mean_vec_0(:,1:totalsamples); %same as u
X_0T = x_meas_vec_0(:,1:totalsamples);
X_1T = x_meas_vec_1(:,1:totalsamples);
rank([u;X_0T])

% define the LMI for computing the gain of the controller 
setlmis([])

[Q,nQ,sQ] = lmivar(2,[steps dim_x]);

lmiterm([-1 1 1 Q],X_0T,1);
lmiterm([-1 1 2 Q],X_1T,1);
lmiterm([-1 2 2 Q],X_0T,1);

lmisys = getlmis;
ndec = decnbr(lmisys);
nlmi = lminbr(lmisys);

[t_min,q_feas] = feasp(lmisys);
Qval = dec2mat(lmisys,q_feas,Q);
K = u * Qval * inv(X_0T*Qval);
end
%------------- END OF CODE --------------
