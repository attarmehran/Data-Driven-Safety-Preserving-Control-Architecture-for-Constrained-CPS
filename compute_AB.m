
% Author:       Mehran Attar
% Written:      10-December-2023
% Last update:  --------------
% Last revision: 10-December-2023 
%---------------------------------------------------------------
% This function computes the set of all system matrices, \mathcal{M}_{AB} that is consistent
% with the data. Moreover, this function computes the set of vertices
% \mathcal{V}_{AB}
%------------- BEGIN CODE --------------

function [V_AB,AB,X_0T,X_1T,u] = compute_AB(sys,X0,U,W)

w = warning ('off','all');
rmpath('folderthatisnotonpath')
warning(w)

rand('seed',1);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

dim_x = size(A,1);
dim_u = size(B,2);

%Number of trajectories
initpoints =1;
%Number of time steps
steps = 1;
initpoints = 3;
%Number of time steps
steps = 2;
totalsamples = initpoints*steps;
%% initial set and input

%Construct matrix zonotpe \mathcal{M}_w
index=1;
for i=1:size(W.generators,2)
    vec = W.Z(:,i+1);
    GW{index}= [vec,zeros(dim_x,totalsamples-1)];
    for j=1:totalsamples-1
        GW{j+index}=  [GW{index+j-1}(:,2:end) GW{index+j-1}(:,1)];
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

X1W_cen =  (X_1T - Wmatzono.center) * pinv([X_0T;u]);
for i=1:size(Wmatzono.generator,2)
   Wmatzonos{i} = Wmatzono.generator{i}* pinv([X_0T;u]);
end

AB = matZonotope(X1W_cen,Wmatzonos);

matrixCenter = [AB.center;zeros(size(B,2),size(A,2)+size(B,2))];
% 
for i=1:AB.gens
  G{i}=[AB.generator{i};zeros(size(B,2),size(A,2)+size(B,2))];  
end

% instantiate matrix zonotope
M_zono = matZonotope(matrixCenter, G);

% obtain result of all vertices 
V_AB = vertices(M_zono);

disp('The rank of data is: ' + string(rank([u;X_0T])))

end

