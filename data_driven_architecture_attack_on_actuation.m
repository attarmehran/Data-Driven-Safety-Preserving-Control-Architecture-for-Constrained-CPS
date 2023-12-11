
% Author:       Mehran Attar
% Written:      10-December-2023
% Last update:  --------------
% Last revision: 10-December-2023
% This codes simulates the proposed data-driven safety preserving control
% architecture in the presence of cyber attacks on the actuation channel. 

%------------- BEGIN CODE --------------

clc
clear all
close all

rand('seed',11);
w = warning ('off','all');
rmpath('folderthatisnotonpath')
warning(w)

% loading ROSC sets and checking the existence of ROSC sets
if isfile('Td.mat') && isfile('Td_aug.mat')
    disp('ROSC sets are available!')
else
    disp('ROSC sets are NOT available!')
    run('computing_ROSC_sets.m')
end

Td = load('Td.mat').Td;   % ROSC sets
V_AB = load('V_AB.mat').V_AB;  % set of vertices \mathcal{V}_{AB} as a matrix zonotope
AB = load('AB.mat').AB;   % set of system matrices \mathcal{M}_{AB} as a matrix zonotope
Td_aug = load('Td_aug.mat').Td_aug;  % set of augmented ROSC sets

% simulation settings
sim_time = 200;    % simulation time

% defining system matrices
A=[0.993 0.003;0.007 0.982];
B=[0.008 -0.003 -0.003;0 0.003 0.003];
C=eye(2);
D=0;

dim_x = size(A,1);
dim_u = size(B,2);

% defining system constraints
X = zonotope(interval([-0.48;-0.48],[0.3;0.3]));  % constraints on states
U = zonotope(interval([-0.7778;-1.25;-1.4765],[0.611;0.75;0.5235])); % constraints on inputs
W = zonotope(zeros(dim_x,1),0.001*eye(dim_x));  % noise zonotope

sys = ss(A,B,C,D);  % defining system

% computing system matrices using set of vertices and data samples
for i=1:size(V_AB,1)
    A_hat{i} = V_AB{i}(1:size(A,1),1:size(A,1));
    B_hat{i} = V_AB{i}(1:size(A,1),size(A,1)+1:size(A,1)+size(B,2));
end

% defining reference signal
for i=1:sim_time
    if i<100
        ref(:,i)=[0.1;0.03];
    elseif i>=100 && i<=200
        ref(:,i)=[0.1;-0.1];
    end
end

% defining attack on the actuation channel
for i=1:sim_time
    if i>=95 && i<=113
        u_a(:,i)=[1;1;2];
    else
        u_a(:,i)=[0;0;0];
    end
end

% computing data-driven tracking controller
K_data = data_driven_controller(A,B,W,X,U);

% initialization
alarm_data(1) = 0;
x_data(:,1) = [0.01;-0.01];  % initial state
flag = 0;
ignore = 0;
emergency(1) = 0;

% computing equlibirium points
for i=1:sim_time
    u_eq(:,i) = pinv(B)*((eye(dim_x)-A)*ref(:,i));
end

% Visualization
close all
f = figure;
ax = axes;
f.Position = [700 70 800 700]

plot(x_data(1,1),x_data(2,1),'*','MarkerSize',4,'MarkerEdgeColor','k')
hold on
annotation('textarrow',[0.84525 0.819015325670498],...
    [0.845428571428572 0.78886229086229],'String','$\hat{\mathcal{T}}^{40}_e$',...
    'Interpreter','latex',...
    'FontSize',14);
hold on
plot(Td{1}.mptPolytope.P,'Alpha',0.7,...
    'color','green','LineWidth',0.01)
hold on
for i=2:40
    plot(Td{i},[1 2],'color',[0.5020    0.5020    0.5020],'LineStyle','-','LineWidth',1)
    hold on
end
xlim([-0.1 0.12])
ylim([-0.15 0.12])
xlabel('$\tilde{x}_1$','interpreter','latex','FontSize',24)
ylabel('$\tilde{x}_2$','interpreter','latex','FontSize',24)
hold on
text(-0.08,0.025,'$\hat\mathcal{T}^0_e$','FontSize',17,'interpreter','latex')
hold on
text(0.103,-0.1,'$r_{100}$','FontSize',17,'interpreter','latex')
hold on
text(0.103,0.03,'$r_1$','FontSize',17,'interpreter','latex')
hold on
text(0.01,0,'$x_0$','FontSize',17,'interpreter','latex')
hold on
handle_ref1 = plot(ref(1,1),ref(2,1),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',0.8)
hold on
handle_ref2 = plot(ref(1,120),ref(2,120),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',1)
hold on
box on
grid off

% simulation of the system -- presence of attack on the actuation channel
for k=1:sim_time
    
    % computing data-driven tracking controller
    ctr_data(:,k) = K_data*(x_data(:,k)-ref(:,k)) + u_eq(:,k);
    ctr_data(:,k) =  min(max(ctr_data(:,k), [-0.7778;-1.25;-1.4765]), [0.611;0.75;0.5235]);
    
    % computing the one-step evolution set, \hat{\mathcal{R}}_k
    x = zonotope(x_data(:,k),0*diag(ones(dim_x,1)));
    u = zonotope(ctr_data(:,k),0*diag(ones(dim_u,1)));
    x_pre_data{k} = AB * (cartProd(x,u))+ W; % \hat{\mathcal{R}}_k as a zonotope
    
    % defining attack on the actuation channel - intelligent attacker
    ctr_data_prime(:,k) = ctr_data(:,k) + u_a(:,k);
    ctr_data_prime(:,k) =  min(max(ctr_data_prime(:,k),...
        [-0.7778;-1.25;-1.4765]), [0.611;0.75;0.5235]);
    
    % attacker changes the flag signal
    if alarm_data(k)==1 && emergency(k-1) == 0
        flag = 0;
    end
    
    %
    if flag == 1 && Td{1}.contains(x_data(:,k)) == 1
        flag = 0;
        ignore = 1;
    else
        ignore = 0;
    end
    %
    
    % safety check
    [x_plus_data{k},safety_data(k)] = data_driven_safety_guard(ctr_data_prime(:,k),...
        x_data(:,k),U,Td{40},AB,W);
    
    
    if safety_data(k) == 1
        flag = 1;
    end
    
    % Swithcing between emergency controller and tracking controller
    if flag == 1
        index_data(k) = set_index(x_data(:,k),Td);
        u_ver(:,k) = one_step_ctrl(3, x_data(:,k), Td_aug, index_data(k));
        emergency(k) = 1;
    else
        u_ver(:,k) = ctr_data_prime(:,k);
        emergency(k) = 0;
    end
    
    % computing the system evolution
    x_data(:,k+1) = A*x_data(:,k) + B*u_ver(:,k) + randPoint(W);
    
    % visualization of system evolution
    plot(x_data(1,k+1),x_data(2,k+1),'o','MarkerSize',2,'MarkerEdgeColor','r','MarkerFaceColor','r');
    
    
    % anomaly detector
    alarm_data(k+1) = detector_data_driven(x_data(:,k+1),x_pre_data{k});
    
    % send a flag=1 in the case of presence of attacks 
    if alarm_data(k+1)==1
        flag = 1;
    end
    
    %
    if norm(u_a(:,k))==[0;0;0]
        attack(k)=false;
    else
        attack(k)=true;
    end
    
    title('Time: '+string(k)+', flag='+string(flag)+', emergency='+string(emergency(k))+...
        ', ignore='+string(ignore)+', D_k='+string(alarm_data(k))+', attack on u_k='+ string(attack(k)),...
        'FontSize',11,...
        'FontName','Courier')
    pause(0.1)
    k
end

alarm = alarm_data(2:201);
active = emergency(2:200);
safety = safety_data(2:200);

%% Visualization of the system evolution 
f = figure;
ax = axes;
f.OuterPosition = [200 50 850 800]
ax.OuterPosition = [-0.06     -0.04     1.15     1.11]

plot(Td{1}.mptPolytope.P,'Alpha',0.7,...
    'color','green','LineWidth',0.01)
hold on
for i=2:40
    plot(Td{i},[1 2],'color',[0.5020    0.5020    0.5020],'LineStyle','-','LineWidth',0.1)
    hold on
end

hold on
x1 = [0.877394636015326 0.896551724137931];
y1 = [0.420045045045045 0.676801801801802];
x2 = [0.699233716475094 0.89272030651341];
y2 = [0.81644144144144 0.757882882882883];

annotation('textarrow',x1,y1,'String','Attack Starting Point','FontSize',12)
hold on
annotation('textarrow',x2,y2,'String','Attack Detection Point','FontSize',12)
hold on
text(-0.08,0.025,'$\hat\mathcal{T}^0_e$','FontSize',17,'interpreter','latex')
hold on
text(0.103,-0.1,'$r_{100}$','FontSize',17,'interpreter','latex')
hold on
text(0.103,0.03,'$r_1$','FontSize',17,'interpreter','latex')
hold on
text(0.098,0.055,'$x_{104}$','FontSize',17,'interpreter','latex')
hold on
text(0.01,0,'$x_0$','FontSize',17,'interpreter','latex')
hold on
annotation('textarrow',[0.921455938697315 0.886015325670498],...
    [0.879630630630629 0.826576576576576],'String','$\hat{\mathcal{T}}^{40}_e$',...
    'Interpreter','latex',...
    'FontSize',14);
hold on
handle_tracking = plot(x_data(1,1:50),x_data(2,1:50),'b-','LineWidth',1.5)
hold on
handle_attack_tracking = plot(x_data(1,95:105),x_data(2,95:105),'r-','LineWidth',1.5)
hold on
handle_em_control = plot(x_data(1,105:115),x_data(2,105:115),'color'...
    ,[0.1294    0.6588    0.0588],'LineWidth',...
    2,'LineStyle','-')
hold on
handle_tracking_after = plot(x_data(1,115:180),x_data(2,115:180),'m-','LineWidth',1.5)
hold on

handle_initial_state = plot(x_data(1,1),x_data(2,1),'O','MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',0.8)
hold on
handle_ref1 = plot(ref(1,1),ref(2,1),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',0.8)
hold on
handle_ref2 = plot(ref(1,120),ref(2,120),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',1)
hold on
plot(x_data(1,105),x_data(2,105),'O','MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.1294    0.6588    0.0588],'LineWidth',0.8)
hold on
plot(Td{40},[1 2], 'k-', 'LineWidth', 1)
xlim([-0.1 0.12])
ylim([-0.15 0.12])
xlabel('$\tilde{x}_1$','interpreter','latex','FontSize',24)
ylabel('$\tilde{x}_2$','interpreter','latex','FontSize',24)
legend([handle_tracking,handle_attack_tracking,handle_em_control,handle_tracking_after],...
    'Tracking Controller Before Attack',...
    'Tracking Controller + Attack',...
    'Emergency Controller',...
    'Tracking Controller After Attack',...
    'Location','NorthWest','FontSize',14)
grid off
box on


ax=axes;
set(ax,'units','normalized','position',[0.15,0.15,0.3,0.3])
box(ax,'on')
hold on

for i=40:40
    plot(Td{i},[1 2],'color','black','LineStyle','-','LineWidth',0.1)
    hold on
end
hold on
handle_tracking = plot(x_data(1,1:27),x_data(2,1:27),'b-','LineWidth',1.5)
hold on
handle_attack_tracking = plot(x_data(1,94:105),x_data(2,94:105),'r-','LineWidth',1.5)
hold on
handle_em_control = plot(x_data(1,105:116),x_data(2,105:116),'color'...
    ,[0.1294    0.6588    0.0588],'LineWidth',...
    2,'LineStyle','-')
hold on
handle_tracking_after = plot(x_data(1,116:180),x_data(2,116:180),'m-','LineWidth',1.5)
hold on

handle_initial_state = plot(x_data(1,1),x_data(2,1),'O','MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r','LineWidth',0.8)
hold on
handle_ref = plot(ref(1,95),ref(2,95),'O','MarkerSize',5,'MarkerEdgeColor','r',...
    'MarkerFaceColor','r','LineWidth',0.8)
hold on
handle_ref1 = plot(ref(1,50),ref(2,50),'O','MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r','LineWidth',1)
hold on
plot(x_data(1,105),x_data(2,105),'O','MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.1294    0.6588    0.0588],'LineWidth',0.8)
hold on
plot(x_plus_data{104},[1 2],'b--','LineWidth',1)
hold on
plot(x_plus_data{105},[1 2],'r--','LineWidth',1)
hold on
plot(x_plus_data{106},[1 2],'Color',[0.1294    0.6588    0.0588],'LineWidth',1,'LineStyle','--')
hold on
annotation('textarrow',[0.276819923371647 0.307471264367816],...
    [0.288414414414414 0.32545045045045],'String','$\hat{\mathcal{S}}^{+}_{103}$','FontSize',...
    14,'interpreter','latex');
hold on
annotation('textarrow',[0.418582375478927 0.404214559386971],...
    [0.287162162162162 0.332207207207207],'String','$\hat{\mathcal{S}}^{+}_{104}$','FontSize',...
    14,'interpreter','latex');
hold on
annotation('textarrow',[0.231800766283524 0.266283524904214],...
    [0.316441441441441 0.336711711711711],'String','$\hat{\mathcal{S}}^{+}_{105}$','FontSize',...
    14,'interpreter','latex');
hold on
annotation('textarrow',[0.401580459770114 0.389846743295019],...
    [0.39724099099099 0.372747747747748],'String','$\hat{\mathcal{T}}^{40}_e$',...
    'Interpreter','latex',...
    'FontSize',14);
hold on
annotation('textbox',...
    [0.228250957854406 0.405405405405405 0.173329501915709 0.0360360360360361],...
    'String','Safety Verification',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
hold on
annotation('line',[0.450191570881226 0.813218390804598],...
    [0.447198198198198 0.800675675675676],...
    'Color',[0.411764705882353 0.411764705882353 0.411764705882353],...
    'LineStyle','--');
hold on
annotation('line',[0.448275862068965 0.814176245210728],...
    [0.150900900900901 0.632882882882883],...
    'Color',[0.411764705882353 0.411764705882353 0.411764705882353],...
    'LineStyle','--');
hold on
annotation('rectangle',...
    [0.811302681992337 0.628378378378378 0.154214559386973 0.173423423423424],...
    'LineStyle','--');

xlim([0.08 0.11])
ylim([-0.01 0.08])
box on

print -depsc -tiff -r300 -painters case_A_state_trajectory.eps
%% Visualization of the alarm signal & safety check signal and emergency signal  

k = 1:201;
f = figure;
f.OuterPosition = [400 400 550 250]
ax = axes;
ax.OuterPosition = [-0.04     -0.03     1.08     1.1]
%%%%%%%%%%%%%%%%%%%%
x = [95 113 113 95];
y = [-0.9 -0.9 1.1 1.1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
%
hold on
colororder({'k','k'})
yyaxis left
xlim([0 200])
handle_alarm = plot(k(80:150),alarm(80:150),'r-','LineWidth',2);
ylim([-0.1 1.1])
xlim([0 200])
names = {'Normal'; 'Anomaly'};
set(gca,'ytick',[0:1],'yticklabel',names)
xticks([])
hold on
yyaxis right
handle_emergency = plot(k(80:150),active(80:150),'b--','LineWidth',2);
hold on
handle_safety = plot(k(80:150),safety(80:150),'m--','LineWidth',1);
ylim([-0.1 1.1])
names = {'False'; 'True'};
set(gca,'ytick',[0:1],'yticklabel',names)
xlabel('$k$','interpreter','latex','FontSize',15)

xticks([80 90 96 104 115 120 130 140 150])
legend([handle_alarm,handle_emergency,handle_safety,handle_attack],...
    '$D_k$','emergency','$\hat{\mathcal{S}}^{+}_k\not \subseteq \mathcal{X}_{\eta} \lor u^{\prime}_k \notin \mathcal{U}$',...
    'attack',...
    'interpreter','latex','Location','NorthEast',...
    'FontSize',11)
xlim([80 150])
box on

print -depsc -tiff -r400 -painters case_A_alarm_safety_revised.eps

%------------- END OF CODE --------------