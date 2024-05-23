
% Author:       Mehran Attar
% Written:      12-Jan-2024
% Last update:  23-May-2024
% Last revision: 23-May-2024
% This codes simulates a data-driven safety preserving control
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
sim_time = 201;    % simulation time

% defining system matrices
A = [0.993 0.003;0.007 0.982];
B = [0.008 -0.003 -0.003;0 0.003 0.003];
C = eye(2);
D = 0;

% computing the dimesions of the system's matrices 
dim_x = size(A,1);
dim_u = size(B,2);

% defining system constraints
X = zonotope(interval([-0.48;-0.48],[0.3;0.3]));  % constraints on states
U = zonotope(interval([-0.7778;-1.25;-1.4765],[0.611;0.75;0.5235])); % constraints on inputs
W = zonotope(zeros(dim_x,1),0.001*eye(dim_x));  % noise zonotope

sys = ss(A,B,C,D);  % defining system in the state space

% computing system matrices using set of vertices and data samples
for i=1:size(V_AB,1)
    A_hat{i} = V_AB{i}(1:size(A,1),1:size(A,1));
    B_hat{i} = V_AB{i}(1:size(A,1),size(A,1)+1:size(A,1)+size(B,2));
end

% defining reference signal
for i=1:sim_time
    if i<100
        ref(:,i)=[0.1;0.03];
    elseif i>=100 && i<=sim_time
        ref(:,i)=[0.1;-0.1];
    end
end

% defining attack on the actuation channel
for i=1:sim_time
    if i>=40 && i<=60
        u_a(:,i)=-[1;2;2];
    elseif i>=145 && i<155
        u_a(:,i)=[2;-2;1];
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

% Visualization of the state space 
close all
f = figure;
ax = axes;
f.Position = [700 70 800 700];

plot(x_data(1,1),x_data(2,1),'*','MarkerSize',4,'MarkerEdgeColor','k');
hold on
text(0.1,0.07,'$\hat{\mathcal{T}}^{40}_e$','FontSize',17,'interpreter','latex');
hold on
plot(Td{1}.mptPolytope.P,'Alpha',0.7,...
    'color','green','LineWidth',0.01);
hold on
for i=2:40
    plot(Td{i},[1 2],'color',[0.5020    0.5020    0.5020],'LineStyle','-','LineWidth',1); % family of ROSC sets
    hold on
end
xlim([-0.1 0.16]);
ylim([-0.15 0.12]);
xlabel('$\tilde{x}_1$','interpreter','latex','FontSize',24);
ylabel('$\tilde{x}_2$','interpreter','latex','FontSize',24);
hold on
text(-0.08,0.025,'$\hat\mathcal{T}^0_e$','FontSize',17,'interpreter','latex');
hold on
text(0.103,-0.1,'$r_{100}$','FontSize',17,'interpreter','latex');
hold on
text(0.103,0.03,'$r_1$','FontSize',17,'interpreter','latex');
hold on
text(0.01,0,'$x_0$','FontSize',17,'interpreter','latex');
hold on
handle_ref1 = plot(ref(1,1),ref(2,1),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',0.8);
hold on
handle_ref2 = plot(ref(1,120),ref(2,120),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',1);
hold on
box on
grid off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% simulation of the system -- presence of attack on the actuation channel %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    if (k>=40 && k<=60) | (k>=145&k<155)
        ctr_data_prime(:,k) =  min(max(ctr_data_prime(:,k),...
        [-0.7778;-1.25;-1.4765]), [0.611;0.75;0.5235]);
    end
    
    % attacker changes the flag signal
    if alarm_data(k)==1 && emergency(k-1) == 0
        flag = 0;
    end
    
    % ignore the flag when the trajectory reach the terminal region 
    if flag == 1 && Td{1}.contains(x_data(:,k)) == 1
        flag = 0;
        ignore = 1;
    else
        ignore = 0;
    end
    %
    
    % safety verification 
    [x_plus_data{k},safety_data(k)] = data_driven_safety_guard(ctr_data_prime(:,k),...
        x_data(:,k),U,Td{40},AB,W);
    
    
    if safety_data(k) == 1
        flag = 1;
    end
    
    % Swithcing between emergency controller and tracking controller
    if flag == 1
        index_data(k) = set_index(x_data(:,k),Td);
        u_ver(:,k) = one_step_ctrl(3, x_data(:,k), Td_aug, index_data(k));
        u_ver(:,k) = min(max(u_ver(:,k),...
        [-0.7778;-1.25;-1.4765]), [0.611;0.75;0.5235]);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% simulation of the system without the safety modules -- presence of attack on the actuation channel %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on 
x_w(:,1) = [0.01;-0.01];  % initial state 

for k=1:sim_time
    u_w(:,k) = K_data*(x_w(:,k)-ref(:,k)) + u_eq(:,k);  % tracking controller 
    u_w(:,k) =  min(max(u_w(:,k), [-0.7778;-1.25;-1.4765]), [0.611;0.75;0.5235]); % saturation module
    u_w_prime(:,k) = u_w(:,k) + u_a(:,k); % attack on the actuation channel 
    
    % system evolution 
    x_w(:,k+1) = A*x_w(:,k) + B*u_w_prime(:,k) + randPoint(W);
    
    if Td{40}.contains(x_w(:,k)) == 1
        safety_vio(k) = 0;
    else
        safety_vio(k) = 1;
    end
    
    % visualization of the system evolution without safety modules 
    plot(x_w(1,k+1),x_w(2,k+1),'o','MarkerSize',2,'MarkerEdgeColor','b','MarkerFaceColor','b');
    pause(0.1)
    hold on 
    k
end
%% Visualization of the system evolution 
% close all

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
annotation('textarrow',[0.885057471264368 0.935823754789272],...
    [0.135135135135135 0.337837837837838],'String','Attack Detection Point');
hold on
annotation('textarrow',[0.516283524904214 0.619731800766283],...
    [0.406531531531531 0.288288288288288],'String','Attack Termination Point');
hold on
text(0.15,-0.055,'${x}_{52}$','FontSize',17,'interpreter','latex')
hold on
text(0.075,-0.07,'${x}_{154}$','FontSize',17,'interpreter','latex')
hold on
text(-0.04,0.03,'$\hat\mathcal{T}^0_e$','FontSize',17,'interpreter','latex')
hold on
text(0.103,-0.1,'$r_{100}$','FontSize',17,'interpreter','latex')
hold on
text(0.103,0.03,'$r_1$','FontSize',17,'interpreter','latex')
hold on
text(0.01,0,'$x_0$','FontSize',17,'interpreter','latex')
hold on
annotation('textarrow',[0.953065134099614 0.931992337164748],...
    [0.530405405405402 0.471846846846843],'String','$\hat{\mathcal{T}}^{40}_e$',...
    'Interpreter','latex',...
    'FontSize',14);
hold on
handle_tracking = plot(x_data(1,1:40),x_data(2,1:40),'b-','LineWidth',1.5)
hold on
handle_attack_tracking = plot(x_data(1,40:52),x_data(2,40:52),'r-','LineWidth',1.5)
hold on
handle_em_control = plot(x_data(1,52:62),x_data(2,52:62),'color'...
    ,[0.1294    0.6588    0.0588],'LineWidth',2.5,'LineStyle','-')
hold on
handle_tracking_after = plot(x_data(1,62:80),x_data(2,62:80),'m--','LineWidth',1)
hold on
handle_tracking = plot(x_data(1,100:125),x_data(2,100:125),'b-','LineWidth',1.5)
hold on 
handle_em_control = plot(x_data(1,145:151),x_data(2,145:151),'color'...
    ,[0.1294    0.6588    0.0588],'LineWidth',2.5,'LineStyle','-')
hold on 
handle_attack_tracking = plot(x_data(1,151:154),x_data(2,151:154),'r-','LineWidth',1.5)
hold on 
plot(x_data(1,154),x_data(2,154),'O','MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r','LineWidth',0.8)
hold on 
handle_tracking_after = plot(x_data(1,154:160),x_data(2,154:160),'m--','LineWidth',1)
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
plot(x_data(1,52),x_data(2,52),'O','MarkerSize',5,'MarkerEdgeColor','black',...
    'MarkerFaceColor','red','LineWidth',0.8)
hold on
plot(Td{40},[1 2], 'k-', 'LineWidth', 1)
xlim([-0.05 0.16])
ylim([-0.12 0.09])
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
set(ax,'units','normalized','position',[0.1,0.1,0.25,0.25])
box(ax,'on')
hold on

for i=40:40
    plot(Td{i},[1 2],'color','black','LineStyle','-','LineWidth',0.1)
    hold on
end
hold on
handle_tracking = plot(x_data(1,1:40),x_data(2,1:40),'b-','LineWidth',1.5)
hold on
handle_attack_tracking = plot(x_data(1,40:52),x_data(2,40:52),'r-','LineWidth',1.5)
hold on
handle_em_control = plot(x_data(1,52:62),x_data(2,52:62),'color'...
    ,[0.1294    0.6588    0.0588],'LineWidth',2.5,'LineStyle','-')
hold on
handle_tracking_after = plot(x_data(1,62:80),x_data(2,62:80),'m-','LineWidth',1)
hold on
handle_ref = plot(x_data(1,52),x_data(2,52),'O','MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r','LineWidth',0.8)

hold on
% plot(x_plus_data{51},[1 2],'b--','LineWidth',1)
% hold on
plot(x_plus_data{52},[1 2],'r--','LineWidth',1)
hold on
plot(x_plus_data{51},[1 2],'Color',[0.1294    0.6588    0.0588],'LineWidth',1,'LineStyle','--')
hold on
annotation('textarrow',[0.241379310344827 0.281609195402298],...
    [0.203828828828828 0.221846846846846],'String','$\hat{\mathcal{S}}^{+}_{51}$','FontSize',...
    14,'interpreter','latex');
hold on
annotation('textarrow',[0.326628352490421 0.328544061302682],...
    [0.158783783783784 0.190315315315315],'String','$\hat{\mathcal{S}}^{+}_{52}$','FontSize',...
    14,'interpreter','latex');
hold on
annotation('textarrow',[0.326867816091953 0.315134099616858],...
    [0.313907657657655 0.289414414414414],'String','$\hat{\mathcal{T}}^{40}_e$',...
    'Interpreter','latex',...
    'FontSize',14);
hold on

hold on 
text(0,0.0,'Safety Verification','FontSize',11,'Position',[0.112 -0.095 0])
xticks([])
yticks([])
xlim([0.1 0.16])
ylim([-0.1 -0.02])
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Box 2 %%%%%%%%%%%%%%%%%%%%%%
ax=axes;
set(ax,'units','normalized','position',[0.745,0.75,0.23,0.23])
box(ax,'on')
xticks([])
yticks([])
xlim([0.06 0.16])
ylim([-0.12 -0.02])

hold on

handle_tracking = plot(x_data(1,100:125),x_data(2,100:125),'b-','LineWidth',1.5)
hold on 
handle_em_control = plot(x_data(1,145:151),x_data(2,145:151),'color'...
    ,[0.1294    0.6588    0.0588],'LineWidth',2.5,'LineStyle','-')
hold on 
handle_attack_tracking = plot(x_data(1,151:155),x_data(2,151:155),'r-','LineWidth',1.5)
hold on 
plot(x_data(1,155),x_data(2,155),'O','MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerFaceColor','r','LineWidth',0.8)
hold on 
handle_tracking_after = plot(x_data(1,155:160),x_data(2,155:160),'m--','LineWidth',1)
hold on 
plot(x_plus_data{145},[1 2],'r--','LineWidth',1)
hold on
plot(x_plus_data{146},[1 2],'Color',[0.1294    0.6588    0.0588],'LineWidth',1,'LineStyle','--')

hold on 
for i=40:40
    plot(Td{i},[1 2],'color','black','LineStyle','-','LineWidth',0.1)
    hold on
end

hold on
annotation('textarrow',[0.817049808429119 0.832375478927203],...
    [0.890765765765765 0.832207207207207],'String','$\hat{\mathcal{S}}^{+}_{144}$','FontSize',...
    14,'interpreter','latex');
hold on
annotation('textarrow',[0.89176245210728 0.870689655172413],...
    [0.865990990990991 0.820945945945946],'String','$\hat{\mathcal{S}}^{+}_{145}$','FontSize',...
    14,'interpreter','latex');
hold on 
text(0,0.0,'Safety Verification','FontSize',11,'Position',[0.085 -0.026 0])
print -depsc -tiff -r300 -painters case_A_state_trajectory.eps
%% Visualization of the alarm signal & safety check signal and emergency signal  
% close all

alarm = alarm_data(1:201);
active = emergency(1:201);
safety = safety_data(1:201);

k = 1:201;
f = figure;
f.OuterPosition = [400 400 600 300]
ax = axes;
ax.OuterPosition = [-0.03  -0.03  1.06   1.1];
%%%%%%%%%%%%%%%%%%%%
x = [40 60 60 40];
y = [-0.9 -0.9 1.1 1.1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
%
hold on
x = [145 154 154 145];
y = [-0.9 -0.9 1.1 1.1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
colororder({'k','k'})
yyaxis left
xlim([0 200])
handle_alarm = plot(k(30:160),alarm(30:160),'r-','LineWidth',1.5);
ylim([-0.1 1.1])
xlim([30 160])
names = {'Normal'; 'Anomaly'};
set(gca,'ytick',[0:1],'yticklabel',names)
% xticks([])
hold on
yyaxis right
handle_emergency = plot(k(30:160),active(30:160),'b--','LineWidth',1.5);
hold on
handle_safety = plot(k(30:160),safety(30:160),'m--','LineWidth',1);
ylim([-0.1 1.1])
names = {'False'; 'True'};
set(gca,'ytick',[0:1],'yticklabel',names)
xlabel('$k$','interpreter','latex','FontSize',15)

xticks([0 20 40 52 60 80 100 120 145 155 180 200])
legend([handle_alarm,handle_emergency,handle_safety,handle_attack],...
    '$D_k$','emergency','$\hat{\mathcal{S}}^{+}_k\not \subseteq \mathcal{X}_{\eta} \lor u^{\prime}_k \notin \mathcal{U}$',...
    'attack',...
    'interpreter','latex','Position',[0.4 0.6 0.28 0.33],...
    'FontSize',11)
box on

print -depsc -tiff -r400 -painters case_A_alarm_safety.eps
%% Visualization of the control signals 

% close all
f2 = figure;
f2.OuterPosition = [400 100 600 700]

f21 = subplot(3,1,1);

x = [40 60 60 40];
y = [-1 -1 1 1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
x = [145 154 154 145];
y = [-1 -1 1 1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
handle_u = plot(u_ver(1,:),'b-','LineWidth',1);
hold on 
handle_c = yline(0.6111,'r--','LineWidth',2);
hold on 
handle_c = yline(-0.7778,'r--','LineWidth',2);
ylim([-0.8,0.65]);
yticks([-0.7778,0,0.6111])
xticks([])
xlim([0 200])
ylabel('$u_p$','interpreter','latex','FontSize',15);
legend([handle_c,handle_attack],...
    '$\mathcal{U}$','attack',...
    'interpreter','latex','Location','NorthEast',...
    'FontSize',11);
box on


%
f22 = subplot(3,1,2)
% ax2 = axes;
% f22.OuterPosition = [-0.02     -0.03     1.1     1.05]
x = [40 60 60 40];
y = [-2 -2 1 1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
x = [145 154 154 145];
y = [-2 -2 1 1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
hold on 
plot(u_ver(2,:),'b-','LineWidth',1)
hold on 
handle_c = yline(0.75,'r--','LineWidth',2);
hold on 
handle_c = yline(-1.25,'r--','LineWidth',2);
ylim([-1.3,0.8])
yticks([-1.25,0,0.75])
xticks([])
xlim([0 200])
ylabel('$u_l$','interpreter','latex','FontSize',15);
legend([handle_c,handle_attack],...
    '$\mathcal{U}$','attack',...
    'interpreter','latex','Location','NorthEast',...
    'FontSize',11);
box on

f23 = subplot(3,1,3)
% ax3 = axes;
% f23.OuterPosition = [-0.02     -0.03     1.1     1.05]
x = [40 60 60 40];
y = [-2 -2 1 1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
x = [145 154 154 145];
y = [-2 -2 1 1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
plot(u_ver(3,:),'b-','LineWidth',1)
hold on 
handle_c = yline(0.5235,'r--','LineWidth',2);
hold on 
handle_c = yline(-1.4765,'r--','LineWidth',2);
xlim([0 200])
ylim([-1.5,0.6])
yticks([-1.4765,0,0.5235])

xlabel('$k$','interpreter','latex','FontSize',15);
ylabel('$u_u$','interpreter','latex','FontSize',15);
legend([handle_c,handle_attack],...
    '$\mathcal{U}$','attack',...
    'interpreter','latex','Location','NorthEast',...
    'FontSize',11);
box on
print -depsc -tiff -r300 -painters case_A_control.eps

%% Comparison: system with the safety modules vs system with the safety modules 

% close all
f = figure;
f.OuterPosition = [400 100 1000 600];

f1 = subplot(2,2,1);
f1.OuterPosition = [0,0.5,0.47,0.45];
x = [40 60 60 40];
y = [-0.1 -0.1 0.5 0.5];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
x = [145 154 154 145];
y = [-0.15 -0.15 0.5 0.5];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
handle_x = plot(x_data(1,:),'Color','blue','LineWidth',2.5,'LineStyle',':')

hold on 
handle_ref = plot(ref(1,:),'r--','LineWidth',2) 
hold on 
handle_x_w = plot(1:200,x_w(1,1:200),'m:','LineWidth',1.7);

xlim([0 200]);
ylim([0 0.26]);
xticks([])
ylabel('$ x_1$','interpreter','latex','FontSize',15);
box on


f2 = subplot(2,2,3);
f2.OuterPosition = [0,0.02,0.47,0.5];
x = [40 60 60 40];
y = [-0.15 -0.15 0.5 0.5];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
box on
x = [145 154 154 145];
y = [-0.15 -0.15 0.5 0.5];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
handle_x = plot(x_data(2,:),'Color','blue','LineWidth',2.5,'LineStyle',':')

hold on 
handle_ref = plot(ref(2,:),'r--','LineWidth',2)
hold on 
handle_x_w = plot(1:200,x_w(2,1:200),'m:','LineWidth',1.7);
hold on 
xlim([0 200]);
ylim([-0.13 0.05]);
xticks([0 20 40 60 80 100 120 140 160 180 200])
xlabel('$k$','interpreter','latex','FontSize',15);
ylabel('$ x_2$','interpreter','latex','FontSize',15);

legend([handle_x,handle_x_w,handle_ref,handle_attack],...
    'With Safety Moduels','Without Safety Moduels','Reference','Attack period',...
    'FontSize',9,'Position',[0.275906382204392 0.357630918490216 0.180803574505565 0.136167714588321]);


f3 = subplot(2,2,[2 4]);
f3.OuterPosition = [0.47,0,0.55,1];
plot(Td{1}.mptPolytope.P,'Alpha',0.7,...
    'color','green','LineWidth',0.01)
hold on
for i=2:40
    plot(Td{i},[1 2],'color',[0.5020    0.5020    0.5020],'LineStyle','-','LineWidth',0.1)
    hold on
end
hold on 
handle_proposed = plot(x_data(1,1:200),x_data(2,1:200),'b:','LineWidth',2)
hold on 
handle_tracking = plot(x_w(1,1:200),x_w(2,1:200),'m:','LineWidth',1.5)
hold on
plot(x_w(1,54),x_w(2,54),'ro','LineWidth',1)
hold on
text(x_w(1,55),x_w(2,55),'$ x_{55}$','FontSize',15,'interpreter','latex')
hold on 
plot(x_w(1,61),x_w(2,61),'ro','LineWidth',1)
hold on
text(0.185,-0.105,'$ x_{60}$','FontSize',15,'interpreter','latex')

xlim([-0.1 0.25])
ylim([-0.13 0.05])


hold on
text(0.001,-0.001,'$ x_0$','FontSize',17,'interpreter','latex')
hold on 
for i=40:40
    plot(Td{i},[1 2],'color','black','LineStyle','-','LineWidth',0.1)
    hold on
end

hold on
text(0.08,-0.105,'$r_{100}$','FontSize',15,'interpreter','latex')
hold on
text(0.1,0.037,'$r_1$','FontSize',15,'interpreter','latex')

hold on
handle_ref1 = plot(ref(1,1),ref(2,1),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',0.8)
hold on
handle_ref2 = plot(ref(1,120),ref(2,120),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',1)
hold on 
plot(x_w(1,146),x_w(2,146),'ro','LineWidth',1)
hold on 
plot(x_w(1,155),x_w(2,155),'ro','LineWidth',1)
hold on 
text(0.23,-0.115,'$ x_{155}$','FontSize',15,'interpreter','latex')
hold on 
plot(x_data(1,62),x_data(2,62),'ro','LineWidth',1)
hold on 
text(0.115,-0.11,'$ x_{147}$','FontSize',15,'interpreter','latex')
hold on 
text(0.055,-0.036,'$ x_{62}$','FontSize',15,'interpreter','latex')
hold on 
plot(x_data(1,151),x_data(2,151),'ro','LineWidth',1)
hold on 
text(0.03,-0.057,'$ x_{151}$','FontSize',15,'interpreter','latex')

for i=40:40
    plot(Td{i},[1 2],'color','black','LineStyle','-','LineWidth',0.1)
    hold on
end

text(0.125,0.02,'$\hat\mathcal{T}^{40}_e$','FontSize',15,'interpreter','latex')

xlabel('${x}_1$','interpreter','latex','FontSize',15)
ylabel('${x}_2$','interpreter','latex','FontSize',15)

legend([handle_proposed,handle_tracking],'With Safety Moduels','Without Safety Moduels',...
    'Location','NorthWest','FontSize',10);



print -depsc -tiff -r400 -painters case_A_comparison.eps
%% Saving variables 

save x_data
save x_w
save u_ver 

%------------- END OF CODE --------------
