
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

rand('seed',20);
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
    else i>=100 && i<=sim_time
        ref(:,i)=[0.1;-0.1];
    end
end

% defining attack on the actuation channel
for i=1:sim_time
    if i>=150 && i<160
        u_a(:,i)=[0;0;0];
    else
        u_a(:,i)=[0;0;0];
    end
end


% defining attack on the measurement channel
for i=1:sim_time
    if i>=95 && i<=113
        y_a(:,i)=[0.0025*(i-95);0];
    else
        y_a(:,i)=[0;0];
    end
end

% computing data-driven tracking controller
K_data = data_driven_controller(A,B,W,X,U);

% initialization
alarm_data(1) = 0;
x_data(:,1) = [0.01;-0.01];  % initial state
x_data_prime(:,1) = [0.01;-0.01]; % initial state
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
f.Position = [700 70 800 700];

plot(x_data(1,1),x_data(2,1),'*','MarkerSize',4,'MarkerEdgeColor','k');
hold on

text(0.1,0.06,'$\hat{\mathcal{T}}^{40}_e$','FontSize',17,'interpreter','latex');

hold on
plot(Td{1}.mptPolytope.P,'Alpha',0.7,...
    'color','green','LineWidth',0.01);
hold on
for i=2:40
    plot(Td{i},[1 2],'color',[0.5020    0.5020    0.5020],'LineStyle','-','LineWidth',1);
    hold on
end
xlim([-0.1 0.5]);
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
x_w_prime(:,1) = [0.01;-0.01];
x_w(:,1) = [0.01;-0.01];

% simulation of the system -- presence of attack on the measurement channel
for k=1:sim_time
    
    % computing data-driven tracking controller
    ctr_data(:,k) = K_data*(x_data_prime(:,k)-ref(:,k)) + u_eq(:,k);
    ctr_data(:,k) =  min(max(ctr_data(:,k), [-0.7778;-1.25;-1.4765]), [0.611;0.75;0.5235]);
    
    
    % computing the one-step evolution set, \hat{\mathcal{R}}_k
    x = zonotope(x_data_prime(:,k),0*diag(ones(dim_x,1)));
    u = zonotope(ctr_data(:,k),0*diag(ones(dim_u,1)));
    x_pre_data{k} = AB * (cartProd(x,u))+ W;
    
    % attack on actuation (Note: in this scenario we dont have any attacks
    % on the actuation channel)
    ctr_data_prime(:,k) = ctr_data(:,k) + u_a(:,k);
    %     u_w_prime(:,k) = u_w(:,k) + u_a(:,k);
    
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
        index_data(k) = set_index(x_data(:,k),Td); % computing set membership index
        u_ver(:,k) = one_step_ctrl(3, x_data(:,k), Td_aug, index_data(k)); % computing D-ST-MPC controller
        u_ver(:,k) = min(max(u_ver(:,k), [-0.7778;-1.25;-1.4765]), [0.611;0.75;0.5235]);
        emergency(k) = 1; % emergency controller activation singal
    else
        u_ver(:,k) = ctr_data_prime(:,k); % applying tracking controller to the plant
        emergency(k) = 0; % emergency controller activation singal
    end
    
    % computing the system evolution
    x_data(:,k+1) = A*x_data(:,k) + B*u_ver(:,k) + randPoint(W);
    
    % attack on measurement channel
    x_data_prime(:,k+1) = x_data(:,k+1) + y_a(:,k);
    
    % visualization of the system evolution
    plot(x_data(1,k+1),x_data(2,k+1),'o','MarkerSize',2,'MarkerEdgeColor','r','MarkerFaceColor','r');
    
    hold on
    %     plot(x_w(1,k+1),x_w(2,k+1),'o','MarkerSize',2,'MarkerEdgeColor','b','MarkerFaceColor','b');
    
    % anomaly detector
    alarm_data(k+1) = detector_data_driven(x_data_prime(:,k+1),x_pre_data{k});
    
    % send a signal from detector to plant in the case of attack on the
    % measurement channel
    if alarm_data(k+1) == 1
        flag = 1;
    end
    %
    if norm(y_a(:,k))==[0;0]
        attack(k)=false;
    else
        attack(k)=true;
    end
    %
    title('Time: '+string(k)+', flag='+string(flag)+', emergency='+string(emergency(k))+...
        ', ignore='+string(ignore)+', D_k='+string(alarm_data(k))+', attack on y_k='+ string(attack(k)),...
        'FontSize',11,...
        'FontName','Courier')
    pause(0.1)
    k
end

        %%% Simulation of the system without the safety modules %%%%%%%%%%
hold on 
x_w(:,1) = [0.01;-0.01]; % system state without the safety modules 
x_w_prime(:,1) = [0.01;-0.01];
for k=1:sim_time
    u_w(:,k) = K_data*(x_w_prime(:,k)-ref(:,k)) + u_eq(:,k);
    u_w(:,k) =  min(max(u_w(:,k), [-0.7778;-1.25;-1.4765]), [0.611;0.75;0.5235]);
    x_w(:,k+1) = A*x_w(:,k) + B*u_w(:,k) + randPoint(W);
    x_w_prime(:,k+1) = x_w(:,k+1) + y_a(:,k);
    
    % attack on measurement channel
    x_w(:,k+1) = x_w(:,k+1) + y_a(:,k);
    
    % visualization of the system evolution
    plot(x_w(1,k+1),x_w(2,k+1),'o','MarkerSize',2,'MarkerEdgeColor','b','MarkerFaceColor','b');
    pause(0.1)
    hold on 
    if Td{40}.contains(x_w(:,k)) == 1
        safety_vio(k) = 0;
    else
        safety_vio(k) = 1;
    end
    k
end
%% Visualization of the system evolution
f = figure;
ax = axes;
f.OuterPosition = [200 50 850 800]
ax.OuterPosition = [-0.06     -0.04     1.15     1.11]

plot(Td{1}.mptPolytope.P,'Alpha',0.7,...
    'color','green','LineWidth',0.01)
hold on
for i=2:40
    plot(Td{i},[1 2],'color',[0.5020    0.5020    0.5020],'LineStyle','-','LineWidth',1)
    hold on
end

hold on
x1 = [0.863984674329502 0.899425287356322];
y1 = [0.809684684684685 0.692567567567567];
x2 = [0.797892720306513 0.909961685823755];
y2 = [0.245495495495495 0.412162162162162];

annotation('textarrow',x1,y1,'String','Attack Starting Point','FontSize',12)
hold on
annotation('textarrow',x2,y2,'String','Attack Detection Point','FontSize',12)
hold on
annotation('textarrow',[0.938697318007663 0.913793103448276],...
    [0.124 0.177927927927928],'String',{'$\hat\mathcal{T}^{40}_e$'},'FontSize',17,'interpreter','latex');
hold on
text(-0.03,0.08,'$\hat\mathcal{T}^0_e$','FontSize',17,'interpreter','latex')
hold on
text(0.01,0,'$x_0$','FontSize',17,'interpreter','latex')
hold on
hold on
text(0.103,-0.1,'$r_{100}$','FontSize',17,'interpreter','latex')
hold on
text(0.103,0.03,'$r_1$','FontSize',17,'interpreter','latex')
hold on
handle_ref1 = plot(ref(1,1),ref(2,1),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',0.8)
hold on
handle_ref2 = plot(ref(1,120),ref(2,120),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',1)
hold on
handle_tracking = plot(x_data(1,1:94),x_data(2,1:94),'b-','LineWidth',1.5)
hold on
handle_attack_tracking = plot(x_data(1,95:111),x_data(2,95:111),'r-','LineWidth',1.5)
hold on
handle_em_control = plot(x_data(1,111:115),x_data(2,111:115),'color'...
    ,[0.1294    0.6588    0.0588],'LineWidth',...
    2,'LineStyle','-')
hold on
handle_tracking_after = plot(x_data(1,115:127),x_data(2,115:127),'m--','LineWidth',1.5)
hold on

handle_initial_state = plot(x_data(1,1),x_data(2,1),'o','MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',0.8)
hold on
plot(x_data(1,111),x_data(2,111),'O','MarkerSize',5,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.1294    0.6588    0.0588],'LineWidth',0.8);
hold on
text(0.105,-0.05,'$x_{110}$','FontSize',17,'interpreter','latex');
hold on
plot(Td{40},[1 2], 'k-', 'LineWidth', 1);
hold on

xlim([-0.11 0.12])
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
%%%%%%%%%%%%%%%%%%%%%% box [anomaly detector] %%%%%%%%
box on
ax=axes;
set(ax,'units','normalized','position',[0.15,0.15,0.24,0.25])
box(ax,'on')
hold on
plot(x_pre_data{110},[1 2],'b','LineWidth',1)
hold on
plot(x_data_prime(1,111),x_data_prime(2,111),'r*')
hold on
text(0.138,-0.05,'$\hat{\mathcal{R}}^+_{109}$','FontSize',14,'interpreter','latex')
hold on
text(0.1428,-0.05,'$x^{\prime}_{110}$','FontSize',14,'interpreter','latex')
box on
xlim([0.133 0.145])
ylim([-0.057 -0.045])
annotation('textbox',...
    [0.203107279693487 0.371058558558559 0.156088122605364 0.027027027027027],...
    'String',{'Anomaly Detector'},...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

print -depsc -tiff -r300 -painters case_B_state_trajectory.eps
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
x = [95 113 113 95];
y = [-1 -1 2.1 2.1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
%
hold on
colororder({'k','k'})
yyaxis left

handle_alarm = plot(k(80:150),alarm(80:150),'r-','LineWidth',1.5);
ylim([-0.1 1.1])

names = {'Normal'; 'Anomaly'};
set(gca,'ytick',[0:1],'yticklabel',names)
xticks([])
hold on
yyaxis right
handle_emergency = plot(k(80:150),active(80:150),'b--','LineWidth',1.5);
hold on
handle_safety = plot(k(80:150),safety(80:150),'m--','LineWidth',1);
ylim([-0.1 1.1])
names = {'False'; 'True'};
set(gca,'ytick',[0:1],'yticklabel',names)
xlabel('$k$','interpreter','latex','FontSize',15)
xticks([80 90 95 105 110 115 130 140 150])
legend([handle_alarm,handle_emergency,handle_safety,handle_attack],...
    '$D_k$','emergency','$\hat{\mathcal{S}}^{+}_k\not \subseteq \mathcal{X}_{\eta} \lor u^{\prime}_k \notin \mathcal{U}$',...
    'attack',...
    'interpreter','latex','Location','NorthEast',...
    'FontSize',11)
xlim([80 150])
box on
print -depsc -tiff -r300 -painters case_B_alarm_safety.eps;
%%
% close all
f2 = figure;
f2.OuterPosition = [400 100 600 700]
f21 = subplot(3,1,1);

x = [95 113 113 95];
y = [-1 -1 1 1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on
handle_u = plot(u_ver(1,:),'b-','LineWidth',1);
hold on
handle_c = yline(0.6111,'r--','LineWidth',2);
hold on
handle_c = yline(-0.7778,'r--','LineWidth',2);
ylim([-0.8,0.65]);
xlim([0 200]);
yticks([-0.7778,0,0.6111])
xticks([])
ylabel('$u_p$','interpreter','latex','FontSize',15);
legend([handle_c,handle_attack],...
    '$\mathcal{U}$','attack',...
    'interpreter','latex','Location','NorthEast',...
    'FontSize',11);
box on


%
f22 = subplot(3,1,2)

x = [95 113 113 95];
y = [-1.3 -1.3 1 1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(u_ver(2,:),'b-','LineWidth',1)
hold on
handle_c = yline(0.75,'r--','LineWidth',2);
hold on
handle_c = yline(-1.25,'r--','LineWidth',2);
ylim([-1.3,0.8])
xlim([0 200]);
yticks([-1.25,0,0.75])
xticks([])
ylabel('$u_l$','interpreter','latex','FontSize',15);
legend([handle_c,handle_attack],...
    '$\mathcal{U}$','attack',...
    'interpreter','latex','Location','NorthEast',...
    'FontSize',11);
box on

f23 = subplot(3,1,3)

x = [95 113 113 95];
y = [-2 -2 1 1];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(u_ver(3,:),'b-','LineWidth',1)
hold on
handle_c = yline(0.5235,'r--','LineWidth',2);
hold on
handle_c = yline(-1.4765,'r--','LineWidth',2);

ylim([-1.5,0.6])
xlim([0 200]);
yticks([-1.4765,0,0.5235])
xticks([0 20 40 60 80 95 113 140 160 180 200])
xlabel('$k$','interpreter','latex','FontSize',15);
ylabel('$u_u$','interpreter','latex','FontSize',15);
legend([handle_c,handle_attack],...
    '$\mathcal{U}$','attack',...
    'interpreter','latex','Location','NorthEast',...
    'FontSize',11);
box on
print -depsc -tiff -r300 -painters case_B_control.eps

%%
% close all
f = figure;
f.OuterPosition = [400 100 1000 600];

f1 = subplot(2,2,1);
f1.OuterPosition = [0,0.5,0.47,0.45];

hold on 
x = [95 113 113 95];
y = [-0.15 -0.15 0.5 0.5];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
handle_x = plot(x_data(1,:),'Color','blue','LineWidth',2.5,'LineStyle',':')
hold on 
handle_ref = plot(ref(1,:),'r--','LineWidth',2)
hold on 
handle_x_w = plot(1:200,x_w(1,1:200),'m:','LineWidth',1.7);
hold on
xlim([0 200]);
ylim([0 0.5]);
xticks([])
ylabel('${x}_1$','interpreter','latex','FontSize',15);
% legend([handle_x,handle_x_w_s,handle_x_w_a,handle_ref,handle_attack],...
%     'proposed method','nominal system (safe)','nominal system (unsafe)','$r_k$',...
%     'interpreter','latex','Location','NorthEast',...
%     'FontSize',10,'Location','NorthWest');
box on
legend([handle_x,handle_x_w,handle_ref,handle_attack],...
    'With Safety Moduels','Without Safety Moduels','Reference','Attack period',...
    'Location','NorthWest',...
    'FontSize',10);

f2 = subplot(2,2,3);
f2.OuterPosition = [0,0.02,0.47,0.5];
x = [95 113 113 95];
y = [-0.15 -0.15 0.5 0.5];
handle_attack = patch(x,y,'red','FaceAlpha',0.2,'EdgeColor','none');
hold on 
handle_x = plot(x_data(2,:),'Color','blue','LineWidth',2.5,'LineStyle',':');
hold on 
handle_ref = plot(ref(2,:),'r--','LineWidth',2);
hold on 
handle_x_w = plot(1:200,x_w(2,1:200),'m:','LineWidth',1.7);
hold on 
xlim([0 200]);
ylim([-0.13 0.05]);
xticks([0 20 40 60 80 90 101 113 140 160 180 200])
xlabel('$k$','interpreter','latex','FontSize',15);
ylabel('${x}_2$','interpreter','latex','FontSize',15);

box on 

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
text(0.11,-0.1,'$r_{100}$','FontSize',17,'interpreter','latex')
hold on
text(0.1,0.037,'$r_1$','FontSize',17,'interpreter','latex')
hold on
text(0.01,0,'$x_0$','FontSize',17,'interpreter','latex')

hold on
handle_ref1 = plot(ref(1,1),ref(2,1),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',0.8)
hold on
handle_ref2 = plot(ref(1,120),ref(2,120),'pentagram','MarkerSize',9,'MarkerEdgeColor','k',...
    'MarkerFaceColor','b','LineWidth',1)
hold on 
plot(x_w(1,101),x_w(2,101),'ro','LineWidth',0.7);
hold on 
text(0.13,0.025,'$x_{101}$','FontSize',15,'interpreter','latex')
hold on 
plot(x_w(1,114),x_w(2,114),'ro','LineWidth',1);
hold on 
plot(x_w(1,173),x_w(2,173),'ro','LineWidth',1);
hold on 
text(0.44,-0.015,'$x_{113}$','FontSize',15,'interpreter','latex')
hold on 
text(0.17,-0.06,'$x_{173}$','FontSize',15,'interpreter','latex')
hold on
for i=40:40
    plot(Td{i},[1 2],'color','black','LineStyle','-','LineWidth',0.1)
    hold on
end
hold on 
text(0.14,-0.01,'$\hat\mathcal{T}^{40}_e$','FontSize',15,'interpreter','latex')

hold on 
plot(x_data(1,111),x_data(2,111),'ro','LineWidth',1);
hold on 
text(0.1,-0.055,'$x_{111}$','FontSize',15,'interpreter','latex')
hold on 
plot(x_data(1,115),x_data(2,115),'ro','LineWidth',1);
hold on 
text(0.02,-0.04,'$x_{115}$','FontSize',15,'interpreter','latex')

xlabel('${x}_1$','interpreter','latex','FontSize',15)
ylabel('${x}_2$','interpreter','latex','FontSize',15)

xlim([-0.1 0.48])
ylim([-0.11 0.045])

legend([handle_proposed,handle_tracking],'With Safety Moduels','Without Safety Moduels',...
    'Location','NorthEast','FontSize',10);

print -depsc -tiff -r400 -painters case_B_comparison.eps
%% Saving data 

x_data_caseB = 'x_data.mat';
save x_data_caseB;

x_w_caseB = 'x_w.mat';
save x_w_caseB;

u_ver_caseB = 'u_ver.mat';
save u_ver_caseB;

%------------- END OF CODE --------------
