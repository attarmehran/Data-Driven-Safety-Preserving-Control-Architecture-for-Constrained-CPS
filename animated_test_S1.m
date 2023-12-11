
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

if isfile('Td.mat')
     disp('ROSC sets ')
else
     disp('no')
end

Td = load('Td.mat').Td;
V_AB = load('V_AB.mat').V_AB;
AB = load('AB.mat').AB;
K = load('K.mat').K;
Td_aug = load('Td_aug.mat').Td_aug;


sim_time = 150;
A=[0.993 0.003;0.007 0.982];
B=[0.008 -0.003 -0.003;0 0.003 0.003];
C=eye(2);
D=0;

dim_x = size(A,1);
dim_u = size(B,2);
X = zonotope(interval([-0.48;-0.48],[0.3;0.3]));
U = zonotope(interval([-0.7778;-1.25;-1.4765],[0.611;0.75;0.5235]));
W = zonotope(zeros(dim_x,1),0.001*eye(dim_x));
% plot(U)
sys = ss(A,B,C,D);
% [V_AB,AB] = compute_AB(sys,X,U,W);


for i=1:size(V_AB,1)
    A_hat{i} = V_AB{i}(1:size(A,1),1:size(A,1));
    B_hat{i} = V_AB{i}(1:size(A,1),size(A,1)+1:size(A,1)+size(B,2));
end
% reference
for i=1:sim_time
    if i<100
        ref(:,i)=[0.1;0.03];
    elseif i>=100 & i<=200
        ref(:,i)=[0.1;-0.1];
    end
end
%
for i=1:sim_time
    if i>=95 & i<=113
        u_a(:,i)=[1;1;2];
%         u_a(:,i)=[0;0;0];
    else
        u_a(:,i)=[0;0;0];
    end
end

SYS = ss(A,B,C,D);
K_data = data_driven_controller(A,B,W,X,U);
x_pred = {};
alarm_data(1) = 0;
x_model_w(:,1) = [0.01;-0.01];
x_data_w(:,1) = [0.01;-0.01];
x_data(:,1) = [0.01;-0.01];
x_prime_model(:,1) = [0.1;-0.1];

for i=1:sim_time
    u_eq(:,i) = pinv(B)*((eye(dim_x)-A)*ref(:,i));
end
%
close all
f = figure;
ax = axes;
f.Position = [100 60 1400 700];

subplot(3,2,[1 5])

plot(x_data(1,1),x_data(2,1),'*','MarkerSize',4,'MarkerEdgeColor','k')
hold on
% hold on 
annotation('textarrow',[0.438964285714286 0.418428571428571],...
    [0.876095238095241 0.829666666666667],'String','$\hat{\mathcal{T}}^{40}_e$',...
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
% hold on
% text(0.098,0.055,'$x_{104}$','FontSize',17,'interpreter','latex')
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
% hold on

flag = 0;
ignore = 0;
emergency(1) = 0;
pause(1)
hold on 
%
h = text(-56,5,'attack free','EdgeColor','green','BackgroundColor','green');
hold on 
h2 = text(-57,-0.6,'tracking controller is active','interpreter','Latex',...
           'FontSize',14,'Color','black','EdgeColor','green','BackgroundColor','green');
hold on 
%
for k=1:sim_time
   
   %
   ctr_data(:,k) = K_data*(x_data(:,k)-ref(:,k)) + u_eq(:,k);
   ctr_data(:,k) =  min(max(ctr_data(:,k), [-0.7778;-1.25;-1.4765]), [0.611;0.75;0.5235]);
   
   x = zonotope(x_data(:,k),0*diag(ones(dim_x,1)));
   u = zonotope(ctr_data(:,k),0*diag(ones(dim_u,1)));
   x_pre_data{k} = AB * (cartProd(x,u))+ W;
   
   % attack on actuation 
   ctr_data_prime(:,k) = ctr_data(:,k) + u_a(:,k);
   ctr_data_prime(:,k) =  min(max(ctr_data_prime(:,k),...
       [-0.7778;-1.25;-1.4765]), [0.611;0.75;0.5235]);
   
   % attacker changes the flag 
   if alarm_data(k)==1 & emergency(k-1) == 0
       flag = 0;
   end
   %
   
   if flag == 1 & Td{1}.contains(x_data(:,k)) == 1
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
    %

    if flag == 1
%         u_ver(:,k) = data_driven_em_ctr(x_data(:,k),Td,A_hat,B_hat,K,W,A,B,U);
        index_data(k) = set_index(x_data(:,k),Td);
        u_ver(:,k) = one_step_ctrl(3, x_data(:,k), Td_aug, index_data(k));
        emergency(k) = 1;
        pause(0.001)
    else
        u_ver(:,k) = ctr_data_prime(:,k);
        emergency(k) = 0;
        pause(0.001)
    end
    %
    % system evolution 
    x_data(:,k+1) = A*x_data(:,k) + B*u_ver(:,k) + randPoint(W);
    
    subplot(3,2,[1 5])
    plot(x_data(1,k+1),x_data(2,k+1),'o','MarkerSize',2,'MarkerEdgeColor','r','MarkerFaceColor','r');
    title('Time: '+string(k-1)+', flag='+string(flag)+...
        ', ignore='+string(ignore),...
        'FontSize',14,...
    'FontName','Courier')
    hold on 
    
                    %%%%%%% BOX %%%%%%%
    if k>=80 & k<=105
        
        box on
        ax=axes;
        set(ax,'units','normalized','position',[0.16 0.15 0.181 0.3])
        box(ax,'on')
        hold on
%         plot(Td{1}.mptPolytope.P,'Alpha',0.7,...
%         'color','green','LineWidth',0.01)
%         hold on
        for i=40:40
            plot(Td{i},[1 2],'color','black','LineStyle','-','LineWidth',0.1)
            hold on
        end
        hold on 
        
        xlim([0.06 0.11])
        ylim([-0.02 0.08])
        
        hold on
   
        if k==105
            plot(x_plus_data{k},[1 2],'r-','LineWidth',1)
            hold on
            plot(x_data(1,k),x_data(2,k),...
            'o','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
            hold on 
            t = text(x_data(1,k)+0.003,x_data(2,k)+0.003,['$\hat{\mathcal{S}}^{+}_{' int2str(k-1) '}$'],...
               'interpreter','Latex','FontSize',11);
            hold on
            pause(0.001)
        else
            plot(x_plus_data{k},[1 2],'b-','LineWidth',1)
            hold on
            plot(x_data(1,k+1),x_data(2,k+1),...
                'o','MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r');
            hold on 
            t = text(x_data(1,k)+0.003,x_data(2,k)+0.003,...
                ['$\hat{\mathcal{S}}^{+}_{' int2str(k-1) '}$'],...
               'interpreter','Latex','FontSize',11);
            pause(0.001)
            hold on
        end         
        hold on
    end
    % detector 
    alarm_data(k+1) = detector_data_driven(x_data(:,k+1),x_pre_data{k});
    
    if alarm_data(k+1)==1
        flag = 1;
    end
    %
    if norm(u_a(:,k))==[0;0;0]
       attack(k)=false;
   else
       attack(k)=true;
    end
   %
   hold on 
   
   %
    
    subplot(3,2,2)
    plot(k,alarm_data(k),'r*','MarkerSize',7)
    xlim([0 sim_time])
    ylim([-0.1 1.1])
    names = {'Normal'; 'Anomaly'};
    set(gca,'ytick',[0:1],'yticklabel',names)
    ylabel('$D_k$','Interpreter','Latex','FontSize',20)
    xlabel('$k$','Interpreter','Latex','FontSize',20)
    hold on 
    title('Anomaly Detector','FontSize',11,'FontName','Courier')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[]%%%%%
    subplot(3,2,4)
    plot(k,emergency(k),'r*','MarkerSize',7)
    xlim([0 sim_time])
    ylim([-0.1 1.1])
    names = {'True'; 'False'};
    set(gca,'ytick',[0:1],'yticklabel',names)
    ylabel('emergency','Interpreter','tex','FontSize',14)
    xlabel('$k$','Interpreter','Latex','FontSize',20)
    title('E-DSTC','FontSize',11,'FontName','Courier')
    
    hold on 
    subplot(3,2,6)
    plot(k,safety_data(k),'r*','MarkerSize',7)
    xlim([0 sim_time])
    ylim([-0.1 1.1])
    names = {'True'; 'False'};
    set(gca,'ytick',[0:1],'yticklabel',names)
    ylabel('$\hat{\mathcal{S}}^{+}_k\not \subseteq \mathcal{X}_{\eta} \lor u^{\prime}_k \notin \mathcal{U}$',...
        'Interpreter','Latex','FontSize',14)
    xlabel('$k$','Interpreter','Latex','FontSize',20)
    title('Safety Verification','FontSize',11,'FontName','Courier')
    
    hold on 
    delete(h)
    %
    if attack(k)==false
        delete(h)
       h = text(-54,5,'attack free','interpreter','Latex',...
           'FontSize',14,'Color','black','EdgeColor','green','BackgroundColor','green');
%         hold on 
%         h = text(-57,-0.6,'tracking controller is active','interpreter','Latex',...
%            'FontSize',14,'Color','black','EdgeColor','green','BackgroundColor','green');
        
    else
       delete(h)
       h = text(-54,5,'attack on $u_k$','interpreter','Latex',...
           'FontSize',14,'Color','black','EdgeColor','red','BackgroundColor','red');

    end 
    hold on 
    if emergency(k) == 1
        delete(h2)
       h2 = text(-57,-0.6,'emergency controller is active','interpreter','Latex',...
           'FontSize',14,'Color','black','EdgeColor','red','BackgroundColor','red');
       hold on 
    else 
       delete(h2)
       h2 = text(-57,-0.6,'tracking controller is active','interpreter','Latex',...
           'FontSize',14,'Color','black','EdgeColor','green','BackgroundColor','green'); 
       hold on 
    end
   
    pause(0.1)
    k
end
%------------- END OF CODE --------------