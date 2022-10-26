%% Time and step
clc; clear; close all;
Step = 0.001;
T_end = 15;
t = 0:Step:T_end;
data = cell(1,length(t));
x_TH = cell(1,3);
u_TH = cell(1,3);
for k = 1:3
    %% Load weight
    load("weight.mat");
    W1 = double(W1);
    W2 = double(W2);
    W3 = double(W3);
    if k == 1
        W = W1;
    elseif k == 2
        W = W2;
    else
        W = W3;
    end
    %% Variable
    x = data;
    u = data;
    ro_a = data;
    delta_w = data;
    %% Parameter
    Bp = 0.8;
    By = 0.318;
    Kpp = 0.204;
    Kyy = 0.072;
    Kpy = 0.0068;
    Kyp = 0.0219;
    Jp = 0.0178;
    Jy = 0.0084;
    l = 0.186;
    m = 1.3872;
    gamma = 5;
    J_Tp = Jp + m*l^2;
    J_Ty = Jy + m*l^2;
    A = [0 0 1 0;
         0 0 0 1;
         0 0 -Bp/J_Tp 0;
         0 0 0 -By/J_Ty];
    B = [0 0;0 0; Kpp/J_Tp Kpy/J_Tp; Kyp/J_Ty Kyy/J_Ty];
    C = [1 0 1 0]';
    R = 0.2*eye(2);
    %% Initial value
    x{1} = [0.2;0.8;0;-0.1];
    %% Simulation
    for i = 1:length(t)
        x_i = x{i};
        ro_a{i} = [4*cos(t(i))*sin(x_i(2))*x_i(1); 5*sin(t(i))*sin(x_i(4))*x_i(2)];
        delta_w{i} = 3*sin(x_i(2))*x_i(1);
        u{i} = -1/2*pinv(R)*B'*gradientActivate(x{i})'*W;
    
        if i == length(t)
            break
        end
    
        % Update
        x{i+1} = x{i} + Step*(A*x{i} + B*(u{i} + ro_a{i}) + C*delta_w{i});
    end
    
    x_TH{k} = cell2mat(x);
    u_TH{k} = cell2mat(u);
end

%% Plot
figure(1);
x_TH1 = x_TH{1};
x_TH2 = x_TH{2};
x_TH3 = x_TH{3};
plot(t,x_TH1(1,:),'r-',t,x_TH1(3,:),'k-',t,x_TH2(1,:),'-.b',t,x_TH2(3,:),'-.c',t,x_TH3(1,:),':m',t,x_TH3(3,:),':g','LineWidth',2);
legend("Q=6I_4:x_1","Q=6I_4:x_3","Q=2I_4:x_1","Q=2I_4:x_3","Q=0.5I_4:x_1","Q=0.5I_4:x_3");
title("Dynamics of the system states x_1 and x_3")
grid on;

figure(2);
plot(t,x_TH1(2,:),'r-',t,x_TH1(4,:),'k-',t,x_TH2(2,:),'-.b',t,x_TH2(4,:),'-.c',t,x_TH3(2,:),':m',t,x_TH3(4,:),':g','LineWidth',2);
legend("Q=6I_4:x_2","Q=6I_4:x_4","Q=2I_4:x_2","Q=2I_4:x_4","Q=0.5I_4:x_2","Q=0.5I_4:x_4");
title("Dynamics of the system states x_2 and x_4")
grid on;

figure(3);
u_TH1 = u_TH{1};
u_TH2 = u_TH{2};
u_TH3 = u_TH{3};
plot(t,u_TH1(1,:),'r-',t,u_TH1(2,:),'k-',t,u_TH2(1,:),'-.b',t,u_TH2(2,:),'-.c',t,u_TH3(1,:),':m',t,u_TH3(2,:),':g','LineWidth',2);
legend("Q=6I_4:u_1","Q=6I_4:u_2","Q=2I_4:u_1","Q=2I_4:u_2","Q=0.5I_4:u_1","Q=0.5I_4:u_2");
title("Dynamics of the control inputs u_1 and u_2")
grid on;
%% Sub function
function a = gradientActivate(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
a = [2*x1 0 0 0;
     x2 x1 0 0;
     x3 0 x1 0;
     x4 0 0 x1;
     0 2*x2 0 0;
     0 x3 x2 0;
     0 x4 0 x2;
     0 0 2*x3 0;
     0 0 x4 x3;
     0 0 0 2*x4];
end