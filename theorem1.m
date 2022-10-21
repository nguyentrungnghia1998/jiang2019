%% Time and step
clc; clear; close all;
Step = 0.001;
T_end = 15;
t = 0:Step:T_end;
data = cell(1,length(t));
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
Q = 2*eye(4);
R = 0.2*eye(2);
P = icare(A,[B C],Q,blkdiag(R,-gamma^2));
%% Initial value
x{1} = [0.2;0;0.8;0];
%% Simulation
for i = 1:length(t)
    x_i = x{i};
    ro_a{i} = [4*cos(t(i))*sin(x_i(2))*x_i(1); 5*sin(t(i))*sin(x_i(4))*x_i(2)];
    delta_w{i} = 3*sin(x_i(2))*x_i(1);
    u{i} = -pinv(R)*B'*P*x{i};

    if i == length(t)
        break
    end

    % Update
    x{i+1} = x{i} + Step*(A*x{i} + B*(u{i} + ro_a{i}) + C*delta_w{i});
end

x_m = cell2mat(x);
plot(t,x_m);