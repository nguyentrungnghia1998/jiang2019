%% Parameter
clc; clear; close all;
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
save('P_matrix.mat','P');