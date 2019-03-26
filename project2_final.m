%% Project 2 - simulated single phase microgrid
% Belongs to MAE284, Robust and Multivariable Control, 2018
% Written by R.A. de Callafon, Dept. of MAE <callafon@ucsd.edu>
clear all
close all
clc
% nominal parameters
R1 = 1;     % [Ohm]
C1 = 0.1;   % [F]
L1 = 1;     % [H] 
R2 = 2;     % [Ohm]
C2 = 0.5;   % [F]
L2 = 1;     % [H]

% nominal model
P11=tf([R1*L1*C1 1],[L1*C1 R1*C1 1]);
P12=tf([1 0],[L1*C1 R1*C1 1]);
P21=tf([L1*C1 0],[L1*C1 R1*C1 1]);
P22=tf([L1*C1 1],[L1*C1 R1*C1 1]);
P=tf([R2 L2*C2 1],[R2*L2*C2 L2*C2 1])*[P11 P12;P21 P22];
Pnom=minreal(ss(P));

model variation/uncertainty
L1nom=L1;
delta=[-1 1];
for k=1:2,
    L1=L1nom*(1+0.25*delta(k));
    P11=tf([R1*L1*C1 1],[L1*C1 R1*C1 1]);
    P12=tf([1 0],[L1*C1 R1*C1 1]);
    P21=tf([L1*C1 0],[L1*C1 R1*C1 1]);
    P22=tf([L1*C1 1],[L1*C1 R1*C1 1]);
    P=tf([R2 L2*C2 1],[R2*L2*C2 L2*C2 1])*[P11 P12;P21 P22];
    P=minreal(ss(P));
    figure(2)
    bodemag(P,'g');
    hold on
end
bodemag(Pnom,'r')
hold off
grid

size(Pnom)

%% Control Design and Stability Analysis

P0 = tf([R2 L2*C2 1],[R2*L2*C2 L2*C2 1])*[P11 P12;P21 P22]; % nominal model
P0_ss = balreal(minreal(ss(P0)));

Wa = 0.02*tf([L2*C2 1],1)*[P11 , P12 ; P21 , P22];
Ws = 0.01*tf([1 , 0 ; 0 , 1]);
% LFT

Q11 = [tf(zeros(2,2)) , tf(zeros(2,2)) ; Ws , Ws];
Q12 = [Wa ; Ws*P0];
Q21 = [tf(eye(2,2)) tf(eye(2,2))];
Q22 = [P0];
Q = [Q11 Q12 ; Q21 Q22];
Q_ss = ss(Q);
Q_ss = minreal(Q_ss);

% H_inf
gamma = 0;
[gopt, K] = hinflmi(ltisys(Q_ss.A,Q_ss.B,Q_ss.C,Q_ss.D),[2 2],gamma,1e-10);
[Ac,Bc,Cc,Dc] = ltiss(K);
eig(Ac);

K_s = tf(ss(Ac,Bc,Cc,Dc));
K_ss = minreal(ss(K_s));

[hb,g] = balreal(K_ss);
g';   % to see what are weakly coupled to the input and output. So that You can then delete these states.

%Singular Values of Sensitivity Function
figure(1)
S = (eye(2,2) - P0*K_s)^-1;
sigma(minreal(ss(S)))
title('Singular Values of S')
grid on
saveas(gcf,'S_sigma.png')

%Singular Values of Complementary Sensitivity Function
figure(2)
T = P0*K_s*(eye(2,2) - P0*K_s)^-1;
sigma(minreal(ss(T)))
title('Singular Values of T')
grid on
saveas(gcf,'T_sigma.png')

%Singular Values of M11
figure(3)
M = lft(Q_ss,K_s);
M11 = M(1:2,1:2);
sigma(minreal(ss(M11)))
grid on
title('Singular Values of M11')
saveas(gcf,'M11_sigma.png')

%Singular Values of M
figure(4)
sigma(minreal(lft(Q,K_s)))
title('Singular Values of M')
grid on
saveas(gcf,'M_sigma.png')

% Robust Stability
rbst_s = Wa*K_s*(eye(2,2) - P0*K_s)^-1;
disp('robust stability')
norm(rbst_s,inf)

% Nominal performance
nmnl_p = Ws*((eye(2)+P0*(eye(2,2) - P0*K_s))^-1);
disp('nominal performance')
norm(nmnl_p,inf)

% Robust performance
rbst_p = [rbst_s;nmnl_p];
disp('robust performance')
norm(rbst_p,inf)

%Prepare for Export
A_k = K_ss.A;
B_k = K_ss.B;
C_k = K_ss.C;
D_k = K_ss.D;
save('Project2_sp18_K','A_k','B_k','C_k','D_k')
return