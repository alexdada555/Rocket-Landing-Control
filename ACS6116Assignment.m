close all;
clear;
clc;

%% Model Parameters

T = 0.5;

I = eye(3,3);
O = zeros(3,3);

theta = 10*pi/180;
phi = 30*pi/180;

%% Model Matrices

Adt = [I T*I;
       O I ];

Bdt = [0.5*(T^2)*I;
       T*I];
   
Cdt = eye(6,6);

Ddt = zeros(6,6);

%%  Pre Analysis

poles = eig(Adt);
figure(1);
plot(poles,'*');
title('OpenLoop Eigenvalues');

Cont = rank(ctrb(Adt,Bdt));

%% Discrete-Time Full State-Feedback Control

Q = diag([0.05,0.05,0.05,1,1,1]); % State penalty
R = eye(3,3);  % Control penalty

Kdt = dlqr(Adt,Bdt,Q,R,0); % DT State-Feedback Controller Gains

Cl_poles = eig((Adt - Bdt*Kdt));
figure(2)
plot(Cl_poles,'*');
title('ClosedLoop Eigenvalues');

%%  LQR Dynamic Simulation

Time = 100;
kT = round(Time/T);

X = zeros(6,kT);
U = zeros(3,kT);

X(:,1) = [600;-600;500;0;0;-5];

for k = 1:kT-1

    %Control    
    U(:,k) = -Kdt*X(:,k);
    
    %Simulation
    X(:,k+1) = Adt*X(:,k) + Bdt*U(:,k);
end

%% Plots

t = (0:kT-1)*T;

figure(3);
subplot(2,1,1);
plot(t,X([1,2,3],:));
legend('x','y','z');
title('Position LQR');
ylabel('Meters(m)');

subplot(2,1,2);
plot(t,X([4,5,6],:));
legend('vx','vy','vz');
title('Velocity LQR');
ylabel('Meters/Seconds m/s');

figure(4);
plot(t,U);
legend('fx','fy','fz - mg');
title('Force/Mass LQR');
ylabel('Newtons/Meter (N/Kg)');

%%   Terminal State Penalty

Phi = (Adt - Bdt*Kdt);  % Closed Loopp System
S = (Q + Kdt' * R * Kdt);
P = dlyap(Phi',S);  % Terminal State Penalty from Lyapunov Function

%%  Definition of MPC Prediction, Cost and Constaint Matrices

N = 5;  % Prediction Horizon

[F,G] = predict_mats(Adt,Bdt,N); % State Parameter Matrices
[H,L,M] = cost_mats(F,G,Q,R,P);  % State Cost Function Matricies

 %Contraints  
 Px = [eye(6); -1*eye(6)];
 qx = [inf;inf;inf;20;20;15;
       inf;inf;inf;20;20;15]; % State Constraints  
 Pxf = Px;
 qxf = qx; % Terminal State Constraints
 Pu = [eye(3); -1*eye(3)];
%% Open Loop Paradigm Dynamic Simulation

Time = 100;
kT = round(Time/T);

X = zeros(6,kT);
U = zeros(3,kT);
Cost = zeros(1,kT);

X(:,1) = [600;-600;500;5;5;-5];
X(:,2) = [600;-600;500;5;5;-5];

Ref = [0;0;0;0;0;0];

for k = 2:kT-1
    
    %Contraints
    qu = [(U(3,k-1)+9.81)*tan(theta);(U(3,k-1)+9.81)*tan(theta);12-9.81;
          (U(3,k-1)+9.81)*tan(theta);(U(3,k-1)+9.81)*tan(theta);9.81]; % Input Constraints

    [Pc, qc, Sc] = constraint_mats(F,G,Pu,qu,Px,qx,Pxf,qxf); % Constraints As Linear Inequality 

    %Control    
    [Useq,Cost(:,k)] = quadprog(H,L*(X(:,k)-Ref),Pc,qc + Sc*X(:,k)); %Solve Quadratic program
    U(:,k) =  Useq(1:3);
    
    %Simulation    
    X(:,k+1) = Adt*X(:,k) + Bdt*U(:,k);
end

%% Plots

t = (0:kT-1)*T;

figure(5);
subplot(2,1,1);
plot(t,X([1,2,3],:));
legend('x','y','z');
title('Position OL LQ-MPC');
ylabel('Meters(m)');

subplot(2,1,2);
plot(t,X([4,5,6],:));
legend('vx','vy','vz');
title('Velocity OL LQ-MPC');
ylabel('Meters/Seconds (m/s)');

figure(6);
plot(t,U);
legend('fx','fy','fz - mg');
title('Force/Mass OL LQ-MPC');
ylabel('Newtons/Meter (N/Kg)');

figure(7);
plot(t,Cost);
title('Cost OL LQ-MPC');

%%  Definition of Closed-Loop MPC Prediction, Cost and Constaint Matrices

% N = 5;  % Prediction Horizon
% 
% [Fxcl,Gxcl,Fycl,Gycl,Fucl,Gucl] = predict_mats_cl(Phi,Bdt,Cdt,Kdt,N); % State Parameter Matrices
% [Sc1,Scx] = cost_mats_cl(Fxcl,Gxcl,Fycl,Gycl,Fucl,Gucl,Q,R,P); % State Cost Function Matricies
% 
% %% Closed Loop Paradigm Dynamic Simulation
% 
% Time = 100;
% kT = round(Time/T);
% 
% X = zeros(6,kT);
% U = zeros(3,kT);
% Y = zeros(6,kT);
% c = zeros(N*3,kT);
% Cost = zeros(1,kT);
% 
% X(:,1) = [600;-600;500;0;0;-5];
% X(:,2) = [600;-600;500;0;0;-5];
% 
% for k = 2:kT-1
%     
%     %Control    
%     [c(:,k),Cost(:,k)] = quadprog(Sc1,Scx*X(:,k)); % Solve Quadratic program
%     U(:,k) =  -Kdt*X(:,k) + c(1:3,k);
%     
%     %Simulation    
%     X(:,k+1) = Adt*X(:,k) + Bdt*U(:,k);
% end
% 
% %% Plots
% 
% t = (0:kT-1)*T;
% 
% figure(8);
% subplot(2,1,1);
% plot(t,X([1,2,3],:));
% legend('x','y','z');
% title('Position CL LQ-MPC');
% ylabel('Meters(m)');
% 
% subplot(2,1,2);
% plot(t,X([4,5,6],:));
% legend('vx','vy','vz');
% title('Velocity CL LQ-MPC');
% ylabel('Meters/Seconds (m/s)');
% 
% figure(9);
% subplot(2,1,1);
% plot(t,U);
% legend('fx','fy','fz - mg');
% title('Force/Mass CL LQ-MPC');
% ylabel('Newtons/Meter (N/Kg)');
% 
% subplot(2,1,2);
% plot(t,c(1:3,:));
% legend('fx~','fy~','fz - mg~');
% title('Force/Mass peterbations CL LQ-MPC');
% ylabel('Newtons/Meter Adjustments (N/Kg)');
% 
% figure(10);
% plot(t,Cost);
% title('Cost CL LQ-MPC');
