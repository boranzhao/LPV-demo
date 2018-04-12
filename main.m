%% Main fucntion for LPV controller design demo
% This code is for two scheduling parameters, but can also be used for one GS. For the later, set Theta2 to 0 and d_Theta2{2}  to 0. 

clear;
% clc;
%% Plant parameters  
% use the example in the paper [Wu 96, Induced L2 norm control for LPV systems with bounded parameter variation rates]
% Note that disturbance signal dn is ignored. 
% generate the generalized plant
% Note that the gamma value computed may be a little different from the one in the paper [Wu 96]. This is due to mainly two reasons
% (1) Practical validity is enforced here, but not in the paper [Wu 96]
% (2) The Lyapunov positivity condition [X I; I Y]> 0 is modified to be [X I; I Y] > eps*I in order to improve the numerical stability. 
% (3) options(3) for mincx is set to be a non-zero number to avoid overly large numbers in the computed decision variables. 
syms theta1 
Gsym.A = [0.75 2 cos(theta1) sin(theta1) 0 0;...
     0 0.5 -sin(theta1) cos(theta1) 0 0;...
     0 0 -10 0 0 0;...
     0 0 0 -10 0 0;...
     0 0 0 0 -0.2 0;
     0 0 0 0 0 -0.2
     ];
Gsym.B1 = [0 0; 0 0; 0 0; 0 0;4 0;0 4];
Gsym.B2 = [0 0; 0 0;10 0; 0 10;0 0; 0 0];
Gsym.C1 = [1 0 0 0 -5 0;...
   0 1 0 0 0 -5;...
   0 0 0 0 0 0;...
   0 0 0 0 0 0];
Gsym.D11 = zeros(4,2);
Gsym.D12 = [0 0; 0 0;1/280 0;0 1/280];
Gsym.C2 = [1 0 0 0 -5 0;0 1 0 0 0 -5];
Gsym.D21 = zeros(2,2);
Gsym.D22 = zeros(2,2);

Theta1 = [-pi pi]; % for non-switching LPV control
Theta2 = 0;
d_Theta = {[-10 10],0}; %Droumax and Droumin

% grid the parameter sets
Theta1Grid  =  [-pi:pi/5:pi];
Theta2Grid = 0;

%% Hinf controller design for comparison purposes for a fixed point theta1= 0
G0 = AugPltEval(Gsym,[0 0]');
sys = ss(G0.A, [G0.B1 G0.B2],[G0.C1;G0.C2], [G0.D11 G0.D12; G0.D21 G0.D22]);
[Khinf,~,Gam_opt,~]= hinfsyn(sys,2,2);
% return;
%% Design parameters for LPV controllers
% note that for designing practically valid LPV controllers, only one of X and Y can be parameter-dependent. Both cases should be tested for best performance.
XY_PD = 2;          % 0 for parameter-independent X and Y, 1 for parameter-dependent X and constant Y; 2 for parameter-dependent Y and constant X
SolverSel = 1;      % 1 for Matlab LMI lab, 2 for CVX(Sedmui); currently the code for CVX is not complete
IsProject = 0;      % boolean, 0 for using the basic conditions (Theorem 2.1 in [Apkarian 1998], 1 for using the projected condition (Theorem 2.2 in [Apkarian 1998]
StrucUncert = 0;    % boolean, 0 for normal LPV synthesis, 1 for LPV synthesis considering the structured uncertainty (Section V in [Apkarian 98])

% Fcn_theta: function handle, Fcn_theta(theta) will give all the scalar functions for the matrix variables. To realize X = X0+ f1(theta)X1+f2(theta)X2, set Fcn_theta(theta) = [1 f1(theta) f2(theta)]
% d_Fcn_theta: function handle for the derivative of the scalar functions
% FthetaNum: 1 by 3 vector in the form of [1 n1 n2], where n1 (n2) denotes the number of scalar functions depending on theta1 (theta2)
Fcn_theta = @(x) [1 cos(x(1)) sin(x(1))]; %function for PD matrices
d_Fcn_theta = @(x) [0 -sin(x(1)) cos(x(1))];
FthetaNum = [1 2 0];

%% solve the LPV  synthesis problem
% collect the parameters 
plant_paras.Theta1Grid = Theta1Grid;
plant_paras.Theta2Grid = Theta2Grid;
plant_paras.d_Theta = d_Theta;
plant_paras.StrucUncert = StrucUncert;

design_paras.XY_PD = XY_PD;
design_paras.Fcn_theta = Fcn_theta;
design_paras.d_Fcn_theta = d_Fcn_theta;
design_paras.FthetaNum = FthetaNum;
design_paras.IsProject = IsProject;
design_paras.SolverSel = SolverSel;

[Gam_opt,Kopt, Xopt,Yopt] = LPV(Gsym,plant_paras,design_paras);
%% for simulation in Simulink
Ah = Kopt.Ah;
Bh = Kopt.Bh;
Ch = Kopt.Ch;
Dh = Kopt.Dh;
