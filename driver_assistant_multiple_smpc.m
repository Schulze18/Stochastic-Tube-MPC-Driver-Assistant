%% Stochastic Predictive Control for Semi-Autonomous Vehicles with an Uncertain Driver Model
%  Andrew Gray ; Yiqi Gao ; Theresa Lin ; J. Karl Hedrick ; Francesco Borrelli 
% DOI: 10.1109/ITSC.2013.6728575
clear all
close all
% clc

% Car Parameters 
mass = 2050;
Iz = 3344;
Calfa = 80000;
Calfa_front = Calfa;
Calfa_rear = Calfa;
lf = 1.43;
lr = 1.47;
u_coef = 1;
Ts = 50*10^-3;
Vx = 15;
Ky = 0.005;
Kpsi = 0.2;
tlp = 2;

%%
Ac = [0 1 0 0;
      0 -(2*Calfa_front + 2*Calfa_rear)/(mass*Vx) (2*Calfa_front + 2*Calfa_rear)/(mass) -(2*Calfa_front*lf + 2*Calfa_rear*lr)/(mass*Vx);
      0 0 0 1;
      0 -(2*Calfa_front*lf - 2*Calfa_rear*lr)/(Iz*Vx) (2*Calfa_front*lf - 2*Calfa_rear*lr)/(Iz) -(2*Calfa_front*lf^2 + 2*Calfa_rear*lr^2)/(Iz*Vx)];

Bc = [0;
      2*Calfa_front/mass;
      0;
       2*Calfa_front/Iz];
   
Ec = [0;
     -(2*Calfa_front*lf - 2*Calfa_rear*lr)/(mass*Vx)-Vx;
     0;
     -(2*Calfa_front*lf^2 + 2*Calfa_rear*lr^2)/(Iz*Vx)];
 
Nstate = size(Ac,1);
Ncontrol = size(Bc,2); 

F = -[Ky 0 Kpsi 0];
Ac_driver = Ac + Bc*F;
Bc_driver = Bc;
Ec_driver = [Ec -Bc*Kpsi];
Dc_driver = Bc;
[Ad,Bd]= c2d(Ac_driver, Bc_driver, Ts);
Dd = Bd;

%% LQR
Qlqr = 0.1*eye(4);
Rlqr = 50*eye(Ncontrol);
Klqr = lqrd(Ac_driver, Bc_driver, Qlqr, Rlqr, Ts);

%% Controler Parameters
N = 15;
w_covar = 0.1;
w_bar = 0;
contraint_p = 0.02;
vmax = 0.4;
vmin = -0.4;
epsilon_max_limits = 2.5;
epsilon_min_limits = -2.5;

g_const = [1 0 0 0];
pd = makedist('Normal','mu', w_bar,'sigma',1);
% qi_p = icdf(pd, 1-contraint_p);
qi_p_1 = icdf(pd, 1-0.01);
qi_p_2 = icdf(pd, 1-0.02);
qi_p_5 = icdf(pd, 1-0.05);
qi_p_8 = icdf(pd, 1-0.08);
w_max_rmpc = 3;
% qi_p = sqrt((1-contraint_p)/contraint_p);
% sigma_contraint = sqrt(g_const*Bd*w_covar*Bd'*g_const')*qi_p;

%% Optimization Problem Matrices
PHI = Ad - Bd*Klqr;
Hz = eye(Nstate);
Hc = zeros(Nstate,Ncontrol*(N));
Hw = zeros(Nstate,Ncontrol*(N));

Kbar = Klqr;
G_const = g_const;
E_prob = [Dd zeros(4,Ncontrol*(N-1))];
for i = 1:N-1
    Hz = [Hz; PHI^i];
    Kbar = blkdiag(Kbar, Klqr);
    G_const = blkdiag(G_const, g_const);
    
    Hc_temp = [];
    Hw_temp = [];
    E_prob_temp = [];
    for j = 1:N
        if i >= j
            Hc_temp = [Hc_temp PHI^(i-j)*Bd];
            Hw_temp = [Hw_temp PHI^(i-j)*Dd];
%             E_prob_temp = [E_prob_temp PHI^(i-j-1)*Dd];
        else
            Hc_temp = [Hc_temp zeros(Nstate,Ncontrol)];
            Hw_temp = [Hw_temp zeros(Nstate,Ncontrol)];
%             E_prob_temp = [E_prob_temp zeros(Nstate,Ncontrol)];
        end
        
        if i >= j-1
            E_prob_temp = [E_prob_temp PHI^(i-j+1)*Dd];
        else
            E_prob_temp = [E_prob_temp zeros(Nstate,Ncontrol)];
        end
    end
    Hc = [Hc; Hc_temp];
    Hw = [Hw; Hw_temp];
    E_prob = [E_prob; E_prob_temp];
end

for i = 1:N
   sigma_contraint_1(i,1) = sqrt(g_const*E_prob((i-1)*Nstate+1:i*Nstate,:)*w_covar*E_prob((i-1)*Nstate+1:i*Nstate,:)'*g_const')*qi_p_1; 
   sigma_contraint_2(i,1) = sqrt(g_const*E_prob((i-1)*Nstate+1:i*Nstate,:)*w_covar*E_prob((i-1)*Nstate+1:i*Nstate,:)'*g_const')*qi_p_2; 
   sigma_contraint_5(i,1) = sqrt(g_const*E_prob((i-1)*Nstate+1:i*Nstate,:)*w_covar*E_prob((i-1)*Nstate+1:i*Nstate,:)'*g_const')*qi_p_5; 
   sigma_contraint_8(i,1) = sqrt(g_const*E_prob((i-1)*Nstate+1:i*Nstate,:)*w_covar*E_prob((i-1)*Nstate+1:i*Nstate,:)'*g_const')*qi_p_8; 
   sigma_contraint_rmpc(i,1) = sqrt(g_const*E_prob((i-1)*Nstate+1:i*Nstate,:)*w_covar*E_prob((i-1)*Nstate+1:i*Nstate,:)'*g_const')*w_max_rmpc; 
end


% % Hqp = (Kbar*Hc + eye(N*Ncontrol))'*(Kbar*Hc + eye(N*Ncontrol));
% % Fqp = [(Hz'*Kbar'*(Hc'*Kbar') + Hz'*Kbar');
% %        (Hw'*Kbar'*(Hc'*Kbar') + Hw'*Kbar')];

Hqp = Hc'*Hc + (Kbar*Hc + eye(N*Ncontrol))'*(Kbar*Hc + eye(N*Ncontrol));
Fqp = [(Hz'*Kbar'*(Hc'*Kbar') + Hz'*Kbar' + Hz'*Hc);
       (Hw'*Kbar'*(Hc'*Kbar') + Hw'*Kbar' + Hw'*Hc)];
 
A_cmax = Kbar*Hc + eye(N*Ncontrol);
b_cmax = vmax*ones(N*Ncontrol,1);
A_cmin = -Kbar*Hc - eye(N*Ncontrol);
b_cmin = -vmin*ones(N*Ncontrol,1);
   
A_epislon_max = G_const*Hc;
b_epislon_max = epsilon_max_limits*ones(N, 1);
A_epislon_min = -G_const*Hc;
b_epislon_min = -epsilon_min_limits*ones(N, 1);

Aineq = [A_cmax;
         A_cmin;
         A_epislon_max;
         A_epislon_min];
bineq = [b_cmax;
         b_cmin;
         b_epislon_max;
         b_epislon_min]; 
    
   
%% Solver Setup
options = optimoptions('quadprog','Display','None');
% Create an OSQP object
prob_smpc_1 = osqp;
prob_smpc_2 = osqp;
prob_smpc_5 = osqp;
prob_smpc_8 = osqp;
prob_rmpc = osqp;
prob_mpc = osqp;


% Setup workspace
prob_smpc_1.setup(Hqp, zeros(1,15), Aineq, [], bineq, 'warm_start', false, 'verbose', false);
prob_smpc_2.setup(Hqp, zeros(1,15), Aineq, [], bineq, 'warm_start', false, 'verbose', false);
prob_smpc_5.setup(Hqp, zeros(1,15), Aineq, [], bineq, 'warm_start', false, 'verbose', false);
prob_smpc_8.setup(Hqp, zeros(1,15), Aineq, [], bineq, 'warm_start', false, 'verbose', false);
prob_rmpc.setup(Hqp, zeros(1,15), Aineq, [], bineq, 'warm_start', false, 'verbose', false);
prob_mpc.setup(Hqp, zeros(1,15), Aineq, [], bineq, 'warm_start', false, 'verbose', false);

%% System Simulation
Nsim = 200;
x_smpc_1 = zeros(4,Nsim);
x_smpc_2 = zeros(4,Nsim);
x_smpc_5 = zeros(4,Nsim);
x_smpc_8 = zeros(4,Nsim);
x_rmpc = zeros(4,Nsim);
x_mpc = zeros(4,Nsim);
x_driver = zeros(4,Nsim);
x_klqr = zeros(4,Nsim);

v_smpc_1 = zeros(1,Nsim);
v_smpc_2 = zeros(1,Nsim);
v_smpc_5 = zeros(1,Nsim);
v_smpc_8 = zeros(1,Nsim);
v_rmpc = zeros(1,Nsim);
v_mpc = zeros(1,Nsim);

c_smpc_1 = zeros(1,Nsim);
c_smpc_2 = zeros(1,Nsim);
c_smpc_5 = zeros(1,Nsim);
c_smpc_8 = zeros(1,Nsim);
c_rmpc = zeros(1,Nsim);
c_mpc = zeros(1,Nsim);
c_max_array = zeros(1,Nsim);
c_min_array = zeros(1,Nsim);

J_smpc_1 = zeros(1,Nsim);
J_smpc_2 = zeros(1,Nsim);
J_smpc_5 = zeros(1,Nsim);
J_smpc_8 = zeros(1,Nsim);
J_rmpc = zeros(1,Nsim);
J_mpc = zeros(1,Nsim);

t = (0:(Nsim-1))*Ts;
x0 = [-1; 0; 6*pi/180; 0.1];
x_smpc_1(:,1) = x0;
x_smpc_2(:,1) = x0;
x_smpc_5(:,1) = x0;
x_smpc_8(:,1) = x0;
x_rmpc(:,1) = x0;
x_mpc(:,1) = x0;
x_driver(:,1) = x0;
x_klqr(:,1) = x0;

%% Constraints
epsilon_max_array = epsilon_max_limits*ones(N + Nsim, 1);
epsilon_min_array = epsilon_min_limits*ones(N + Nsim, 1);

dist_constraint = 30;
time_constraint = dist_constraint/Vx;
n_constraint = time_constraint/Ts;
epsilon_min_array(n_constraint:(n_constraint+30)) = 0.0;
% % 
dist_constraint = 70;
time_constraint = dist_constraint/Vx;
n_constraint = time_constraint/Ts;
epsilon_max_array(n_constraint:(n_constraint+30)) = 0.0;
% % 
% % 
dist_constraint = 120;
time_constraint = dist_constraint/Vx;
n_constraint = time_constraint/Ts;
epsilon_min_array(n_constraint:(n_constraint+30)) = 0.0;
% % 
% % 
% % dist_constraint = 1750;
% % time_constraint = dist_constraint/Vx;
% % n_constraint = time_constraint/Ts;
% % epsilon_max_array(n_constraint:(n_constraint+150)) = 0.0;


% Variables
cont_violations_smpc_1 = 0;
cont_violations_smpc_2 = 0;
cont_violations_smpc_5 = 0;
cont_violations_smpc_8 = 0;
cont_violations_rmpc = 0;
cont_violations_mpc = 0;
cont_violations_lqr = 0;
cont_cant_u_smpc_1 = 0;
cont_cant_u_smpc_2 = 0;
cont_cant_u_smpc_5 = 0;
cont_cant_u_smpc_8 = 0;
cont_cant_u_rmpc = 0;
cont_cant_u_mpc = 0;


w_value = 0;
for i = 1:Nsim
    w_bar = 0;%F*x(:,i);
    
    Epsilon_max = epsilon_max_array(i:(i+N-1),1);
    Epsilon_min = epsilon_min_array(i:(i+N-1),1);
    
    
    W_array = w_bar*ones(N*Ncontrol,1);
    
    %%%%%%%%%%%%%%%%%%%% SMPC - 1%  %%%%%%%%%%%%%%%%  
    b_cmax = vmax*ones(N*Ncontrol,1) - [-Kbar*Hz Kbar*Hw]*[x_smpc_1(:,i); W_array];
    b_cmin = -vmin*ones(N*Ncontrol,1) + [-Kbar*Hz Kbar*Hw]*[x_smpc_1(:,i); W_array];
    b_epislon_max = Epsilon_max - [G_const*Hz G_const*Hw]*[x_smpc_1(:,i); W_array] - sigma_contraint_1;
    b_epislon_min = -Epsilon_min + [G_const*Hz G_const*Hw]*[x_smpc_1(:,i); W_array] - sigma_contraint_1;
    bineq = [b_cmax;
             b_cmin;
             b_epislon_max;
             b_epislon_min]; 
    prob_smpc_1.update('q',[x_smpc_1(:,i); W_array]'*Fqp,'u', bineq);
    res_smpc_1 = prob_smpc_1.solve();
    
    % Check solver status
    if ~strcmp(res_smpc_1.info.status, 'solved')
%         error('OSQP did not solve the problem!')
        cont_cant_u_smpc_1 = cont_cant_u_smpc_1 + 1;
        c_smpc_1(:,i) = 0;%c_smpc(:,i-1);
    else
        c_smpc_1(:,i) = res_smpc_1.x(1,1);
    end
    J_smpc_1(i) = res_smpc_1.x'*Hqp*res_smpc_1.x + [x_smpc_1(:,i); W_array]'*Fqp*res_smpc_1.x;
    
    W_array = w_bar*ones(N*Ncontrol,1);
    
    %%%%%%%%%%%%%%%%%%%% SMPC - 2%  %%%%%%%%%%%%%%%%  
    b_cmax = vmax*ones(N*Ncontrol,1) - [-Kbar*Hz Kbar*Hw]*[x_smpc_2(:,i); W_array];
    b_cmin = -vmin*ones(N*Ncontrol,1) + [-Kbar*Hz Kbar*Hw]*[x_smpc_2(:,i); W_array];
    b_epislon_max = Epsilon_max - [G_const*Hz G_const*Hw]*[x_smpc_2(:,i); W_array] - sigma_contraint_2;
    b_epislon_min = -Epsilon_min + [G_const*Hz G_const*Hw]*[x_smpc_2(:,i); W_array] - sigma_contraint_2;
    bineq = [b_cmax;
             b_cmin;
             b_epislon_max;
             b_epislon_min]; 
    prob_smpc_2.update('q',[x_smpc_2(:,i); W_array]'*Fqp,'u', bineq);
    res_smpc_2 = prob_smpc_2.solve();
    
    % Check solver status
    if ~strcmp(res_smpc_2.info.status, 'solved')
%         error('OSQP did not solve the problem!')
        cont_cant_u_smpc_2 = cont_cant_u_smpc_2 + 1;
        c_smpc_2(:,i) = 0;
    else
        c_smpc_2(:,i) = res_smpc_2.x(1,1);
    end
    J_smpc_2(i) = res_smpc_2.x'*Hqp*res_smpc_2.x + [x_smpc_2(:,i); W_array]'*Fqp*res_smpc_2.x;
    
    
    %%%%%%%%%%%%%%%%%%%% SMPC - 5%  %%%%%%%%%%%%%%%%  
    b_cmax = vmax*ones(N*Ncontrol,1) - [-Kbar*Hz Kbar*Hw]*[x_smpc_5(:,i); W_array];
    b_cmin = -vmin*ones(N*Ncontrol,1) + [-Kbar*Hz Kbar*Hw]*[x_smpc_5(:,i); W_array];
    b_epislon_max = Epsilon_max - [G_const*Hz G_const*Hw]*[x_smpc_5(:,i); W_array] - sigma_contraint_5;
    b_epislon_min = -Epsilon_min + [G_const*Hz G_const*Hw]*[x_smpc_5(:,i); W_array] - sigma_contraint_5;
    bineq = [b_cmax;
             b_cmin;
             b_epislon_max;
             b_epislon_min]; 
    prob_smpc_5.update('q',[x_smpc_5(:,i); W_array]'*Fqp,'u', bineq);
    res_smpc_5 = prob_smpc_5.solve();
    
    % Check solver status
    if ~strcmp(res_smpc_5.info.status, 'solved')
%         error('OSQP did not solve the problem!')
        cont_cant_u_smpc_5 = cont_cant_u_smpc_2 + 1;
        c_smpc_5(:,i) = 0;
    else
        c_smpc_5(:,i) = res_smpc_5.x(1,1);
    end
    J_smpc_5(i) = res_smpc_5.x'*Hqp*res_smpc_5.x + [x_smpc_5(:,i); W_array]'*Fqp*res_smpc_5.x;
    
    
    %%%%%%%%%%%%%%%%%%%% SMPC - 8%  %%%%%%%%%%%%%%%%  
    b_cmax = vmax*ones(N*Ncontrol,1) - [-Kbar*Hz Kbar*Hw]*[x_smpc_8(:,i); W_array];
    b_cmin = -vmin*ones(N*Ncontrol,1) + [-Kbar*Hz Kbar*Hw]*[x_smpc_8(:,i); W_array];
    b_epislon_max = Epsilon_max - [G_const*Hz G_const*Hw]*[x_smpc_8(:,i); W_array] - sigma_contraint_8;
    b_epislon_min = -Epsilon_min + [G_const*Hz G_const*Hw]*[x_smpc_8(:,i); W_array] - sigma_contraint_8;
    bineq = [b_cmax;
             b_cmin;
             b_epislon_max;
             b_epislon_min]; 
    prob_smpc_8.update('q',[x_smpc_8(:,i); W_array]'*Fqp,'u', bineq);
    res_smpc_8 = prob_smpc_8.solve();
    
    % Check solver status
    if ~strcmp(res_smpc_8.info.status, 'solved')
%         error('OSQP did not solve the problem!')
        cont_cant_u_smpc_8 = cont_cant_u_smpc_8 + 1;
        c_smpc_8(:,i) = 0;
    else
        c_smpc_8(:,i) = res_smpc_8.x(1,1);
    end
    J_smpc_8(i) = res_smpc_8.x'*Hqp*res_smpc_8.x + [x_smpc_8(:,i); W_array]'*Fqp*res_smpc_8.x;
    
    
     %%%%%%%%%%%%%%%%%%%% RMPC %%%%%%%%%%%%%%
    b_cmax = vmax*ones(N*Ncontrol,1) - [-Kbar*Hz Kbar*Hw]*[x_rmpc(:,i); W_array];
    b_cmin = -vmin*ones(N*Ncontrol,1) + [-Kbar*Hz Kbar*Hw]*[x_rmpc(:,i); W_array];
    b_epislon_max = Epsilon_max - [G_const*Hz G_const*Hw]*[x_rmpc(:,i); W_array] - sigma_contraint_rmpc;
    b_epislon_min = -Epsilon_min + [G_const*Hz G_const*Hw]*[x_rmpc(:,i); W_array] - sigma_contraint_rmpc;
    bineq = [b_cmax;
             b_cmin;
             b_epislon_max;
             b_epislon_min]; 
    prob_rmpc.update('q',[x_rmpc(:,i); W_array]'*Fqp,'u', bineq);
    res_rmpc = prob_rmpc.solve();
    
    % Check solver status
    if ~strcmp(res_rmpc.info.status, 'solved')
%         error('OSQP did not solve the problem!')
        cont_cant_u_rmpc = cont_cant_u_rmpc + 1;
        c_rmpc(:,i) = 0;%c_smpc(:,i-1);
    else
        c_rmpc(:,i) = res_rmpc.x(1,1);
    end
    J_rmpc(i) = res_rmpc.x'*Hqp*res_rmpc.x + [x_rmpc(:,i); W_array]'*Fqp*res_rmpc.x;
    
    %%%%%%%%%%%%%%%%%%%% MPC %%%%%%%%%%%%%%
%     b_cmax = vmax*ones(N*Ncontrol,1) - [Kbar*Hz Kbar*Hw]*[x_mpc(:,i); W_array];
%     b_cmin = -vmin*ones(N*Ncontrol,1) + [Kbar*Hz Kbar*Hw]*[x_mpc(:,i); W_array];
    b_cmax = vmax*ones(N*Ncontrol,1) - [-Kbar*Hz Kbar*Hw]*[x_mpc(:,i); W_array];
    b_cmin = -vmin*ones(N*Ncontrol,1) + [-Kbar*Hz Kbar*Hw]*[x_mpc(:,i); W_array];
    b_epislon_max = Epsilon_max - [G_const*Hz G_const*Hw]*[x_mpc(:,i); W_array]; 
    b_epislon_min = -Epsilon_min + [G_const*Hz G_const*Hw]*[x_mpc(:,i); W_array];
    bineq = [b_cmax;
             b_cmin;
             b_epislon_max;
             b_epislon_min]; 
    prob_mpc.update('q',[x_mpc(:,i); W_array]'*Fqp,'u', bineq);
    res_mpc = prob_mpc.solve();
    
    % Check solver status
    if ~strcmp(res_mpc.info.status, 'solved')
%         error('OSQP did not solve the problem!')
        cont_cant_u_mpc = cont_cant_u_mpc + 1;
        c_mpc(:,i) = 0;%c_mpc(:,i-1);
    else
        c_mpc(:,i) = res_mpc.x(1,1);
    end
    J_mpc(i) = res_mpc.x'*Hqp'*res_mpc.x + [x_mpc(:,i); W_array]'*Fqp*res_mpc.x;
    
   %%%%%%%%%%%%%%%%%%%% Model Update %%%%%%%%%%%%%%
   w_value(i) = normrnd(w_bar,(w_covar));

   v_smpc_1(:,i) = -Klqr*x_smpc_1(:,i) + c_smpc_1(:,i);
   v_smpc_2(:,i) = -Klqr*x_smpc_2(:,i) + c_smpc_2(:,i);
   v_smpc_5(:,i) = -Klqr*x_smpc_5(:,i) + c_smpc_5(:,i);
   v_smpc_8(:,i) = -Klqr*x_smpc_8(:,i) + c_smpc_8(:,i);
   v_rmpc(:,i) = -Klqr*x_rmpc(:,i) + c_rmpc(:,i);
   v_mpc(:,i) = -Klqr*x_mpc(:,i) + c_mpc(:,i);
   
   
   % SMPC State
   x_smpc_1(:,i+1) = Ad*x_smpc_1(:,i) + Bd*v_smpc_1(:,i) + Dd*w_value(i);
   x_smpc_2(:,i+1) = Ad*x_smpc_2(:,i) + Bd*v_smpc_2(:,i) + Dd*w_value(i);
   x_smpc_5(:,i+1) = Ad*x_smpc_5(:,i) + Bd*v_smpc_5(:,i) + Dd*w_value(i);
   x_smpc_8(:,i+1) = Ad*x_smpc_8(:,i) + Bd*v_smpc_8(:,i) + Dd*w_value(i);
   
   
   % SMPC State
   x_rmpc(:,i+1) = Ad*x_rmpc(:,i) + Bd*v_rmpc(:,i) + Dd*w_value(i);
   
   % MPC State
   x_mpc(:,i+1) = Ad*x_mpc(:,i) + Bd*v_mpc(:,i) + Dd*w_value(i);
   
   % KLR State
   x_klqr(:,i+1) = Ad*x_klqr(:,i) - Bd*Klqr*x_klqr(:,i) + Dd*w_value(i);
   
   % Driver
   x_driver(:,i+1) = Ad*x_driver(:,i);
   
   if(x_smpc_1(1,i) <= Epsilon_min(1,1) || x_smpc_1(1,i) >= Epsilon_max(1,1))
        cont_violations_smpc_1 = cont_violations_smpc_1 + 1; 
   end
   if(x_smpc_2(1,i) <= Epsilon_min(1,1) || x_smpc_2(1,i) >= Epsilon_max(1,1))
        cont_violations_smpc_2 = cont_violations_smpc_2 + 1; 
   end
   if(x_smpc_5(1,i) <= Epsilon_min(1,1) || x_smpc_5(1,i) >= Epsilon_max(1,1))
        cont_violations_smpc_5 = cont_violations_smpc_5 + 1; 
   end
   if(x_smpc_8(1,i) <= Epsilon_min(1,1) || x_smpc_8(1,i) >= Epsilon_max(1,1))
        cont_violations_smpc_8 = cont_violations_smpc_8 + 1; 
   end
   if(x_rmpc(1,i) <= Epsilon_min(1,1) || x_rmpc(1,i) >= Epsilon_max(1,1))
        cont_violations_rmpc = cont_violations_rmpc + 1; 
   end
   if(x_mpc(1,i) <= Epsilon_min(1,1) || x_mpc(1,i) >= Epsilon_max(1,1))
        cont_violations_mpc = cont_violations_mpc + 1; 
   end
   if(x_klqr(1,i) <= Epsilon_min(1,1) || x_klqr(1,i) >= Epsilon_max(1,1))
        cont_violations_lqr = cont_violations_lqr + 1; 
   end
   
end


%% Plot
x_array = Vx*t;

% % % %%
% % % figure
% % % plot(x_array, x_smpc(1,1:Nsim))
% % % hold on
% % % plot(x_array, x_mpc(1,1:Nsim),'k')
% % % plot(x_array, x_rmpc(1,1:Nsim),'g')
% % % plot(x_array, x_klqr(1,1:Nsim),'m')
% % % plot(x_array, zeros(Nsim),'--r')
% % % ylabel('Y [m]','FontSize',10)
% % % xlabel('X [m]','FontSize',10)
% % % legend('SMPC','MPC','RMPC','LQR')
% % % plot(x_array, epsilon_max_array(1:Nsim) ,'--r');
% % % plot(x_array, epsilon_min_array(1:Nsim) ,'--r');
% % % % plot(x_array, epsilon_max_array(1:Nsim)-sigma_contraint ,'--k');
% % % % plot(x_array, epsilon_min_array(1:Nsim)+sigma_contraint ,'--k');
% % % epsilon_max_sigma = epsilon_max_array;
% % % epsilon_min_sigma = epsilon_min_array;
% % % epsilon_max_sigma(30:29+N) = epsilon_max_array(30:29+N) - sigma_contraint;
% % % epsilon_max_sigma(95:94+N) = epsilon_max_array(95:94+N) - sigma_contraint;
% % % epsilon_min_sigma(30:29+N) = epsilon_min_array(30:29+N) + sigma_contraint;
% % % epsilon_min_sigma(45:44+N) = epsilon_min_array(45:44+N) + sigma_contraint;
% % % epsilon_min_sigma(100:99+N) = epsilon_min_array(100:99+N) + sigma_contraint;
% % % epsilon_min_sigma(160:159+N) = epsilon_min_array(160:159+N) + sigma_contraint;
% % % % epsilon_max_sigma(430:430+N-1) = epsilon_max_array(430:430+N-1) - sigma_contraint;
% % % plot(x_array, epsilon_max_sigma(1:Nsim) ,'--k');
% % % plot(x_array, epsilon_min_sigma(1:Nsim) ,'--k');
% % % grid on
% % % title(['Probabilidade: ',num2str(contraint_p*100),' %'])

%%
figure
hold on
plot(x_array, J_smpc_1)
plot(x_array, J_smpc_2)
plot(x_array, J_smpc_5)
plot(x_array, J_smpc_8)
plot(x_array, J_mpc,'k')
plot(x_array, J_rmpc,'g')
legend('SMPC - 1%','SMPC - 2%','SMPC - 5%','SMPC - 8%','MPC','RMPC')

cont_violations_smpc_1
cont_violations_smpc_2
cont_violations_smpc_5
cont_violations_smpc_8
cont_violations_mpc
cont_violations_rmpc
cont_violations_lqr

cont_cant_u_smpc_1
cont_cant_u_smpc_2
cont_cant_u_smpc_5
cont_cant_u_smpc_8
cont_cant_u_mpc
cont_cant_u_rmpc


%%
figure
hold on;
plot(x_array,v_smpc_1)
grid on
plot(x_array,v_smpc_2)
plot(x_array,v_smpc_5)
plot(x_array,v_smpc_8)
plot(x_array,v_mpc)
plot(ones(Nsim,1)*vmax,'--r')
plot(ones(Nsim,1)*vmin,'--r')
legend('SMPC - 1%','SMPC - 2%','SMPC - 5%','SMPC - 8%', 'MPC');
axis([0 150 -0.5 0.5])


%%
figure 
hold on;
plot(x_array, x_driver(1,1:Nsim))
plot(x_array, x_smpc_1(1,1:Nsim))
plot(x_array, x_smpc_2(1,1:Nsim))
plot(x_array, x_smpc_5(1,1:Nsim))
plot(x_array, x_smpc_8(1,1:Nsim))
plot(x_array, x_mpc(1,1:Nsim),'k')
plot(x_array, zeros(Nsim),'--','Color', uint8([100 100 100]));
ylabel('$e_y$ [m]','FontSize',12)
xlabel('x [m]','FontSize', 12)
h = legend('Condutor','SMPC - 1%','SMPC - 2%','SMPC - 5%','SMPC - 8%', 'MPC');
set(h,'Location','southeast');
plot(x_array, epsilon_max_array(1:Nsim) , '--','Color', uint8([255 100 100]));
plot(x_array, epsilon_min_array(1:Nsim) ,'--','Color', uint8([255 100 100]));
grid on
axis([0 150 -2.5 2.5])

%%
%%255,69,0
Nplot = 45;
start_plot = 25;
figure
hold on;
plot(x_array(start_plot:Nplot), x_smpc_2(1,start_plot:Nplot),'b')
% plot(x_array(start_plot:Nplot), x_smpc_2(1,start_plot:Nplot))
plot(x_array(start_plot:Nplot), x_smpc_5(1,start_plot:Nplot), 'Color', uint8([255 69 0]));
% plot(x_array(start_plot:Nplot), x_smpc_8(1,start_plot:Nplot))
plot(x_array(start_plot:Nplot), zeros(1,Nplot-start_plot+1),'--k')
ylabel('$e_y$ [m]','FontSize',12)
xlabel('x [m]','FontSize',12)
h = legend('SMPC - 2%','SMPC - 5%');
set(h,'Location','southeast');
% plot(x_array(start_plot:Nplot), epsilon_max_array(start_plot:Nplot) , 'Color', uint8([17 17 17]));
plot(x_array(start_plot:Nplot), epsilon_min_array(start_plot:Nplot) ,'Color', uint8([255 0 0]));
epsilon_max_sigma_2 = epsilon_max_array;
epsilon_min_sigma_2 = epsilon_min_array;
epsilon_max_sigma_5 = epsilon_max_array;
epsilon_min_sigma_5 = epsilon_min_array;
% epsilon_max_sigma_2(30:29+N) = epsilon_max_array(30:29+N) - sigma_contraint_2;
% epsilon_max_sigma_5(30:29+N) = epsilon_max_array(30:29+N) - sigma_contraint_5;
epsilon_min_sigma_2(30:29+N) = epsilon_min_array(30:29+N) + sigma_contraint_2;
epsilon_min_sigma_5(30:29+N) = epsilon_min_array(30:29+N) + sigma_contraint_5;
% plot(x_array(start_plot:Nplot), epsilon_max_sigma_2(start_plot:Nplot) ,'--b');
plot(x_array(start_plot:Nplot), epsilon_min_sigma_2(start_plot:Nplot) ,'--b');
% plot(x_array(start_plot:Nplot), epsilon_max_sigma_5(start_plot:Nplot) ,'--r');
plot(x_array(start_plot:Nplot), epsilon_min_sigma_5(start_plot:Nplot) , '--','Color', uint8([255 69 0]));
grid on
% xlim([26 32])
axis([26 32 -2.5 2])
% title('Comparação Compensação nas Restrições')

%%
rms(J_mpc(~isnan(J_mpc)))
rms(J_smpc_1(~isnan(J_smpc_1)))
rms(J_smpc_2(~isnan(J_smpc_2)))
rms(J_smpc_5(~isnan(J_smpc_5)))
rms(J_smpc_8(~isnan(J_smpc_8)))
