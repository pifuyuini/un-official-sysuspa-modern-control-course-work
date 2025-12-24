%% ===== 2) Output transform: y -> z (average/difference) =====
clc;
% y = [y1; y2; y3; y4]
% z1 = yv      = (y1 + y3)/2
% z2 = yv_diff = (y1 - y3)/2
% z3 = yh      = (y2 - y4)/2   (since your C makes y4 ≈ -y2, this "fixes" sign)
% z4 = yh_sum  = (y2 + y4)/2   (should be ~0, just for sanity check)
Ty = [ ...
 0.5  0    0.5  0
 0.5  0   -0.5  0
 0    0.5  0   -0.5];

% ===== 3) Input transform: u -> v (common/differential for u1,u2) =====
% v = [uSigma; uDelta; u3]
% u1 = uSigma + uDelta
% u2 = uSigma - uDelta
% u3 = u3
Tu = [ 1  1  0
       1 -1  0
       0  0  1 ];

sys_zv = ss(A, Bu*Tu, Ty*C, zeros(3,3)); % (z,v) plant, "virtual" channels

Gv  = minreal(sys_zv(1,1));  
Gvd = minreal(sys_zv(2,2));  
Gh  = minreal(sys_zv(3,3)); 


% %% SVD Decouple
% clc;
% 
% w = 2 * 2 * pi;
% Gjw = freqresp(ss(A, Bu, C, Du), w);
% Gjw = Gjw(:,:,1);   % 4x3
% G0 = dcgain(ss(A, Bu, C, Du));
% 
% [Usvd,Ssvd,Vsvd] = svd(G0);
% 
% V_in  = Vsvd;          % 输入变换（3x3）
% U_out = Usvd(:,1:3);   % 只保留前三个输出方向（4x3）
% 
% B_new = Bu * V_in;
% D_new = Du * V_in;
% 
% C_new = U_out' * C;
% D_new = U_out' * D_new;
% 
% sys_decoupled = ss(A, B_new, C_new, D_new);
% 
% Gv  = sys_decoupled(1,1);  
% Gvd = sys_decoupled(2,2);  
% Gh  = sys_decoupled(3,3); 
% 
% Ty = U_out';
% Tu = V_in;

%% KF

% 开环输出典型量级
Vmag = 1e-4;    % y1,y3 竖直
Hmag = 1e-5;    % y2,y4 水平

sigma_wx = 1e-6;   
sigma_wz = 1e-6;

Qw = eye(size(Gv.A)) * sigma_wz^2;

sigma_v = 1e-9;
Rv = (sigma_v^2) * eye(1);
beta_Qn = 1;
beta_Rn = 1;

W = Qw * beta_Qn;
V = Rv * beta_Rn;

L = lqe(Gv.A, eye(size(Gv.A)), Gv.C, W, V); % Kalman Gain

%% LQR
% 把输出按“典型量级”归一化：y_norm = Sy * y
Sy = diag([1/Vmag]);   

% 对归一化后的输出 z = Sy*y，统一罚 1
Wy = Sy.' * Sy;        % = diag([1/Vmag^2, 1/Hmag^2, ...])

% 状态权重 Q：从输出权重回推
Q_base = (Gv.C).' * Wy * Gv.C;    % 8x8

R_base = eye(1);    % 三个控制力一视同仁，后面整体缩放

alpha_Q = 1;        % 输出惩罚刻度
alpha_R = 0.1;        % 力惩罚刻度（越大 => 力越小、越温柔）

Q = alpha_Q * Q_base;
R = alpha_R * R_base;

K = lqr(Gv.A, Gv.B, Q, R);

%% LQG Synthsis
Alqg = Gv.A - Gv.B*K - L*Gv.C;
Blqg = L;
Clqg = -K;
Dlqg = zeros(1, 1);

%% LQG 2 Gvd
% KF

Qw2 = eye(size(Gvd.A)) * sigma_wz^2;
Rv2 = (sigma_v^2) * eye(1);
beta_Qn2 = 1;
beta_Rn2 = 1;

W2 = Qw2 * beta_Qn2;
V2 = Rv2 * beta_Rn2;

L2 = lqe(Gvd.A, eye(size(Gvd.A)), Gvd.C, W2, V2); % Kalman Gain

% LQR

% 状态权重 Q：从输出权重回推
Q_base2 = (Gvd.C).' * Wy * Gvd.C;    % 8x8

R_base2 = eye(1);    % 三个控制力一视同仁，后面整体缩放

alpha_Q2 = 0.0005;        % 输出惩罚刻度
alpha_R2 = 0.1;        % 力惩罚刻度（越大 => 力越小、越温柔）

Q2 = alpha_Q2 * Q_base2;
R2 = alpha_R2 * R_base2;

K2 = lqr(Gvd.A, Gvd.B, Q2, R2);

% LQG Synthsis
Alqg2 = Gvd.A - Gvd.B*K2 - L2*Gvd.C;
Blqg2 = L2;
Clqg2 = -K2;
Dlqg2 = zeros(1, 1);

%% LQG3 Gh
% KF

Qw3 = eye(size(Gh.A)) * sigma_wx^2;

beta_Qn3 = 1;
beta_Rn3 = 1;

W3 = Qw3 * beta_Qn3;
V3 = Rv * beta_Rn3;

L3 = lqe(Gh.A, eye(size(Gh.A)), Gh.C, W3, V3); % Kalman Gain

% LQR
% 把输出按“典型量级”归一化：y_norm = Sy * y
Sy3 = diag([1/Hmag]);   

% 对归一化后的输出 z = Sy*y，统一罚 1
Wy3 = Sy3.' * Sy3;        % = diag([1/Vmag^2, 1/Hmag^2, ...])

% 状态权重 Q：从输出权重回推
Q_base3 = (Gh.C).' * Wy3 * Gh.C;    % 8x8

R_base3 = eye(1);    % 三个控制力一视同仁，后面整体缩放

alpha_Q3 = 1;        % 输出惩罚刻度
alpha_R3 = 0.1;        % 力惩罚刻度（越大 => 力越小、越温柔）

Q3 = alpha_Q3 * Q_base3;
R3 = alpha_R3 * R_base3;

K3 = lqr(Gh.A, Gh.B, Q3, R3);

% LQG Synthsis
Alqg3 = Gh.A - Gh.B*K3 - L3*Gh.C;
Blqg3 = L3;
Clqg3 = -K3;
Dlqg3 = zeros(1, 1);
