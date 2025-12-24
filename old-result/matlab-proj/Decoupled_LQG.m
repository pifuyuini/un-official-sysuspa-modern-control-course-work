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

%%

% 开环输出典型量级
Vmag = 1e-4;    % y1,y3 竖直
Hmag = 1e-5;    % y2,y4 水平

%% KF
sigma_wx = 1e-6;   
sigma_wz = 1e-6;

Qw = diag([sigma_wx^2, sigma_wz^2]);   % 2x2

sigma_v = 1e-9;
Rv = (sigma_v^2) * eye(3);             % 4x4
beta_Qn = 1;
beta_Rn = 1;

W = Qw * beta_Qn;
V = Rv * beta_Rn;

L = lqe(sys_zv.A, Bw, sys_zv.C, W, V); % Kalman Gain

%% LQR
% 把输出按“典型量级”归一化：y_norm = Sy * y
Sy = diag([1/Vmag, 1/Vmag, 1/Hmag]);   % 4x4

% 对归一化后的输出 z = Sy*y，统一罚 1
Wy = Sy.' * Sy;        % = diag([1/Vmag^2, 1/Hmag^2, ...])

% 状态权重 Q：从输出权重回推
Q_base = (sys_zv.C).' * Wy * sys_zv.C;    % 8x8

R_base = eye(3);    % 三个控制力一视同仁，后面整体缩放

alpha_Q = 1;        % 输出惩罚刻度
alpha_R = 0.1;        % 力惩罚刻度（越大 => 力越小、越温柔）

Q = alpha_Q * Q_base;
R = alpha_R * R_base;

K = lqr(sys_zv.A, sys_zv.B, Q, R);

%% LQG Synthsis
Alqg = sys_zv.A - sys_zv.B*K - L*sys_zv.C;
Blqg = L;
Clqg = -K;
Dlqg = zeros(3, 3);