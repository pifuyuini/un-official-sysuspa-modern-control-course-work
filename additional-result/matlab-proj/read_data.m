%% 查看地面扰动数据文件里有些啥
file_gnd = 'rec1_142_without_control_OK.mf4';
m_gnd    = mdf(file_gnd);

% 看看里面有哪些通道
names_gnd = m_gnd.ChannelNames;   % 这是个嵌套 cell
chan_gnd = names_gnd{1,1};        % 把真正的名字列表拎出来，方便查看
disp(chan_gnd)                    % 在命令行看看通道名

%% 查看噪声数据文件里有些啥
file_noise = 'rec1_005_ADC_noise.mf4';
m_noise    = mdf(file_noise);

names_noise = m_noise.ChannelNames;
chan_noise  = names_noise{1,1};
disp(chan_noise)   % 看看这个文件里都记录了什么

%% 十分难绷地读取数据（通常只需直接运行这个就好了）
% 找到地面通道的索引
idx_HL_gnd = find(strcmp(chan_gnd, 'Model Root/H_L_ground2/In1'));
idx_VL_gnd = find(strcmp(chan_gnd, 'Model Root/V_L_ground/In1'));

% 安全检查一下有没有找到
if isempty(idx_HL_gnd) || isempty(idx_VL_gnd)
    error('找不到 H_L_ground 或 V_L_ground 这些通道名，检查一下文件或名字是否拼错');
end

% 读出地面输入信号（向量）
H_L_ground = read(m_gnd, 1, chan_gnd{idx_HL_gnd}, 'OutputFormat','vector');
V_L_ground = read(m_gnd, 1, chan_gnd{idx_VL_gnd}, 'OutputFormat','vector');

% 找噪声通道
idx_noise = find(strcmp(chan_noise, 'Model Root/OUTPUT NOISE/In1'));
if isempty(idx_noise)
    error('找不到 OUTPUT NOISE 通道，请检查文件和通道名');
end

% 读取原始噪声
noise_raw = read(m_noise, 1, chan_noise{idx_noise}, 'OutputFormat','vector');

% 你原来是除以 1500 做标定，这里保持不变
vel1 = noise_raw / 1500;

%% 构建用于仿真的时间序列
clc;
fs = 10000;     % 10 kHz
Ts = 1/fs;
N  = length(H_L_ground);   % 数据点数（以水平地面信号为参考）

t  = (0:N-1)' * Ts;         % 列向量

H_L_ts = timeseries(H_L_ground, t);
V_L_ts = timeseries(V_L_ground, t);

t2 = (0:length(vel1)-1)' * Ts;
NOISE_ts = timeseries(vel1, t2);     % 假设噪声采样率一样

%% 如果要导出数据

% writematrix([NOISE_ts.time, NOISE_ts.data], 'noise_ts.csv');
% writematrix([H_L_ts.time, H_L_ts.data, V_L_ts.data], 'ground_ts.csv');


