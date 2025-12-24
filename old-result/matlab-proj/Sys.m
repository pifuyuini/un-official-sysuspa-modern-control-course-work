
clear;clc;
load('Original+model.mat');
fs = 10000;
Ts = 1/fs;

%% CT System
A = sys_3dof_w2v.A;
B = sys_3dof_w2v.B;
C = sys_3dof_w2v.C;
D = sys_3dof_w2v.D;
Bu = B(:, 1:3);
Bw = B(:, 4:5);
Du = D(:, 1:3);
Dw = D(:, 4:5);

%% CT Grd Dis
Kg = 150/3/1e6;

zg = [-1.3*2*pi*(1+1i), -1.3*2*pi*(1-1i)];
dzg = [-1.3*2*pi*(1+1i), -1.3*2*pi*(1-1i), 0];
pg = [-0.03*2*pi + 1i*(0.18*2*pi), -0.03*2*pi - 1i*(0.18*2*pi), ...
      -20*2*pi, -20*2*pi];

Gw = zpk(zg, pg, Kg);
[Ag, Bg, Cg, Dg] = ssdata(Gw);
Gdw = zpk(dzg, pg, Kg);
[Adg, Bdg, Cdg, Ddg] = ssdata(Gdw);

%% CT Sensor Noise
Ks = 3000*2*pi*db2mag(-200);

zs = [-2, -2, -2, -2];
ps = [-0.1, -0.1, -0.1, -0.1, -3000*2*pi];

Gns = zpk(zs, ps, Ks);
[Ans, Bns, Cns, Dns] = ssdata(Gns);
