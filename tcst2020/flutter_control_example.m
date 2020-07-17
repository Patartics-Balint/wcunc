%% initialise
clear;
close all;
rng(0);
load('flutter_control.mat');
addpath('../', '../wcgainlbgrid');

%% add delay to the uncertain aircraft model
delay = pade(tf(1, 'IOdelay', 0.015), 4);
plant = usys * delay;
plant.y = usys.y;
plant.u = usys.u;

%% construct MIMO flutter controller
yblend_sym = [-2, 1, 1];
ublend_sym = [1; 1];
yblend_asym = [0, 1, -1];
ublend_asym = [1; -1];
Kf = [ublend_sym, ublend_asym] * blkdiag(K_sym, K_asym) * [yblend_sym; yblend_asym];
Kf.y = {'Aileron4L', 'Aileron4R'};
Kf.u = {'q', 'IMU_L6_ry', 'IMU_R6_ry'};

%% closed loop
Lu = [ublend_sym, ublend_asym];
Ly = [yblend_sym; yblend_asym];
K = blkdiag(K_sym, K_asym);
Ldes = Ly * plant(3 : end, :) * Lu * K;
S = feedback(eye(2), Ldes);
Td = feedback(plant, Kf, 1:2, 3:5);
Td.u = Kf.y;
Tu = feedback(Lu * K, Ly * plant(3 : end, :));

%% worst-case uncertainty
wcfreq = [61.0876, 26.5919];
[wcu1, wcTT1] = wcunc(S, wcfreq(1));
[wcu2, wcTT2, info2] = wcunc(S, wcfreq);
wcTd1 = usubs(Td, wcu1);
wcTd2 = usubs(Td, wcu2);
wcS1 = usubs(S, wcu1);
wcS2 = usubs(S, wcu2);

%% plot closed loop transfer function
fp = {1e0, 1e3};
defl = deg2rad(1);
% singular values of the sensitivity function
figure;
sigma(S.NominalValue, 'g',	wcTT1, 'k', wcTT2, 'r--', fp);
legend('nominal', 'single-peak', 'multi-peak', 'Location', 'SouthWest');
vline(info2.freq, 'r:');
% magnitude of the transfer function from input disturbance to output
figure;
bodemag(Td.NominalValue, 'g',	wcTd1, 'k', wcTd2, 'r--', fp);
legend('nominal', 'single-peak', 'multi-peak', 'Location', 'NorthWest');
vline(info2.freq, 'r:');
% step response from all input to all output
figure;
step(Td.NominalValue * defl, 'g', wcTd1 * defl, 'k', wcTd2 * defl, 'r--', 0.5);
legend('nominal', 'single-peak', 'multi-peak', 'Location', 'SouthEast');

%% plot step response from aileron to relative acceleration
figure;
step(Td(1, 1) * defl, 'y', Td(1, 1).NominalValue * defl, 'g', wcTd1(1, 1) * defl, 'k', wcTd2(1, 1) * defl, 'r', 0.5);
legend('random', 'nominal', 'single-peak', 'multi-peak', 'Location', 'SouthEast');