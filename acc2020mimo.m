clear;
close all;

%% create a robustly stable uncertain system
load('acc2020mimo.mat', 'P');
% pick frequency points
freq = [0.1, 2, 22, 200];

%% classical construction
opt = wcOptions('VaryFrequency', 'on');
[wcg, wcu_bp, wcginfo] = wcgain(P, opt);
for delname = fieldnames(wcu_bp)'
	delname = delname{1};
	del = wcu_bp.(delname);
	[u, ~, v] = svd(freqresp(del, wcg.CriticalFrequency));
	u = u(:, 1);
	vH = v(:, 1)';
	u_ss = [];
	for kk = 1 : numel(u)
		u_ss = [u_ss; cnum2sys(u(kk), wcg.CriticalFrequency)];
	end
	vH_ss = [];
	for kk = 1 : numel(vH)
		vH_ss = [vH_ss, cnum2sys(vH(kk), wcg.CriticalFrequency)];
	end
	del_ss = u_ss * vH_ss;
	wcu_s.(delname) = del_ss;
end

%% calculate worst-case uncertainty
[wcu_m, ~, info] = wcunc(P, freq);

%% plot result
[~, ~, pointsinfo] = wcgain(P, freq);
sigma(P.NominalValue, 'g', usubs(P, wcu_m), 'r', usubs(P, wcu_s), 'b--',...
	frd(pointsinfo.Bounds(:, 1), pointsinfo.Frequency), 'kx',...
	frd(wcginfo.Bounds(:, 1), wcginfo.Frequency), 'k:');
legend('nominal', 'wcu_m', 'wcu_s', 'points', 'LB');
figure;
step(P.NominalValue, 'g', usubs(P, wcu_m), 'r', usubs(P, wcu_s), 'b--', 200);
legend('nominal', 'wcu_m', 'wcu_s');

%% Recompute wcu_s with BNP
wcu_s2 = wcunc(P, wcg.CriticalFrequency);
del_s2 = blkdiag(wcu_s2.del1, wcu_s2.del2, wcu_s2.del3); % has 6 states