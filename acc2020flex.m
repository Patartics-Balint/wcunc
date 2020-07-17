clear;
close all;

%% create vibration control example
w1 = 1;
G1 = tf(w1^2, [1 2*0.4*w1 w1^2]); 
w2 = 4;
G2 = tf(w2^2, [1 2*0.05*w2  w2^2]);
del1 = ultidyn('del1', [1, 1]);
del2 = ultidyn('del2', [1, 1]);
Gflex = G1 * (1 + 0.3 * del1) * G2 * (1 + 0.3 * del2);
P = feedback(Gflex, 1);

%% construct band bass uncertainty
opt = wcOptions;
opt.VaryFrequency = 'on';
[wcg, wcu_bp, wcginfo] = wcgain(P, logspace(-1, 1, 100));

%% construct all pass uncertainty
del1_crit = freqresp(wcu_bp.del1, wcg.CriticalFrequency);
del2_crit = freqresp(wcu_bp.del2, wcg.CriticalFrequency);
wcu_s.del1 = cnum2sys(del1_crit, wcg.CriticalFrequency);
wcu_s.del2 = cnum2sys(del2_crit, wcg.CriticalFrequency);

%% interpolate the peaks of the upper bound
[~, sortidx] = sort(wcginfo.Bounds(:, 1), 1, 'descend');
maxidx = sortidx([1, 6]);
[wcu_m, ~, info] = wcunc(P, [wcginfo.Frequency(maxidx)]);

%% plot gain and step response
figure;
sigma(P.NominalValue, 'g', usubs(P, wcu_s), 'c', usubs(P, wcu_m), 'r',...
			frd(wcginfo.Bounds(:, 1), wcginfo.Frequency), 'k',...
			frd(wcginfo.Bounds(maxidx, 1), wcginfo.Frequency(maxidx)), 'xk');
xlim([1e-1, 1e1]);
legend('nominal', 'wcu_s', 'wcu_m', 'LB', 'points');
figure;
step(P.NominalValue, 'g', usubs(P, wcu_s), 'c', usubs(P, wcu_m), 'r');
legend('nominal', 'wcu_s', 'wcu_m');