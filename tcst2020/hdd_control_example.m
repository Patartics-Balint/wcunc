% initialise
clear;
close all;
rng(0);
load('HDDdata.mat');
freq = reshape(wexp, [1, numel(wexp)]);
addpath('../', '../wcgainlbgrid');

%% fit nominal model
% remove the double integrator
resp = RespData .* (1i * freq').^2;
nom_resp = mean(resp, 2);
bw_fit = 7; % the nominal model shall be accurate up to this frequency
% construct low frequency weighting
weight_nom0 = [logspace(2, 0, numel(freq(freq <= bw_fit))), zeros(size(freq(freq > bw_fit)))];
weight_nom0 = frd(weight_nom0, freq);
hdd_nom0 = fitfrd(frd(nom_resp, freq), 2, 0, weight_nom0);
hdd_nom0_resp = squeeze(freqresp(hdd_nom0, freq));
wnom0_resp = squeeze(weight_nom0.ResponseData);
% plot amplitudes
loglog(freq, abs(resp), 'b', freq, abs(nom_resp), 'g', freq, abs(wnom0_resp), 'r',...
	freq, abs(hdd_nom0_resp), 'k--', 'LineWidth', 2);
title('Frequency responses without the double integrator');
xlabel('frequency [rad/s]');
ylabel('magnitude');
axis tight;
drawnow;
% add the double intergrator back in
hdd_nom = hdd_nom0 * tf(1, conv([1, 0], [1, 0]));

%% Bode plots of the nominal system and measured data
hdd_nom_resp = frd(hdd_nom, freq);
hdd_nom_gain = abs(squeeze(hdd_nom_resp.ResponseData));
hdd_nom_phase = angle(squeeze(hdd_nom_resp.ResponseData));
nom_resp = frd(hdd_nom, freq);
nom_gain = abs(squeeze(nom_resp.ResponseData));
nom_phase = angle(squeeze(nom_resp.ResponseData));
% find resonance frequencies
ii = find(freq > 2); % exclude low frequencies because of the double integrator
gain_exp = max(abs(RespData(ii, :)), [], 2);
[peak_exp, pf_exp] = findpeaks(gain_exp, freq(ii),...
	'MinPeakHeight', 7, 'MinPeakDistance', 2);
pf_ind = find(arrayfun(@(f)(any(f == pf_exp)), freq));

% plot gain
figure;
subplot(2, 1, 1);
loglog(freq, abs(RespData), 'b');
hold on;
	loglog(freq, hdd_nom_gain, 'r', freq, nom_gain, 'k', pf_exp, peak_exp, 'rx', 'LineWidth', 2);
hold off;
ylabel('magnitude');
xlabel('frequency [rad/s]');
axis tight;
% plot phase
subplot(2, 1, 2);
semilogx(freq, rad2deg(angle(RespData)), 'b');
hold on;
	semilogx(freq, rad2deg(hdd_nom_phase), 'r', freq, rad2deg(nom_phase), 'k', 'LineWidth', 2);
hold off;
ylabel('phase [^\circ]');
xlabel('frequency [rad/s]');
axis tight;
drawnow;

%% plot error and uncertainty weight
rel_err = (RespData - squeeze(hdd_nom_resp.ResponseData)) ./...
	squeeze(hdd_nom_resp.ResponseData);
der = tf([1, 0], 1)^2; % double integrator
der_resp = freqresp(der, freq);
der_resp = abs(squeeze(der_resp))';
rel_err_ub = max(abs(rel_err)') ./ der_resp;
% plot the error and the upper bound
figure;
loglog(freq, abs(rel_err), 'b', freq, abs(rel_err_ub), 'k');
data = frd(rel_err_ub, freq);
data = genphase(data);
% find frequency weighting to make the fit more accurate
freqweight = 1 ./ rel_err_ub;
m = max(freqweight);
freqweight(pf_ind([1, 2, 3, 5])) = m;
freqweight(pf_ind(4)) = 0.08 * m;
freqweight(freq < 1.1) = 0;
freqweight = frd(freqweight, freq);
% fit dynamic uncertainty weight
nxWu = 25;
Wu0 = fitfrd(data, nxWu, 2, freqweight, 0);
if max(real(eig(Wu0))) > 0 % reflect the right half plane poles on the imaginary axis 
	try
		Wu0 = spectralfact(Wu0, []);
	catch er
		warning(er.message);
		[z, p, k] = zpkdata(Wu0);
		p = p{1};
		p(real(p) > 0) = -p(real(p) > 0);
		Wu0 = zpk(z, p, k);
	end
end
% add double integrator back in
Wu = ss(tf(Wu0) * der);
% plot the gain of the uncertainty weight
Wu_resp = frd(Wu, freq);
hold on;
	loglog(freq, abs(squeeze(Wu_resp.ResponseData)), 'r', 'LineWidth', 2);
	vline(freq(pf_ind), 'c');
hold off;
title('Uncertainty weight vs. experimental data');
axis tight;
drawnow;

%% construct uncertain model
del_dyn = ultidyn('del_dyn', [1, 1]); % dynamic uncertainty block
del_theta = ureal('del_theta', 0); % uncertainty of the rotational inertia
hdd_unc = hdd_nom / (1 + 0.1 * del_theta) * (1 + Wu * del_dyn);
figure;
sigma(hdd_unc, 'b', hdd_unc.NominalValue, 'g', freq);
drawnow;

%% analyse sensitivity and complementary sensitivity
freq = logspace(log10(freq(1)) - 1, log(freq(end)) - 1, 4 * numel(freq))';
beta = 2;
wB = 0.45;
K = tf([beta, wB], [1, beta * wB]) / abs(freqresp(hdd_nom, wB));
L = hdd_unc * K; % open loop
S = feedback(1, L); % closed-loop sensitivity
% check robust stability of the closed-loop
rs = robstab(S);
if rs.UpperBound < 1 || rs.LowerBound < 1
	error('The closed loop is not robustly stable');
end

%% calculate worst-case uncertainties
gain = @(sys)(calc_gain(sys, freq));
wcgT = wcgain(S);
lb = wcgainlbgrid(S, freq);
[wcpeaks, wcfreq] = findpeaks(lb, freq, 'MinPeakHeight', 1.23, 'MinPeakDistance', 2);
wcfreq = reshape(wcfreq, [1, numel(wcfreq)]);
[~, imax] = max(wcpeaks);
wcfreq(imax) = wcgT.CriticalFrequency;
wcpeaks(imax) = wcgT.LowerBound;
% worst-case uncertainty at the peak
[wcuT1, Twc1] = wcunc(S, wcfreq(imax));
[wcuT2, Twc2, info2] = wcunc(S, wcfreq);

%% plot refquency domain results
figure;
sigma(S.NominalValue, 'g', frd(lb, freq), 'k', Twc1, 'b', Twc2, 'r--', {0.2, 50});
vline(info2.freq, 'r-.');
vline(wcgT.CriticalFrequency, 'b:');
legend('nom', 'LB', 'single-peak wcu', 'multi-peak wcu', 'Location', 'NorthWest');
xlim([0.2, 50]);
ylim([-5, inf]);
drawnow;

%% analyse time domain simulations
figure;
t = linspace(0, 600, 2e5);
msfreq = wcfreq;
msamp = 0.01 * ones(size(msfreq));
msphase = linspace(0, 2 * pi, numel(msfreq));
sins = arrayfun(@(a, w, p)(a * sin(w * t + p)), msamp, msfreq, msphase, 'UniformOutput', false);
u = sum(cell2mat(sins')) + 1;
sim_nom = lsim(S.NominalValue, u, t);
sim_wc1 = lsim(Twc1, u, t);
sim_wc2 = lsim(Twc2, u, t);
plot(t, sim_nom, 'g', t, sim_wc1, 'b', t, sim_wc2, 'r');
legend('nom', 'single-peak wcu', 'multi-peak wcu');
axis('tight');
xlim([0, 40]);
drawnow;

% analyse stead-state response
anal_ind = t > 300;
fprintf('\tmax\tstd. dev.\n');
fprintf('nom.\t%.4f\t%.4f\n', max(abs(sim_nom(anal_ind))), std(sim_nom(anal_ind)));
fprintf('sp.\t%.4f\t%.4f\n', max(abs(sim_wc1(anal_ind))), std(sim_wc1(anal_ind)));
fprintf('mp.\t%.4f\t%.4f\n', max(abs(sim_wc2(anal_ind))), std(sim_wc2(anal_ind)));
samp = usample(S, 100);
std_dev = zeros(nmodels(samp), 1);
amp = zeros(nmodels(samp), 1);
for kk = 1 : nmodels(samp)
	sim_samp = lsim(samp(:, :, kk), u, t);
	std_dev(kk) = std(sim_samp(anal_ind));
	amp(kk) = max(abs(sim_samp(anal_ind)));
end
[~, ind_amp_max] = max(amp);
[~, ind_dev_max] = max(std_dev);
fprintf('mamp.\t%.4f\t%.4f\n', amp(ind_amp_max), std_dev(ind_amp_max));
fprintf('mdev.\t%.4f\t%.4f\n', amp(ind_dev_max), std_dev(ind_dev_max));

%% local functions
function gain = calc_gain(sys, freq)
	gain = sigma(sys, freq);
	gain = gain(1, :)';
end

function val = eval_obj(usys, wcfreq)
	fgrid = [0, logspace(-5, 4, 100)];
	wcu = wcunc(usys, wcfreq);
	sys = usubs(usys, wcu);
	g = calc_gain(sys, fgrid);
	val = trapz(fgrid, log(g));
end