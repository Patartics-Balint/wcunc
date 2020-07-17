%% initialise
close all;
clear;
addpath('../', '../wcgainlbgrid');

%% define uncertain system
F1 = tf([2, 0], conv([1, 1], [1, 1]));
F2 = tf([2 * 1 / 100, 0], conv([1/100, 1], [1/100, 1]));
alpha = 2;
delr = ureal('delr', 0);
Wdyn = 0.1 * tf([1 / 7, 1], [1 / 2 / 7, 1]);
usys = (alpha * F1 * (0.2 + delr) + F2 * (0.9 - delr)) * (1 + Wdyn * ultidyn('deld', [1, 1]));

%% specify frequencies and find the worst-case uncertainty at the specified points
pfreq = [1, 100];
% save lower bound at these frequencies for plotting
wcglb = wcgainlbgrid(usys, pfreq);
points = frd(wcglb, pfreq);

%% find worst-case uncertainty
[wcu1, wcsys1] = wcunc(usys, pfreq(1)); % the single peak uncertainty is also calculated using our algorithm
[wcu2, wcsys2] = wcunc(usys, pfreq);

%% print results
msv1 = sigma(wcsys1, pfreq);
msv2 = sigma(wcsys2, pfreq);
fprintf('frequencies \t %.0f \t %.0f \t sum\n', pfreq(1), pfreq(2));
fprintf('wcg. lb. \t %.2f \t %.2f \t %.2f\n', wcglb(1), wcglb(2), sum(wcglb));
fprintf('singular values\n');
fprintf('single-peak \t %.2f \t %.2f \t %.2f\n', msv1(1), msv1(2), sum(msv1));
fprintf('multi-peak \t %.2f \t %.2f \t %.2f\n', msv2(1), msv2(2), sum(msv2));

%% plot results
fr = unique([logspace(-1, 3, 200), pfreq]);
lb = wcgainlbgrid(usys, fr);
lb = frd(lb, fr);
sigma(usys.NominalValue, 'g', lb, 'k', wcsys1, 'b', wcsys2, 'r--', fr);
hold on;
	sigma(points, 'kx');
hold off;
legend('nominal', 'LB', 'single-peak', 'multi-peak', 'points', 'Location', 'southwest');
drawnow;

%% plot objective function
obj = @(del)(sum(wcgainlbgrid(usubs(usys, 'delr', del), pfreq)));
ddel = linspace(-1, 1, 101);
for kk = 1 : numel(ddel)
	J(kk) = obj(ddel(kk));
end
figure;
plot(ddel, J, 'r', [-1, 1], [1, 1] * sum(wcglb), 'b');
xlabel('delr');
ylabel('objective function');
legend('objective function', 'upper bound');