%% initialise
close all;
clear;
rng(0);
addpath('../', '../wcgainlbgrid');

%% define uncertain system
del_w0 = ureal('del_w0', 0);
del_xi = ureal('del_xi', 0);
w0 = 1 + 0.3 * del_w0;
K = 1;
xi = (1 / sqrt(2)) * (1 + 0.5 * del_xi);
P = tf(K * w0^2, [1, 2 * xi * w0, w0^2]);

%% define frequencies and objective
freq = [0.75, 1];
obj = @(del)(-sum(abs(squeeze(freqresp(usubs(P, 'del_w0', del(1), 'del_xi', del(2)), freq)))));

%% evaluate objective function over the parameter set
dd = linspace(-1, 1, 20);
dels = allcomb(dd, dd);
for kk = 1 : size(dels, 2)
	objs(kk) = obj(dels(:, kk));
end
% contour plot of the objective function
[X, Y] = meshgrid(dd, dd);
Z = -reshape(objs, size(X))';
contourf(X, Y, Z);
drawnow;

%% rund cube refinement search
[delopt, objopt, dd] = hypercube_interval_search(2, obj, [], [], 100);
% plot the points where the objective function was evaluated
hold on;
plot(dd(1, :), dd(2, :), 'rx', 'LineWidth', 3);
xlabel('del_{w0}');
ylabel('del_{xi}');
legend('obj. fun.', 'obj. fun. evaluated');
drawnow;
hold off;

%% local function definitions
function [delopt, objopt, dd] = hypercube_interval_search(n_dim, obj, con, delcon, n_eval_max)
	if nargin < 3
		con = [];
	end
	n_eval = 0;
	% initial grid
	grid = cell(n_dim, 1);
	grid(:) = {linspace(-1, 1, 2)};
	dels = allcomb(grid{:});
	dels = [mean(dels, 2), dels];
	n_dels = size(dels, 2);
	objval = nan(1, n_dels);
	ind = 1 : n_dels;
	% split cubes
	dd = dels;
	while n_eval <= n_eval_max
		for kk = ind
			delk = dels(:, kk);
			if isempty(con) || con(delk) <= 0
				objval(kk) = obj(delk);
				n_eval = n_eval + 1;
				if n_eval == n_eval_max
					break;
				end
			end
		end
		if all(isnan(objval))
			distances = diag((dels - delcon)' *(dels - delcon));
			[~, distsortind] = sort(distances);
			closestind = distsortind(1);
			dels(:, closestind) = delcon;
			continue;
		else
			[~, sortind] = sort(objval);
			minind = sortind(1);
			if minind == 1 % the minimum is in the center
				minind = sortind(2);
			end
			del_min = dels(:, minind);
			del_cen = dels(:, 1);
			objval_prev = objval;
			objval = inf(size(objval));
			what = mat2cell([del_min, del_cen]', 2, ones(1, n_dim));
			dels = allcomb(what{:});
			dels = [mean(dels, 2), dels];
			d1ind = find(all(dels == del_min));
			d2ind = find(all(dels == del_cen));
			objval(d1ind) = objval_prev(minind);
			objval(d2ind) = objval_prev(1);
			ind = 1 : n_dels;
			ind([d1ind, d2ind]) = [];
			dd = [dd, dels];
		end
	end
	delopt = dels(:, minind);
	objopt = objval_prev(minind);
end

function combs = allcomb(varargin)
	% calculate all combinations of a vector
	n = numel(varargin);
	ndgrid_output = cell(n, 1);
	[ndgrid_output{:}] = ndgrid(varargin{:});
	combs = cellfun(@(m)(reshape(m, [1, numel(m)])), ndgrid_output, 'UniformOutput', false);
	combs = cell2mat(combs);
end