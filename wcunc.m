function [wcu, wcsys, info] = wcunc(usys, freq)
	% Calculate worst-case uncertainty which maximises the gain of an
	% uncertain system at specified frequnecies.
	%
	% Usage
	%   [wcu, wcsys, info] = WCUNC(usys, freq)
	% Inputs:
	%   usys: uss
	%   freq: Vector of frequencies where the gain of usys is to be maximised.
	% Outputs:
	%   wcu: Worst-case uncertainty. The gain of usubs(usys, wcu) is maximal
	%   at the frequencies in freq.
	%   wcsys: usubs(usys, wcu)
	%   info: Structure containing information about the interpolation for
	%   each uncertainty block of the dynamic uncertainty.
	%
	% See also wcgain, bnpinterp
	
	[dnames, rnames, nr] = process_system(usys);
	[susys, cfreq, info] = scale_for_robust_stability(usys);
	if nargin < 2
		freq = pick_freq_grid(susys, info);
	end
	[freq, info] = process_freq(freq, susys, info);

	if ~isempty(rnames) % there is parametric uncertainty in the system
		if numel(freq) == 1 && isempty(cfreq)
			[~, wcu_full] = wcgainlbgrid(usys, freq);
			for kk = 1 : numel(rnames)
				wcu.(rnames{kk}) = wcu_full.(rnames{kk});
			end
		else
			[obj, con, delrdestab] = pick_obj_and_con(susys, freq, cfreq, rnames, dnames);
			n_eval_max = 100;
			delropt = hypercube_interval_search(nr, obj, con, delrdestab, n_eval_max);
			if all(isnan(delropt))
				delropt = delrdestab;
			end
			for kk = 1 : numel(rnames)
				wcu.(rnames{kk}) = delropt(kk);
			end
		end
	end
	
	if ~isempty(dnames) % there is dynamic uncertainty in the system
		if exist('wcu')
			susys_dyn = usubs(susys, wcu);
		else
			susys_dyn = susys;
		end
		if ~isempty(cfreq) % draw destabilising sample
			[~, ~, robinfo] = robstab(susys_dyn, cfreq);
		end
	  % calculate worst-case dynamic uncertainty block by block.
		[~, wcpert] = wcgainlbgrid(susys_dyn, freq);
		for kblk = 1 : numel(dnames)
			blkname = dnames{kblk};
			delsamp = [];
			lvec = [];
			rvec = [];
			for kfr = 1 : numel(freq)
				lvec(:, kfr) = wcpert(kfr).(blkname){1};
				rvec(:, kfr) = wcpert(kfr).(blkname){2}';
			end
			if ~isempty(cfreq)
				stabsamp = freqresp(robinfo.WorstPerturbation.(blkname), cfreq);
				[l, s, r] = svd(stabsamp);
				l = l(:, 1) * sqrt(s(1));
				r = r(:, 1) * sqrt(s(1));
				if ~isempty(freq)
					lvec(:, end + 1) = l;
					rvec(:, end + 1) = r;
				else
					lvec = l;
					rvec = r;
				end
			end
			[wcublk, infoblk] = bnpinterp({lvec, rvec}, [freq, cfreq]);
			info.(blkname) = infoblk;
			wcu.(blkname) = wcublk;			
		end
		check_result(susys_dyn, wcu, freq);
	end
	% scale uncertainty block
	if ~isempty(cfreq)
		for name = [rnames; dnames]'
			name = name{1};
			wcu.(name) = info.uscale * wcu.(name);
		end
	end
	% save uncertainty block as a dynamic system
	ublk = [];
	for name = [rnames; dnames]'
		name = name{1};
		ublk = blkdiag(ublk, wcu.(name));
	end
	info.ublk = ublk;
	% substitute worst-case uncertainty
	wcsys = usubs(usys, wcu);
end
	
function freq = pick_freq_grid(susys, info)	
	freqs = [0, logspace(-6, 6, 150)]';
	if info.uscale == 1
		wcg = wcgain(susys);
		if isfinite(wcg.CriticalFrequency)
			freqs = unique([wcg.CriticalFrequency; freqs]);
		end
	end
	if info.uscale < 1
		[~, ~, robinfo] = robstab(susys, freqs);
		freqs(robinfo.Bounds(:, 1) <= 1) = [];	
	end
	wcglb = wcgainlbgrid(susys, freqs);
	sig = sigma(susys.NominalValue, freqs);
	sig = sig(1, :);
	sig = reshape(sig, [numel(sig), 1]);
	freqs = [-flip(freqs); freqs];
	wcglb = [flip(wcglb); wcglb];
	sig = [flip(sig); sig];
	[freqs, finds] = unique(freqs);
	wcglb = wcglb(finds);
	sig = sig(finds);
	dif = wcglb - sig;
	if info.uscale == 1
		[~, maxind] = max(wcglb);
		freq = freqs(maxind);
		[~, maxind] = max(dif);
		freq = [freq; freqs(maxind)];
	else
		freq = [];
	end
	range = max(dif) - min(dif);
	if info.uscale < 1
		range = min([range, 100]);
	end
	[pks, pfreq] = findpeaks(dif, freqs, 'Threshold', 0.01 * range);
	if isempty(pfreq)
		[pks, pfreq] = findpeaks(dif, freqs, 'NPeaks', 1);
	end
	if numel(pks) < 4
		freq = [freq; pfreq];
	else
		[~, sortinds] = sort(pks, 'descend');
		freq = [freq; pfreq(sortinds(1 : 4))];
	end
	freq = unique(abs(freq));
end

function check_result(usys, wcu, freq)
	% Check the accuracy of the worst-case dynamic uncertainty
	% construction.
	lb1 = wcgainlbgrid(usys, freq);
	wcsys = usubs(usys, wcu);
	resp = freqresp(wcsys, freq);
	for kk = 1 : size(resp, 3)
		lb2(kk, 1) = norm(resp(:, :, kk), 2);
	end
	tol = 1e-3;
	den = lb1;
	den(den < 0.1 * tol) = 1;
	err = abs(lb2 - lb1) ./ den;
	if max(err) > tol
		error('The gain of the system does not match the desired lower bound at the specified frequnecies.');
	end
end

function [obj, con, delrdestab] = pick_obj_and_con(susys, freq, cfreq, rnames, dnames)
	nr = numel(rnames);
	if isempty(dnames) % no dynamic uncertainty
		obj = @(dels)(eval_obj_par(susys, freq, rnames, dels));
	else % mixed uncertainty
		obj = @(dels)(eval_obj_mixed(susys, freq, rnames, dels));
	end
	if isempty(cfreq)
		con = []; % no constraints
		delrdestab = [];
	else
		if isempty(dnames) % no dynamic uncertainty
			con = @(dels)(eval_con_unstab_par(susys, cfreq, rnames, dels));
			[~, destab_unc] = robstab(susys); % gives the wrong result when called with the critical
			% frequency specified (probably due to a bug)
		else % mixed uncertainty
			con = @(dels)(eval_con_unstab_mixed(susys, cfreq, rnames, dels));
			[~, destab_unc] = robstab(susys, cfreq);
		end
		delrdestab = zeros(nr, 1);
		for kk = 1 : nr
			delrdestab(kk) = destab_unc.(rnames{kk});
		end
	end
end

function obj = eval_obj_mixed(susys, freq, rnames, dels)
	% The objective function is the sum of the worst-case gain
	% lower bounds against the dynamic uncertainty at the
	% specified frequencies.
	for kk = 1 : numel(rnames)
		reals.(rnames{kk}) = dels(kk);
	end
	lb = wcgainlbgrid(usubs(susys, reals), freq);
	obj = -sum(lb);
	if ~isfinite(obj)
		obj = -1e10;
	end
end

function obj = eval_obj_par(susys, freq, rnames, dels)
	% The objective function is the sum of the largest singular
	% values at the specified frequencies.
	for kk = 1 : numel(rnames)
		reals.(rnames{kk}) = dels(kk);
	end
	sys = usubs(susys, reals);
	gain = sigma(sys, freq);
	gain = gain(1, :);
	obj = -sum(gain);
	if ~isfinite(obj)
		obj = -1e10;
	end
end

function [c, ceq] = eval_con_unstab_mixed(usys, cfreq, rnames, dels)
	ceq = 0; % no equality constraints
	target_mu = 1;
	for kk = 1 : numel(rnames)
		reals.(rnames{kk}) = dels(kk);
	end
	sm = robstab(usubs(usys, reals), cfreq);	
	c = target_mu - 1 / sm.LowerBound; % mu(cfreq) >= target_mu
end

function [c, ceq] = eval_con_unstab_par(usys, cfreq, rnames, dels)
	ceq = 0; % no equality constraints
	target_mu = 1;
	for kk = 1 : numel(rnames)
		reals.(rnames{kk}) = dels(kk);
	end
	sys = usubs(usys, reals);
	% the real part of at least one pole must be non-negative
	c = -max(real(eig(sys))); % max(Re(p)) >= 0
end

function [dnames, rnames, nr] = process_system(usys)
	% check usys
	if ~isa(usys, 'uss')
		error('The first input argument must be ''uss''');
	end
	% find uncertainty names
	unames = fieldnames(usys.Uncertainty);
	if isempty(unames)
		error('The system has no uncertainty blocks.');
	end
	nr = 0;
	nd = 0;
	rnames = {};
	dnames = {};
	for kk = 1 : numel(unames)
		uname = unames{kk};
		if isa(usys.Uncertainty.(uname), 'ureal')
			nr = nr + 1;
			rnames{nr, 1} = uname;
		elseif isa(usys.Uncertainty.(uname), 'ultidyn')
			nd = nd + 1;
			dnames{nd, 1} = uname;
		end
	end
end

function [susys, cfreq, info] = scale_for_robust_stability(usys, info)
	% check robust stability
	sm = robstab(usys);
	if sm.LowerBound <= 1 % the system is not robustly stable => scale the unceratinty block
		[M, delta] = lftdata(usys);
		info.uscale = 1.01 * sm.UpperBound;
		susys = lft(delta * info.uscale, M);
		cfreq = sm.CriticalFrequency;
	else
		info.uscale = 1;
		cfreq = [];
		susys = usys;
	end
	% save frequencies
	info.susys = susys;
	info.cfreq = cfreq;
end

function [freq, info] = process_freq(freq, susys, info)
	% check frequencies
	if ~isnumeric(freq)
		error('The array of frequencies must be a numeric array.');
	end
	% replace zero frequency with a small positive number
	min_freq = 1e-6;
	freq(freq < min_freq) = min_freq;
	if numel(freq) ~= numel(unique(freq))
		freq = unique(freq);
		warning('Repeated frequencies removed.');
	end
	freq = reshape(freq, [1, numel(freq)]);
	freq = sort(freq, 'ascend');
	% remove the frequency points where mu >= 1 from the set of frequencies
	[~, ~, rinfo] = robstab(susys, freq);
	freq(rinfo.Bounds(:, 1) <= 1) = [];
	info.freq = freq;
end

function delopt = hypercube_interval_search(n_dim, obj, con, delcon, n_eval_max)
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
% 			error('No feasible point found.');
			distances = diag((dels - delcon)' *(dels - delcon));
			[~, distsortind] = sort(distances);
			closestind = distsortind(1);
% 			if closestind == 1
% 				closestind = distsortind(2);
% 			end
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
		end
	end
	delopt = dels(:, minind);
end

function combs = allcomb(varargin)
	n = numel(varargin);
	ndgrid_output = cell(n, 1);
	[ndgrid_output{:}] = ndgrid(varargin{:});
	combs = cellfun(@(m)(reshape(m, [1, numel(m)])), ndgrid_output, 'UniformOutput', false);
	combs = cell2mat(combs);
end