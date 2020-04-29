function [wcu, gain, info] = wcunc(usys, freq)
	% Calculate worst-case uncertainty which maximises the gain of an
	% uncertain system at specified frequnecies.
	%
	% Usage
	%   [wcu, gain, info] = WCUNC(usys, freq)
	% Inputs:
	%   usys: uss
	% Outputs:
	%   wcu: Worst-case uncertainty. The gain of usubs(usys, wcu) is maximal
	%   at the frequencies in freq.
	%   gain: hinfnorm(usubs(usys, wcu))
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
		[obj, con, delr_init] = pick_obj_con_and_init_val(susys, freq, cfreq, rnames, dnames);
		fminopt = optimset('Display', 'off');
		[delropt, ~, exitflag] = fmincon(obj, delr_init, [], [], [], [],...
			-ones(nr, 1), ones(nr, 1), con, fminopt);
		if exitflag == -2
			error('The worst-case value of the parametric uncertainty could not be found.');
		end
		for kk = 1 : numel(rnames)
			wcu.(rnames{kk}) = delropt(kk);
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
			for kfr = 1 : numel(freq)
				delsamp(:, :, kfr) = wcpert(kfr).(blkname);
			end
			if ~isempty(cfreq)
				stabsamp = freqresp(robinfo.WorstPerturbation.(blkname), cfreq);
				if ~isempty(freq)
					delsamp(:, :, end + 1) = stabsamp;
				else
					delsamp = stabsamp;
				end
			end
			[wcublk, infoblk] = bnpinterp(delsamp, [freq, cfreq]);
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
	% compute the gain
	gain = hinfnorm(usubs(usys, wcu));
end
	
function freq = pick_freq_grid(susys, info)	
	freqs = [0, logspace(-6, 6, 150)]';
	if info.uscale == 1
		wcg = wcgain(susys);
		freqs = unique([wcg.CriticalFrequency; freqs]);
	else
		[~, ~, robinfo] = robstab(susys, freqs);
		freqs(robinfo.Bounds(:, 1) <= 1) = [];	
	end
	wcglb = wcgainlbgrid(susys, freqs);
	if info.uscale == 1
		[~, maxind] = max(wcglb);
		freq = freqs(maxind);
	else
		freq = [];
	end
	range = max(wcglb) - min(wcglb);
	if info.uscale < 1
		range = min([range, 100]);
	end
	[pks, pfreq] = findpeaks(wcglb, freqs, 'Threshold', 0.01 * range);
	if isempty(pfreq)
		[pks, pfreq] = findpeaks(wcglb, freqs, 'NPeaks', 1);
	end
	if isempty(pfreq)
		keyboard;
	end
	if ~ismember(0, pfreq) && ismember(0, freqs)
		% add zero
		pks = [wcglb(1); pks];
		pfreq = [freqs(1); pfreq];
	end
	if numel(pks) < 4
		freq = [freq; pfreq];
	else
		[~, sortinds] = sort(pks, 'descend');
		freq = [freq; pfreq(sortinds(1 : 4))];
	end
	freq = unique(freq);
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
	if max(abs(lb2 - lb1)) > 1e-3
		error('The gain of the system does not match the desired lower bound at the specified frequnecies.');
	end
end

function [obj, con, delr_init] = pick_obj_con_and_init_val(susys, freq, cfreq, rnames, dnames)
	if isempty(dnames) % no dynamic uncertainty
		obj = @(dels)(eval_obj_par(susys, freq, rnames, dels));
	else % mixed uncertainty
		obj = @(dels)(eval_obj_mixed(susys, freq, rnames, dels));
	end
	if isempty(cfreq)
		con = @(dels)(eval_con_stab(dels)); % no constraints (always satisfied)
		delr_init = zeros(numel(rnames), 1);
	else
		if isempty(dnames) % no dynamic uncertainty
			con = @(dels)(eval_con_unstab_par(susys, cfreq, rnames, dels));
			[~, destab_unc] = robstab(susys); % gives the wrong result when called with the critical
			% frequency specified (probably due to a bug)
		else % mixed uncertainty
			con = @(dels)(eval_con_unstab_mixed(susys, cfreq, rnames, dels));
			[~, destab_unc] = robstab(susys, cfreq);
		end			
		for kk = 1 : numel(rnames)
			delr_init(kk) = destab_unc.(rnames{kk});
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

function [c, ceq] = eval_con_stab(dels)
	c = 0;
	ceq = 0;
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