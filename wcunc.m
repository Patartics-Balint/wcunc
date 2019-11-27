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
	
	[usys, freq, dnames, rnames, nr, info] = parse_input(usys, freq);
	[susys, freq, cfreq, info] = scale_for_robust_stability(usys, freq, info);
	
	if ~isempty(rnames) % there is parametric uncertainty in the system
		if isempty(cfreq)
			con = @(dels)(default_con(dels)); % default constraints (always satisfied)
			delr_init = zeros(nr, 1);
		else
			con = @(dels)(robstab_con(susys, cfreq, rnames, dels));
			[~, destab_unc] = robstab(susys, cfreq);
			for kk = 1 : numel(rnames)
				delr_init(kk) = destab_unc.(rnames{kk});
			end
		end
		obj = @(dels)(eval_obj(susys, freq, rnames, dels));
		fminopt = optimset('Display', 'off');
		[delropt, ~, exitflag] = fmincon(obj, delr_init, [], [], [], [],...
			-ones(nr, 1), ones(nr, 1), con, fminopt);
		if exitflag == -2
			error('Something went wrong. The worst-case value of the parametric uncertainty could not be found.');
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
		[~, ~, wcginfo] = wcgain(susys_dyn, freq);
		for kblk = 1 : numel(dnames)
			blkname = dnames{kblk};
			delsamp = [];
			for kfr = 1 : numel(freq)
				delsamp(:, :, kfr) = freqresp(wcginfo.WorstPerturbation(kfr).(blkname), wcginfo.Frequency(kfr));
			end
			if ~isempty(cfreq)
				stabsamp = freqresp(robinfo.WorstPerturbation.(blkname), cfreq);
				if ~isempty(freq)
					delsamp(:, :, end + 1) = stabsamp;
				else
					delsamp = stabsamp;
				end
			end
			[wcublk, infoblk] = bnpinterp(delsamp, [freq, cfreq], 0);
			info.(blkname) = infoblk;
			wcu.(blkname) = wcublk;
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
end

function check_result(usys, wcu, freq)
	% Check the accuracy of the worst-case dynamic uncertainty
	% construction.
	[~, ~, wcginfo] = wcgain(usys, freq);
	wcsys = usubs(usys, wcu);
	lb1 = wcginfo.Bounds(:, 1);
	resp = freqresp(wcsys, freq);
	for kk = 1 : size(resp, 3)
		lb2(kk, 1) = norm(resp(:, :, kk), 2);
	end
	if max(abs(lb2 - lb1)) > 1e-3
		error('Something went wrong. The gain of the system does not match the desired lower bound at the specified frequnecies.');
	end
end

function obj = eval_obj(usys, freq, rnames, dels)
	% The objective function is the sum of the worst-case gain
	% lower bounds against the dynamic uncertainty at the
	% specified frequencies.
	for kk = 1 : numel(rnames)
		reals.(rnames{kk}) = dels(kk);
	end
	[~, ~, info] = wcgain(usubs(usys, reals), freq);
	obj = -sum(info.Bounds(:, 1));
	if ~isfinite(obj)
		obj = -1e10;
	end
end

function [c, ceq] = default_con(dels)
	c = 0;
	ceq = 0;
end

function [c, ceq] = robstab_con(usys, cfreq, rnames, dels)
	ceq = 0; % no equality constraints
	target_mu = 1;
	for kk = 1 : numel(rnames)
		reals.(rnames{kk}) = dels(kk);
	end
	sm = robstab(usubs(usys, reals), cfreq);	
	c = target_mu - 1 / sm.LowerBound; % mu(cfreq) >= target_mu
end

function [usys, freq, dnames, rnames, nr, info] = parse_input(usys, freq)
	% check usys
	if ~isa(usys, 'uss')
		error('The first input argument must be ''uss''');
	end
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
	info.freq = freq;
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

function [susys, freq, cfreq, info] = scale_for_robust_stability(usys, freq, info)
	cfreq = [];
	% check robust stability
	sm = robstab(usys);
	if sm.LowerBound <= 1 % the system is not robustly stable => scale the unceratinty block
		[M, delta] = lftdata(usys);
		info.uscale = 1.01 * sm.LowerBound;
		susys = lft(delta * info.uscale, M);
		info.susys = susys;
		cfreq = sm.CriticalFrequency;
		% remove the frequency points where mu >= 1 from the set of frequencies
		[~, ~, rinfo] = robstab(susys, freq);
		freq(rinfo.Bounds(:, 1) <= 1) = [];
		% save frequencies
		info.cfreq = cfreq;
		info.freq = freq;
	end
end