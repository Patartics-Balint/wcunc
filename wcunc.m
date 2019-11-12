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
	%   each uncertainty block.
	%
	% See also wcgain, bnpinterp
	
	[usys, freq, dnames, rnames, nr, info] = parse_input(usys, freq);
	
	if ~isempty(rnames) % there is real uncertainty in the system
		obj = @(dels)eval_obj(usys, freq, rnames, dels);
		delropt = fmincon(obj, zeros(nr, 1), [], [], [], [], -ones(nr, 1), ones(nr, 1));
		for kk = 1 : numel(rnames)
			wcur.(rnames{kk}) = delropt(kk);
		end
		usys_dyn = usubs(usys, wcur);
		wcu = wcunc(usys_syn, freq);
		for kk = 1 : numel(rnames)
			wcu.(rnames{kk}) = wcur.(rnames{kk});
		end
	end
	
	% calculate worst-case uncertainty block by block.
	[~, ~, wcginfo] = wcgain(usys, freq);
	for kblk = 1 : numel(dnames)
		blkname = dnames{kblk};
		delsamp = [];
		for kfr = 1 : numel(freq)
			delsamp(:, :, kfr) = freqresp(wcginfo.WorstPerturbation(kfr).(blkname), wcginfo.Frequency(kfr));
		end
		[wcublk, infoblk] = bnpinterp(delsamp, freq);
		info.(blkname) = infoblk;
		wcu.(blkname) = wcublk;
	end		
	gain = hinfnorm(usubs(usys, wcu));
	
	% check result
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
	info.frequency = freq;
	%find uncertainty names
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