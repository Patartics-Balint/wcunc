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
	
	% parse inputs
	if ~isa(usys, 'uss')
		error('The first input argument must be ''uss''');
	end
	if ~isnumeric(freq)
		error('The array of frequencies must be a numeric array.');
	end
	if numel(freq) ~= numel(unique(freq))
		freq = unique(freq);
		warning('Repeated frequencies removed.');
	end
	unames = fieldnames(usys.Uncertainty);
	if isempty(unames)
		error('The system has no uncertainty blocks.');
	end
	freq = reshape(freq, [1, numel(freq)]);
	freq = sort(freq, 'ascend');
	info.frequency = freq;
	
	% calculate worst-case uncertainty block by block.
	[~, ~, wcginfo] = wcgain(usys, freq);
	for kblk = 1 : numel(unames)
		blkname = unames{kblk};
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