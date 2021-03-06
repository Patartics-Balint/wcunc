function [F, info] = bnpinterp(varargin)
	% Construct LTI system F such that
	%   1) F is stable
	%   2) hinfnorm(F) <= 1
	%   3) F(j * freq(i)) * v(i) = u(i) for u(i) * v(i)' = data(:, :, i)
	%      for scalar data this is equivalent to F(j * freq(i)) = data(i)
	%
	% Usage:
	%   [F, info] = BNPINTERP(data, freq)
	%   [F, info] = BNPINTERP(data, freq, G)
	% Inputs:
	%   data: Data to be interpolated. Eigther a vector of scalars with
	%   modulus one or a 3D array of rank-one matrices with largest
	%   singular value equal to one.
	%   freq: Frequencies at which data is to be interpolated. Vector of
	%   positive real numbers without duplicated elements.
	%   G: Stable system with hinfnorm(G) <= 1. If data is an array of
	%   matrices, than G is n-by-n where n is the bigger of the comlumn and
	%   row dimensions of data(:, :, 1). When ommited G = eye(n) is used.
	% Outputs:
	%   F: Interpolant satisfying conditions 1)-3).
	%   info: Structure containing information corresponding to the
	%   interpolation procedure.
	%
	% The algorithm is the tangential matricial boundary Nevalinna-Pick
	% inerpolation in Corollary 21.4.2. on page 481 and Example 21.3.2 on
	% page 476 in Ball-Gohberg-Rogman Interpolation of Rational Matix
	% Functions.
	%
	% See also cnum2sys, wcgain
	
	[C0minus, C0plus, freq, G, nr, nc, n, info] = parse_input(varargin);
	
	% add complex conjugates
	z = freq * 1i;
	z = [conj(z), z];
	C0minus = [conj(C0minus), C0minus];
	C0plus = [conj(C0plus), C0plus];
	zero_inds = find(z == 0);
	if ~isempty(zero_inds) % do not duplicate zero
		z(zero_inds(1)) = [];
		C0minus(:, zero_inds(1)) = [];
		C0plus(:, zero_inds(1)) = [];
	end

	% Pick matrix and the minimisation of derivatives
	H = zeros(numel(z));
	for ii = 1 : numel(z)
		for jj = 1 : numel(z)
			if ii ~= jj
				H(ii, jj) = (C0minus(:, ii)' * C0minus(:, jj) - C0plus(:, ii)' * C0plus(:, jj))...
					/ (conj(z(ii)) + z(jj));
			end
		end
	end % at this point the main diagonoal of H is zero
	rho = sdpvar(1, numel(freq), 'full');
	rho_peak = sdpvar(1);
	if isempty(zero_inds) % zero is not among the frequency points
		H = H + diag([rho, rho]);
	else
		H = H + diag([rho(2 : end), rho]);
	end
	epsilon = 1e-6;
	opt = sdpsettings('solver', 'lmilab', 'verbose', 0);
	diagn = optimize([rho_peak >= rho,...
		rho >= 0,...
		epsilon <= H],...
		rho_peak + sum(rho), opt);
	rho = value(rho);
	H = value(H);
	if strfind(diagn.info, 'Infeasible')
		error('The boundary Pick matrix could not be found.');
	end
	if any(eig(H) <= 0)
		error('The boundary Pick matrix is not positive definite.');
	end
	info.Pick = H;
	info.derivatives = rho;

	% construct Theta
	A0 = diag(z);
	B0 = inv(H) * [-C0plus', C0minus'];
	C0 = [C0plus; C0minus];
	D0 = eye(2 * n);
	[Ar, Br, Cr] = transform_ss_matrices_to_real(A0, B0, C0, z);
	Theta = ss(Ar, Br, Cr, D0);
	info.Theta = Theta;

	Theta11 = Theta(1 : n, 1 : n);
	Theta12 = Theta(1 : n, n + 1 : end);
	Theta21 = Theta(n + 1 : end, 1 : n);
	Theta22 = Theta(n + 1 : end, n + 1 : end);
	
	% construct interpolant
	Anum = [Theta11.a, Theta11.b * G.c; zeros(order(G), order(Theta11)), G.a];
	Bnum = [Theta12.b + Theta11.b * G.d; G.b];
	Cnum = [Theta11.c, G.c];
	Dnum = G.d;
	Cden = [Theta21.c, zeros(size(Theta21, 1), order(G))];
	F = ss(Anum - Bnum * Cden, Bnum, Cnum - Dnum * Cden, Dnum);
	F = sminreal(F);

	% check result
	if max(real(eig(F))) > 0
		error('The interpolant is unstable.');
	end
	inds = find(imag(z) >= 0);
	zchk = imag(z(inds));
	F0 = freqresp(F, zchk);
	u = C0minus(:, inds);
	v = C0plus(:, inds);
	for kk = 1 : size(u, 2)
		uu = u(:, kk);
		vv = v(:, kk);
		err(kk) = norm(F0(:, :, kk) * uu - vv);
	end
	info.error = err;
	if max(err) > 1e-5
		error('The fit is innacurate.');
	end
	if hinfnorm(F) > 1 + 1e-2
		error('The H-inf norm of the interpolant is greater than 1.');
	end
	
	% remove rows or colums corresponding to the zeros added to u and v
	F = F(1 : nr, 1 : nc);
end

function [C0minus, C0plus, freq, G, nr, nc, n, info] = parse_input(input)
	% sort frequencies
	freq = input{2};
	freq = reshape(freq, [1, numel(freq)]);
	if any(freq < 0)
		error('The frequencies must be non-negative.');
	end
	if numel(freq) ~= numel(unique(freq))
		error('Repeated frequencies are not allowed.');
	end
	[freq, sortind] = sort(freq);
	
	% sort out the interpolation data
	data = input{1};
	if iscell(data) % data contains left and right vectors
		C0plus = data{1};
		C0minus = data{2};
		C0plus = C0plus(:, sortind);
		C0minus = C0minus(:, sortind);
	elseif isnumeric(data) % data contains the matrices to be interpolated
		if size(data, 3) == numel(freq)
			data = data(:, :, sortind);
		elseif numel(data) == numel(freq)
			data = data(sortind);
		else
			error('The number of data samples and frequencies must be the same.');
		end
		data_is_scalar = (size(data, 3) == 1 &&...
			(size(data, 1) == 1 || size(data, 2) == 1) && numel(freq) > 1) ||...
			numel(data) == numel(freq);
		info.frequency = freq;
		info.data = data;
		if data_is_scalar
			data = reshape(data, [1, numel(data)]);		
		end
		if data_is_scalar
			C0minus = ones(size(data));
			C0plus = data;
		else		
			C0minus = zeros(size(data, 1), numel(freq));
			C0plus = zeros(size(data, 2), numel(freq));
			for kk = 1 : numel(freq)
				if rank(data(:, :, kk)) > 1
					warning('The data to be interpolated is not rank-one. Only the diad corresponding to the larges singular value is used.');
				end
				[U, S, V] = svd(data(:, :, kk));
				u = U(:, 1) * S(1, 1);
				v = V(:, 1);
				C0minus(:, kk) = v;
				C0plus(:, kk) = u;
			end
		end
	else
		error('The first input is either a numeric or a cell array.');
	end
	% add zeros when numel(u) ~= numel(v)
	nr = size(C0plus, 1);
	nc = size(C0minus, 1);
	n = max([nr, nc]);
	C0minus = [C0minus; zeros(n - nc, numel(freq))];
	C0plus = [C0plus; zeros(n - nr, numel(freq))];
	
	% find G
	if numel(input) < 3
		G = ss(eye(n));
	else
		G = ss(input{3});
		if ~isstable(G)
			error('G must be stable.');
		end
		if hinfnorm(G) > 1 + 1e-2
			error('The H-inf norm of G must be at most one.');
		end
	end
	info.G = G;
end

function [Ar, Br, Cr] = transform_ss_matrices_to_real(A0, B0, C0, z)
	nfreq = numel(z(imag(z) >= 0));
	if even(numel(z)) % no zero frequency
		It = eye(nfreq);
		T0 = [It, It; 1i * It, - 1i * It] / 2;	
		T0inv = [It, - 1i * It; It, 1i * It];		
	else % zero frequnecy
		It = eye(nfreq - 1);
		zc = zeros(nfreq - 1, 1);
		T0 = [It/2, zc, It/2;...
					zc', 1, zc';
					1i * It/2, zc, -1i * It/2];
		T0inv = [It, zc, -1i * It;...
						 zc', 1, zc';
						 It, zc, 1i * It];
	end
	Ar = T0 * A0 * T0inv;
	Br = T0 * B0;
	Cr = C0 * T0inv;
	maximag = @(M)(max(abs(imag(M(:)))));
	maxreal = @(M)(max(abs(real(M(:)))));
	tol = 1e-6;
	if maximag(Ar) > tol * maxreal(Ar) ||...
		 maximag(Br) > tol * maxreal(Br) ||...
		 maximag(Cr) > tol * maxreal(Cr)
		error('The state-space matrices of Theta could not be transformed to real matrices.');
	end
	Ar = real(Ar);
	Br = real(Br);
	Cr = real(Cr);
end