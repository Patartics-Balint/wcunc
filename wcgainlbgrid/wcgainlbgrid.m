function [wcglb, wcu] = wcgainlbgrid(sys,w)

	%% Get LFT Data
	% blk is nblk-by-3 matrix describing the uncertainty where 
	%    blk(i,:) = [rdim cdim nrep] 
	% describes the i^th uncertainty as rdim-by-cdim with nrep repetitions.
	% nrep<0 for repeated real parameters and nrep>0 for complex.
	[M,~,BLKSTRUCT] = lftdata(sys);
	nblk = numel(BLKSTRUCT);
	blk = zeros(nblk,3);
	for i=1:nblk
			blk(i,:) = [BLKSTRUCT(i).Size BLKSTRUCT(i).Occurrences];
			if isequal(BLKSTRUCT(i).Type,'ureal')
					blk(i,3) = -blk(i,3);
			end
	end   
	index = blk2idx(blk);

	% Call low-level ulftdata
	% tmp.Data_ is a ltipack.lftdataSS object. There is probably a way to
	% by-pass accessing this low level data/function but it is fine for now.
	warning off; 
	tmp = struct(sys);
	tmp = tmp.Data_;
	warning on;
	[~,sysB,~,sysmuB] = ulftdata(tmp);

	%% Set Options
	aMXGAIN = 1e7;
	mMXGAIN = 0;
	NTIMES = 4;  % # of times to alternate between power iter and coord ascent
	MAXCNT = 6;  % Max # of times to do complete coord ascent iter
	RNG = RandStream('mt19937ar');


	%% Compute lower bound at each frequency
	Nw = numel(w);
	wcglb = zeros(Nw,1);
	PertData = cell(Nw,1);
	for i=1:Nw
			% Get frequency response
			m = freqresp(M,w(i));

			% Call wclowc:
			% This performs coordinatewise maximization over real parameters and
			% power iteration on complex blocks.
			[wcglb(i),pertsreal,pertrepreal,pertLvec,pertRvec,pertrepcomp] = ...
					wclowc(m,index,NTIMES,MAXCNT,mMXGAIN,aMXGAIN,RNG);

			if nargout >= 2
				% Package Perturbation Data
				PertData{i} = localGetCriticalPert(index,nblk,...
						pertsreal,pertrepreal,pertLvec,pertRvec,pertrepcomp);
				wcu(i) = getWorstCasePerturbation(sysB, sysmuB, PertData{i}, w(i));
			end
	end
end


function PertData = localGetCriticalPert(bidx,nblk,...
    pertsreal,pertrepreal,pertLvec,pertRvec,pertrepcomp)
	% Packages worst-case perturbation data in cell-array format for
	% use by getWorstCasePerturbation.
	PertData = cell(nblk,1);

	for k=1:bidx.full.num
		 korig = bidx.full.origloc(k);
		 PertData{korig} = {pertLvec{k},pertRvec{k}};
	end
	for k=1:bidx.sreal.num
		 korig = bidx.sreal.origloc(k);
		 PertData{korig} = pertsreal(k);
	end
	for k=1:bidx.repreal.num
		 korig = bidx.repreal.origloc(k);
		 PertData{korig} = pertrepreal(k);
	end
	for k=1:bidx.repcomp.num
		 korig = bidx.repcomp.origloc(k);
		 PertData{korig} = pertrepcomp(k);
	end
end

function WCP = getWorstCasePerturbation(B,muB,PertData,freq)
	% Reconstructs worst-case block values from output of MUSSV. 
	% This function takes the block vector B, the MUSSV block description  
	% MUB, the worst-case perturbations PERTDATA, and optionally the 
	% frequency point FREQ (used to construct worst-case to construct 
	% worst-case values for ULTIDYN blocks). PERTDATA is either a block-
	% diagonal matrix of block values (see MUSSVUNWRAP) or a cell vector
	% of block values (see WCGBNB). This function can be used in batch 
	% mode by specifying a 3D array PERTDATA and a vector of frequencies 
	% FREQ. If NS=SIZE(PERTDATA,3), WCP is a NS-by-1 struct array with 
	% fields named after the blocks and containing the worst-case values 
	% for each block (matrix for static blocks, ltipack.ssdata for 
	% ULTIDYN blocks).
	%
	% Note: Assumes the blocks are sorted with sortBlocks.
	nblk = numel(muB.BlkInfo);
	ns = size(PertData,3);
	WCP = cell2struct(cell(nblk,ns),{muB.BlkInfo.BlockName},1);
	CellFlag = iscell(PertData);
	% Create BSIDX so that B(BSIDX(k)) is the "first" instance of 
	% the k'th block (BSidx(k) is not equal to k because of repeated 
	% blocks).
	BSidx = cumsum([1 muB.BlkInfo.NCopies]);  % (NBLK+1)-by-1
	for k=1:nblk
		Bk = B(BSidx(k));      % normalized LFTBlockWrapper
		BlockData = Bk.getBlockList{1}; % original atom (ureal, ucomplex, ultidyn)
		BlockName = getName(BlockData);
		DynamicFLag = isa(BlockData,'DynamicSystem');
		% Get perturbation values for this block
		if CellFlag
			 % WCGBNB format
			 nprt = PertData(k,:);
		else
			 % MUSSV format
			 szD = iosize(BlockData);
			 nprt = PertData(muB.BlkInfo(k).RowStart:muB.BlkInfo(k).RowStart+szD(1)-1,...
					muB.BlkInfo(k).ColStart:muB.BlkInfo(k).ColStart+szD(2)-1,:);
		end
		for ct=1:ns
			 if CellFlag
					prt = nprt{ct};
					if freq == 0 && DynamicFLag
						prtmat = prt{1} * prt{2};
						prtmat = real(prtmat);
						[u, s, v] = svd(prtmat);
						prt{1} = u(:, 1) * sqrt(s(1));
						prt{2} = v(:, 1) * sqrt(s(1));
					end
			 else
					prt = nprt(:,:,ct);
			 end
			 % Un-normalize perturbation
			 WCP(ct).(BlockName) = evalInverseBlockTransform(Bk,prt);
		end
	end
end