function [wcglb,wcu,critfreq] = wcgainlb(sys,w)

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
% AbsMax = 5;
% RelMax = 20;
% aMXGAIN = RelMax*norm(sys.Nominal,inf) + AbsMax;
% mMXGAIN = 0;
% NTIMES = 2;
% MAXCNT = 3;
aMXGAIN = 1e7;
mMXGAIN = 0;
NTIMES = 4;  % # of times to alternate between power iter and coord ascent
MAXCNT = 6;  % Max # of times to do complete coord ascent iter
RNG = RandStream('mt19937ar');


%% Compute lower bound at each frequency
Nw = numel(w);
wcv = zeros(Nw,1);
PertData = cell(Nw,1);
for i=1:Nw
    % Get frequency response
    m = freqresp(M,w(i));
    
    % Call wclowc:
    % This performs coordinatewise maximization over real parameters and
    % power iteration on complex blocks.
    [wcv(i),pertsreal,pertrepreal,pertLvec,pertRvec,pertrepcomp] = ...
        wclowc(m,index,NTIMES,MAXCNT,mMXGAIN,aMXGAIN,RNG);
    
    % Package Perturbation Data
    PertData{i} = localGetCriticalPert(index,nblk,...
        pertsreal,pertrepreal,pertLvec,pertRvec,pertrepcomp);    
end


%% Return largest lower bound found
[wcglb,idx] = max(wcv);
wcu = getWorstCasePerturbation(sysB,sysmuB,PertData{idx},w(idx));
critfreq = w(idx);


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
