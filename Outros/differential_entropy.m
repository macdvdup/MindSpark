function [H] = differential_entropy(x,support,method,kNN)
%https://www.mathworks.com/matlabcentral/fileexchange/105945-differential-continuous-entropy-estimator/
%Evaluates the entopy for continuos random variables
%Input:
% x: Samples. Each row is an IID sample.
% support (optional): The support of the distribution.
%   A 2xD matrix. First row are lower limits.
%                  Second row are upper limits.
%                  [] if not known, or do not specify.
% method (optional): Estimation mthod.
% kNN  (optional): The number of neighbors or bins (when applicable)
%
%Output:
% H: Estimated entropy.

% Simple usage using default methods
%differential_entropy(x)
% or
%differential_entropy(x,support);


% Detailed usage:
% For univariate continuous distributions (D=1)
% x may be a row or column vector.

% If the support of the distribution is not-known or infinite. 
%    Using the default method (sample spacings):
%differential_entropy(x);
%    Choose method='spacings'|'bins'|'plugin'|'leave1out'|'average'|'vanEs'|'Correa':
%differential_entropy(x,method);

% For multivariate continuous distributions (D>1)
% If the support of the distribution finite and known, specify it as explained above. 
%    Using the default method (copula):
%differential_entropy(x,support);
%    Choose method='copula'|'KL'|'kNN'
%differential_entropy(x,support,method);

% Setting the spacings/number of nearest neighbors when applicable (methods 'spacings','vanEs','Correa','kNN'):
%differential_entropy(x,method,k);
% Or
%differential_entropy(x,[],method,k);
% 
% 
% Setting the number of bins when applicable (methods 'bins','plugin', 'leave1out','average'):
%differential_entropy(x,support,method,k);


%check input
if nargin==0
    error('differential_entropy: no input');
end
if ~isa(x,'double')
    error('differential_entropy: first input should be double');
end
if nargin==1
    support=[];
    method='default';
end
if nargin==2 
    if isa(support,'char') || isa(support,'string')
        method=support;
        support=[];
    elseif isa(support,'double')
        method='default';
    else
        error('differential_entropy: second input is invalid');
    end
end
if nargin==3 
    if isa(support,'char') || isa(support,'string')
        kNN=method;
        method=support;
        support=[];
    elseif isa(support,'double')
        kNN=[];
    else
        error('differential_entropy: third input is invalid');
    end
end
if any(support(:)==Inf) || any(support(:)==-Inf)|| any(isnan(support(:)))
    support=[];
end
if size(x,1)==1
    x=x';
end
if size(x,1)<size(x,2)
    error('differential_entropy: the number of samples (rows) has to be larger than the dimension (columns).');
end
if strcmp(method,'auto')
    method='default';
end

N=size(x,1); %number of samples
D=size(x,2); %dimension

%Algorithm parameters
%for univariate variables
params.nMaxBins=1000;       % Max number of bins for 1D entropy.
params.nMinBins=2;          % Min number of bins
params.nBinsExponent=0.4;   % Min number of bins for 1D entropy = n^params.nBinsExponent.
params.nBinsFraction=0.1;   % Min number of bins for 1D entropy = n*params.nBinsFraction.
params.nSpacingExponent=0.3;   % m_n for 1D entropy using spacings = n^params.nSpacingExponent.
if ~exist('kNN')
    params.kNN=0;
else
    if isempty(kNN)
        params.kNN=0;
    elseif kNN<0
        params.kNN=0;
    else
        params.kNN=round(kNN);
    end
end

%copula method
params.randomLevels=false;                       % Levels change randomly. Otherwise sequentially.
params.pValuesForIndependence=0.05;              % p-Value for independence using correlations.
params.exponentForH2D=0.62;                      % Scaling exponent for pair-wise entropy.
params.acceptenceThresholdForIndependence=-0.7;  % Threshold for pair-wide entropy to accept hypothesis of independence with p-Value 0.05.
params.minSamplesToProcessCopula=5*D;            % With fewer samples, assume the dimensions are independent.
params.numberOfPartitionsPerDim=2;               % For recursion: partition copula (along the chosen dimension) into numberOfPartitionsPerDim equal parts.


if D==1
    %univariate distributions
    %methods: spacings (default w/o support), bins, plugin, leave1out, average (default with support)
    if strcmp(method,'default')
        method='average';
    end
    if ~strcmp(method,'spacings') && ~strcmp(method,'bins') ...
             && ~strcmp(method,'plugin')  && ~strcmp(method,'leave1out') ...
              && ~strcmp(method,'average') && ~strcmp(method,'vanEs') ...
              && ~strcmp(method,'Correa')
          error('differential_entropy: the chosed 1D method is not implemented.');
    end
    
    if strcmp(method,'spacings') || isempty(support)
        %use method of sample spacings
        H=differential_entropy_H1Dspacings(sort(x),params);
    elseif strcmp(method,'bins')
        H=differential_entropy_H1Dbins(x,params,support);
    elseif strcmp(method,'plugin')
        H=differential_entropy_H1Dplugin(x,params,support);
    elseif strcmp(method,'leave1out')
        H=differential_entropy_H1Dleave1out(x,params,support);
    elseif strcmp(method,'vanEs')
        H=differential_entropy_H1DvanEs(sort(x),params);
    elseif strcmp(method,'Correa')
        H=differential_entropy_H1DCorrea(sort(x),params);
    elseif strcmp(method,'average') || strcmp(method,'default')
        H1=differential_entropy_H1Dbins(x,params,support);
        H2=differential_entropy_H1Dleave1out(x,params,support);
        H=(H1+H2)/2;
    else
        error('differential_entropy: 1D method not available.');
    end
    
    return;
    
end %if D==1

%D>1
if strcmp(method,'default')
    method='copula';
end
if ~strcmp(method,'copula') && ~strcmp(method,'KL') ...
        && ~strcmp(method,'K-L')  && ~strcmp(method,'kNN') ...
        && ~strcmp(method,'average')
    error('differential_entropy: the chosed 1D method is not implemented.');
end

if strcmp(method,'copula')
    H = differential_entropy_copulasH(x,params,support);
elseif strcmp(method,'KL') || strcmp(method,'K-L')
    H = differential_entropy_KL(x);
elseif strcmp(method,'kNN')
    H = differential_entropy_kNN(x,params);
else
    error('differential_entropy: 1D method not available.');
end

end


% copula method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H] = differential_entropy_copulasH(xs,params,support,level)
%Evaluates the entopy using recursive copula splitting.
%Input:
% xs: Samples. Each row is an iid sample.
% range: The support of the distribution.
%   A 2xD matrix. First row are lower limits.
%                  Second row are upper limits.
%                  [] if not known.
% level: Depth in the tree. Default is 0.
%
%Output:
% H: Estimated entropy.

if nargin<=2 % Default range value (unknown).
    support=[];
end
if nargin<=3 % Default level (=0).
    level=0;
end
D=size(xs,2);   % The dimension.
n=size(xs,1);   % Number of samples.
if n==0
    H=0;
    return;
end

%Calculate the empirical cumulative distribution along each dimention
%   (i.e., the integral transform of marginals).
H1s=0;
ys=zeros(n,D);
y1=zeros(n,1);
for d=1:D
    x1=xs(:,d);
    [x1s,idx]=sort(x1);
    
    % 1D entropy of marginals
    if isempty(support)
        % The range is not known, use m_n spacing method.
        H1=differential_entropy_H1Dspacings(x1s,params);
    else
        % The range is given, use binning.
        H1=(differential_entropy_H1Dbins(x1s,params,[support(1,d) support(2,d)])+differential_entropy_H1Dleave1out(x1s,params,[support(1,d) support(2,d)]))/2;
    end
    H1s=H1s+H1;
    
    % Rank of marginal.
    y1(idx)=( (0:n-1)+0.5 )/n;
    ys(:,d)=y1';
end

clear x1 xs; %No need for these any longer.
if n<params.minSamplesToProcessCopula || D==1
    H=H1s;
    return;
end

[R,P]=corr(ys); % Calculate the Peterson correlation coefficient of all pairs. 
                % Since the ys are the ranks, it is equal to the Spearman correlation of the xs.
                % P is the P-value matrix of the hypothesis that dimension pairs are independent.
P(eye(D)==1)=0;
isCorrelated=P<params.pValuesForIndependence; % A Boolean matrix. ==true is the corresponding pair is correlated.

if sum(isCorrelated(:))<D^2
    % Do more checks for pairs that are not correlated
    [ci,cj]=find(isCorrelated==0);
    nFactor=n^params.exponentForH2D;
    for k=1:length(ci)
        i=ci(k);
        j=cj(k);
        if i>j
            %Calculate pair-wide entropy = mutual information (because marginals are U(0,1).
            H2=differential_entropy_H2D([ys(:,i),ys(:,j)],params);
            isCorrelated(i,j)=H2*nFactor<params.acceptenceThresholdForIndependence;
            isCorrelated(j,i)=isCorrelated(i,j);
        end
    end
end %if sum(isCorrelated(:))<D^2

% Partition isCorrelated into blocks
nCorrelated=sum(isCorrelated);
L=isCorrelated-diag(nCorrelated);
Z=null(L,'r');
nBlocks=size(Z,2);
if nBlocks>=2
    H=H1s;
     %Split into blocks
    for c=1:nBlocks
        clusterSize=sum(Z(:,c));
        if clusterSize>1 % Blocks of size 1 are a marginal. The distribution is U(0,1), so its entropy is 0.
            dimsInCluster=find(Z(:,c)==1);
            ysmall=ys(:,dimsInCluster);     % Samples of the block
            smallD=length(dimsInCluster);   % The dimension of the block
            unitCube=[zeros(1,smallD) ; ones(1,smallD)];  % The support is always [0,1]^smallD.
            Hpart=differential_entropy_copulasH(ysmall,params,unitCube,level+1);      % Calculate the entropy of the reduced block.
            
            % Add the entropy of the block
            H=H+Hpart;
        end
    end
    
    return;
end

if n<params.minSamplesToProcessCopula
    H=H1s;
    return;
end

% Find which dimension to is most correlated ith others.
R2=R.^2;
R2(eye(D)==1)=0;
Rsum=sum(R2);
largeDims=find(max(Rsum)-Rsum<0.1 & Rsum>0); % List of dimenstions which are highly correlated with others (can be 0, 2, 3, ..., but not 1).
if isempty(largeDims)
    % Correlations are small, yet variables are dependent. All dims are equally problematic.
    largeDims=1:D;
end
nLargeDims=length(largeDims);
% Pick one of the dims in largeDims
if params.randomLevels
    % Choose randomly.
    maxCorrsDim=largeDims(randi(nLargeDims));
else
    % Choose sequentially.
    maxCorrsDim=largeDims(mod(level,nLargeDims)+1);
end

unitCube=[zeros(1,D) ; ones(1,D)];
Hparts=zeros(params.numberOfPartitionsPerDim,1);
%Split the data along Dim maxCorrsDim into params.numberOfPartitionsPerDim equal parts.
for prt=1:params.numberOfPartitionsPerDim
    % Range of data is [f,l)
    f=(prt-1)/params.numberOfPartitionsPerDim;
    if prt==params.numberOfPartitionsPerDim
        l=2;
    else
        l=prt/params.numberOfPartitionsPerDim;
    end
    
    mask=ys(:,maxCorrsDim)>=f & ys(:,maxCorrsDim)<l;
    nIn1=sum(mask);
    maskM=repmat(mask,[1 D]);
    
    % Subset of data.
    y1=reshape(ys(maskM),[nIn1,D]);
    % Scale back to [0,1].
    y1(:,maxCorrsDim)=(y1(:,maxCorrsDim)-f)*params.numberOfPartitionsPerDim;
    % Entropy of subset.
    Hparts(prt)=differential_entropy_copulasH(y1,params,unitCube,level+1);
end

% Add the entropies of all subsets.
H=H1s+sum(Hparts)/params.numberOfPartitionsPerDim;

end

function [H] = differential_entropy_H2D(xs,params)
% Calculate the entropy of 2D data compatly supported in [0,1]^2.
n=size(xs,1);
if n<4
    H=0;
    return;
end
nbins=max([min([params.nMaxBins floor(n^(params.nBinsExponent/2)) floor(params.nBinsFraction*n)]) params.nMinBins]);  % Number of bins.
pys=floor(xs*nbins);         % The partition number in each dimension.
py1=pys(:,1)+pys(:,2)*nbins; % Enumerate the partitions from 1 to npartitions^2.
edges=-0.5:(nbins^2-0.5);
counts=histcounts(py1,edges);% 1D histogram
mask=counts>0;
H=-sum( counts(mask).*log(counts(mask)) )/n + log(n/nbins/nbins);
end


% The original Kozachenko-Leonenko estimator

function [H] = differential_entropy_KL(xs)
    N=size(xs,1);
    D=size(xs,2);
    [~,NNdist]=knnsearch(xs,xs,'K',2);
    Ri=NNdist(:,2);
    H=log(N-1)+D*mean(log(Ri))+log(pi^(D/2)/gamma(D/2+1)) + 0.577;
end

% kNN - similar to Kozachenko-Leonenko 

function [H] = differential_entropy_kNN(xs,params)
    N=size(xs,1);
    D=size(xs,2);
    if params.kNN==0
        mn=max([floor(1*N^params.nSpacingExponent) 1]);
    else
        mn=params.kNN;
    end
    [~,NNdist]=knnsearch(xs,xs,'K',mn+1,'IncludeTies',false);
    %M=cell2mat(NNdist);
    Ri=max(NNdist,[],2);
    H=psi(N)-psi(mn)+log(pi^(D/2)/gamma(D/2+1)) + D*mean(log(Ri));
end

% 1D estimatiors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H] = differential_entropy_H1Dspacings(xs,params)
%Calculate the entropy of 1D data using m_N spacings.
%   Assume xs is ordered.
n=length(xs);
if params.kNN==0
    mn=max([floor(0.2*n^params.nSpacingExponent) 1]);
else
    mn=params.kNN;
end
mnSpacings=xs(mn+1:end)-xs(1:end-mn);
H=( sum(log(mnSpacings)) + log(n/mn)*(n-mn) )/n + log(mn) - psi(mn); %+ 1/2/mn + 1/12/mn^2;
end %differential_entropy_H1Dspacings

function [H] = differential_entropy_H1Dbins(xs,params,support)
%Calculate the entropy of 1D data using bins.
n=length(xs);
if params.kNN==0
    nbins=max([min([params.nMaxBins floor(n^params.nBinsExponent) floor(params.nBinsFraction*n)]) params.nMinBins]);  % Number of bins.
else
    nbins=params.kNN;
end
if nargin==2
    % No support.
    [p,edges]=histcounts(xs,nbins,'Normalization','pdf');
else
    % Range given.
    edges=linspace(support(1),support(2),nbins+1);
    p=histcounts(xs,edges,'Normalization','pdf');
end
mask=p>0;
H=-(edges(2)-edges(1))*sum( p(mask).*log(p(mask)) );
end %differential_entropy_H1Dbins

function [H] = differential_entropy_H1Dplugin(xs,params,support)
%Calculate the entropy of 1D data using bins.
n=size(xs,1);
if params.kNN==0
    nbins=max([min([params.nMaxBins floor(n^params.nBinsExponent) floor(params.nBinsFraction*n)]) params.nMinBins]);  % Number of bins.
else
    nbins=params.kNN;
end
if nargin==2
    support=[0 1];
end
edges=linspace(support(1),support(2),nbins+1);
mask=randi(2,n,1)==1;
n1=sum(mask);
n2=n-n1;
x1(:)=xs(mask);
x2(:)=xs(~mask);
p=histcounts(x1,edges,'Normalization','pdf');
c=histcounts(x2,edges);
p(p==0)=1;
H=-sum( c.*log(p) )/n2;
end %differential_entropy_H1DbinsPlugin

function [H] = differential_entropy_H1Dleave1out(xs,params,support)
n=size(xs,1);
if params.kNN==0
    nbins=max([min([params.nMaxBins floor(n^params.nBinsExponent) floor(params.nBinsFraction*n)]) params.nMinBins]);  % Number of bins.
else
    nbins=params.kNN;
end
if nargin==2
    support=[0 1];
end
edges=linspace(support(1),support(2),nbins+1);
de=edges(2)-edges(1);
c=histcounts(xs,edges);
c(c==0)=2;
c(c==1)=2;
H=-sum(c.*log(c-1))/n+log((n-1)*de);
end %differential_entropy_H1Dleave1out

function [H] = differential_entropy_H1DvanEs(xs,params)
%Calculate the entropy of 1D data using m_N spacings.
%   Assume xs is ordered.
n=length(xs);
if params.kNN==0
    mn=max([floor(0.2*n^params.nSpacingExponent) 1]);
else
    mn=params.kNN;
end
mnSpacings=xs(mn+1:end)-xs(1:end-mn);
H= mean(log(mnSpacings)) + sum((mn:n).^(-1));
end %differential_entropy_H1Dspacings

function [H] = differential_entropy_H1DCorrea(xs,params)
%Calculate the entropy of 1D data using m_N spacings.
%   Assume xs is ordered.
n=length(xs);
if params.kNN==0
    mn=max([floor(0.1*n^params.nSpacingExponent) 1]);
else
    mn=params.kNN;
end
mn2=2*mn+1;
xExt1=[ones(mn,1)*xs(1) ; xs ; ones(mn,1)*xs(1)];
xExt2=[ones(mn,1)*xs(1) ; xExt1 ; ones(mn,1)*xs(1)];
Xbar = movmean(xExt2,mn2,'Endpoints','discard');
Xcentered=xExt1-Xbar;
Xvar = movmean(Xcentered.^2,mn2,'Endpoints','discard');
top=conv(xExt1,-mn:mn,'valid'); 
bi=top./Xvar/n;
H=-mean(bi);
end %differential_entropy_H1Dspacings
