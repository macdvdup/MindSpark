% QQDIAGRAM - Empirical quantile-quantile diagram.
%
% Description:
%               The quantiles (percentiles) of the input distribution Y are plotted (Y-axis)
%               against the corresponding quantiles of the input distribution X.
%               If only X is given, the corresponding quantiles are plotted (Y-axis)
%               against the quantiles of a Gaussian distribution ('Normal plot').
%               Two black dots indicate the lower and upper quartiles.
%               If the data in X and Y belong the same distribution the plot will be linear.
%               In this case,the red and black reference lines (.-.-.-.-) will overlap.
%               This will be true also if the data in X and Y belong to two distributions with
%               the same shape, one distribution being rescaled and shifted with respect to the
%               other.
%               If only X is given, a line is plotted to indicate the mean of X, and a segment
%               is plotted to indicate the standard deviation of X. If the data in X are normally
%               distributed, the red and black reference lines (.-.-.-.-) will overlap.
%
% Usage:
%   >>  ah  =  qqdiagram( x, y, pk );
%
% Inputs:
%   x       - vector of observations
%
% Optional inputs:
%   y       - second vector of observation to compare the first to
%   pk      - the empirical quantiles will be estimated at the values in pk [0..1]
%
% Author: Luca Finelli, CNL / Salk Institute - SCCN, 20 August 2002
%
% Reference: Stahel W., Statistische Datenanalyse, Vieweg, Braunschweig/Wiesbaden, 1995
%
% See also: 
%   QUANTILE, SIGNALSTAT, EEGLAB 

% Copyright (C) 2002 Luca Finelli, Salk/SCCN, La Jolla, CA
%
% Reference: Stahel, W. Statistische Datenanalyse, Vieweg, Braunschweig/Wiesbaden 1995
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function qqdiagram( x , y, pk )

if nargin < 1
	help qqdiagram;
	return;
end;	

if (nargin == 3 && (any(pk > 1) || any(pk < 0)))
    error('qqdiagram(): elements in pk must be between 0 and 1');
end

if nargin==1
	y=x;    
	nn=max(1000,10*length(y))+1;
	x=randn(1,nn);
end

if nargin < 3
	nx=sum(~isnan(x));
	ny=sum(~isnan(y));
   	k=min(nx,ny);
    pk=((1:k) - 0.5) ./ k;  % values to estimate the empirical quantiles at 
else 
    k=length(pk);
end

if nx==k
    xx=sort(x(~isnan(x)));
else
    xx=quantile(x(~isnan(x)),pk);
end

if ny==k
    yy=sort(y(~isnan(y)));
else
    yy=quantile(y(~isnan(y)),pk);
end

% QQ diagram
plot(xx,yy,'+')
hold on

% x-axis range
maxx=max(xx);
minx=min(xx);
rangex=maxx-minx;
xmin=minx-rangex/50;
xmax=maxx+rangex/50;

% Quartiles
xqrt1=quantile(x,0.25); xqrt3=quantile(x,0.75);
yqrt1=quantile(y,0.25); yqrt3=quantile(y,0.75);

plot([xqrt1 xqrt3],[yqrt1 yqrt3],'k-','LineWidth',2); % IQR range

% Drawing the line
sigma=(yqrt3-yqrt1)/(xqrt3-xqrt1);
cy=(yqrt1 + yqrt3)/2;
	
if nargin ==1
    maxy=max(y);
    miny=min(y);
    rangey=maxy-miny;
	ymin=miny-rangey/50;
	ymax=maxy+rangey/50;
	
	plot([(miny-cy)/sigma (maxy-cy)/sigma],[miny maxy],'r-.') % the line
    % For normally distributed data, the slope of the plot line
    % is equal to the ratio of the standard deviation of the distributions
	plot([0 (maxy-mean(y))/std(y)],[mean(y) maxy],'k-.') % the ideal line
	
	xlim=get(gca,'XLim');
	plot([1 1],[ymin  mean(y)+std(y)],'k--')
	plot([1 1],[mean(y)  mean(y)+std(y)],'k-','LineWidth',2)
        % textx = 1.0;
        % texty = mean(y)+3.0*rangey/50.0;
	% text(double(textx), double(texty),' St. Dev.','horizontalalignment','center')
    set(gca,'xtick',get(gca,'xtick'));  % show that vertical line is at 1 sd
	plot([0 0],[ymin  mean(y)],'k--')
	plot(xlim,[mean(y) mean(y)],'k--')
	% text(double(xlim(1)), double(mean(y)+rangey/50),'Mean X')
	plot([xqrt1  xqrt3],[yqrt1 yqrt3],'k.','MarkerSize',10)
	set(gca,'XLim',[xmin xmax],'YLim',[ymin ymax])
	xlabel('Standard Normal Quantiles')
	ylabel('X Quantiles')
else
    cx=(xqrt1 + xqrt3)/2;
    maxy=cy+sigma*(max(x)-cx);
	miny=cy-sigma*(cx-min(x));
	
	plot([min(x) max(x)],[miny maxy],'r-.'); % the line
    xlabel('X Quantiles');
    ylabel('Y Quantiles');
end
