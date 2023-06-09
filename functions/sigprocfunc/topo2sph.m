% TOPO2SPH - convert a TOPOPLOT style 2-D polar-coordinate
%              channel locations file to a 3-D spherical-angle
%              file for use with HEADPLOT
% Usage: 
%   >> [c h] = topo2sph('eloc_file','eloc_outfile', method, unshrink);
%   >> [c h] = topo2sph( topoarray, method, unshrink );
%
% Inputs:
%   'eloc_file'    = filename of polar 2-D electrode locations file used by 
%                    TOPOPLOT. See >> topoplot example or CART2TOPO
%   'eloc_outfile' = output file of 3-D electrode locations in spherical angle 
%                    coords. for use in HEADPLOT.
%   topoarray      = polar array of 2-D electrode locations, with polar angle
%                    in the first column and radius in the second one.
%   method         = [1|2] 1 is for Besa compatibility, 2 is for
%                    compatibility with Matlab function CART2SPH. {default: 2}
%   unshrink       = [0<real<1] unshrink factor. Enter a shrink factor used
%                    to convert spherical to topo (see SPH2TOPO). Only 
%                    implemented for 'method' 1 (above). Electrode 'shrink' 
%                    is now deprecated. See >> help topoplot
% Outputs:
%   c = coronal rotation
%   h = horizontal rotation
%
% Author: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1999 
%
% See also: SPH2TOPO, CART2TOPO

% Copyright (C) 1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 3-16-00 changed name to TOPO2SPH for compatibility with CART2TOPO -sm
% 01-25-02 reformated help & license -ad 
% 03-22-02 complete remodeling for returning arguments and taking arrays -ad 

function [c, h] = topo2sph(eloc_locs,eloc_angles, method, unshrink)

MAXCHANS = 1024;

if nargin < 1
    help topo2sph;
    return;
end
if nargin > 1 && ~ischar(eloc_angles)
	if nargin > 2
		unshrink = method;
    end
	method = eloc_angles;
else
	method = 2;
end

if ischar(eloc_locs)
	fid = fopen(eloc_locs);
	if fid<1
	    fprintf('topo2sph()^G: cannot open eloc_loc file (%s)\n',eloc_locs)
	    return
	end
	E = fscanf(fid,'%d %f %f  %s',[7 MAXCHANS]);
	E = E';
	fclose(fid);
else
    E = eloc_locs;
    E = [ ones(size(E,1),1) E ];
end
    
if nargin > 1 && ischar(eloc_angles)
	if exist(eloc_angles)==2
	   fprintf('topo2sph: eloc_angles file (%s) already exists and will be erased.\n',eloc_angles);
	end

	fid = fopen(eloc_angles,'a');
	if fid<1
	    fprintf('topo2sph()^G: cannot open eloc_angles file (%s)\n',eloc_angles)
	    return
	end
end

if method == 2
	t = E(:,2); % theta
	r = E(:,3); % radius
	h = -t;  % horizontal rotation
	c = (0.5-r)*180;
else
	for e=1:size(E,1)
		% (t,r) -> (c,h)
		
		t = E(e,2); % theta
		r = E(e,3); % radius
		r = r*unshrink;
		if t>=0
			h(e) = 90-t; % horizontal rotation
		else
			h(e) = -(90+t);
		end
		if t~=0
			c(e) = sign(t)*180*r; % coronal rotation
		else
			c(e) = 180*r;
		end
    end
	t = t';
	r = r';
end

for e=1:size(E,1)
   if nargin > 1 && ischar(eloc_angles)
        chan = E(e,4:7);
        fprintf('%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
        fprintf(fid,'%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
   end 
end

