% ERRORDLG2 - Makes a popup dialog box with the specified message and (optional)
%               title.
%
% Usage:
%   errordlg2(Prompt, Title);
%
% Example:
%   errordlg2('Explanation of error','title of error');
%
% Input:
%   Prompt  -   A text string explaining why the user is seeing this error message.
%   Title   _   A text string that appears in the title bar of the error message.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 12 August 2002
%
% See also: INPUTDLG2, QUESTDLG2

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, arno@salk.edu
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

function errordlg2(Prompt, Title);

if exist('beep') == 5
	beep;
else
	disp(char(7));
end
if nargin <2
	Title = 'Error';
end
questdlg2(Prompt, Title, 'OK', 'OK');
