function avi_compress(avifilename)
%AVI_COMPRESS compresses avifile on Mac systems
%
% USAGE:
%   avi_compress(avifilename)
%
% Compress avi file if possible. Function uses Mencoder which must be
% installed with mplayer from web.
%
% This piece of code was adapted form a suggestion on a user forum
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/170918
% Andrea Tagliassachi, 16 June 2008
%
% It Took some effort to install mplayer (which includes mencoder)
% but it worked tremendously well, thanks Andrea !! (800 MB --> 10 MB in a
% couple of seconds).
%
% TO 110325

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

outputfilename='output.avi';

if length(avifilename)<4 || ~strcmpi(avifilename(end-3:end),'.avi')
    avifilename=[avifilename '.avi'];
end

fprintf('Please wait, compressing avi file ...');

command = sprintf('/opt/local/bin/mencoder %s -ovc lavc -o %s',avifilename,outputfilename);

[STATUSVAR,STATUSMESSAGE] = unix(command);
    
if STATUSVAR ~= 0
   error( STATUSMESSAGE );
else
   delete(avifilename); % remove temp file
   movefile(outputfilename,avifilename);
   fprintf('compression completed!\n');
end