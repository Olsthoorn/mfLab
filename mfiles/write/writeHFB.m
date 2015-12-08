function writeHFB(basename,hfb)
%WRITEHFB writes the input file for the HFB6 package (mf2k)
%
% Example:
%   writeHFB(basename,hfb)  --- writing discretization file
%
%   This function is called from mf_setup when HFB package is on in
%   worksheet NAM
%
% TO 101020 110429

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.',hfb.ext],'wt');

%0.
fprintf(fid,'# MATLAB  writeHFB %s\n',datestr(now));
fprintf(    '# MODFLOW writeHFB %s\n',datestr(now));

%1.
fprintf(fid,'%10d%10d%10d      %s\n',...
    hfb.NPHFB,hfb.MXFB,hfb.NHFBNP,'     NPHFB MXFP NHFBNP');
%
% NHFBNP is number of HFP lines. This number is a single value
% per simulationa and so has the same number of HFB lines has to 
% be repeated in every stress period.
% This is a very primitive way of input with tons of redundancy
% a bit of a shame in the programmers !
%
%2 skip
%3 skip
%4

I=1:hfb.NHFBNP;  % counter, HFB requires NHFBNP lines per stress period

if ~isempty(I)
    for iPer=1:hfb.NPER
        if size(hfb.HFB,2)==6
            fprintf(fid,'%10d%10d%10d%10d%10d %12g     L R C R C hydchr\n',hfb.HFB(I,:)');
        else
            I=find(hfb.HFB(:,1)==iPer);
            fprintf(fid,'%10d%10d%10d%10d%13d %12g     L R C R C hydchr\n',hfb.HFB(I,2:end)');
        end
        if hfb.NHFBNP>I(end),
            I=I+hfb.NHFBNP;
        else
            % repeat same HFB for next stress period
        end
    %5 Number of active HFB paramreters
        fprintf(fid,'%10d\n',hfb.NACTHFB);  % always zero because no parameters are used in mfLab
    %6
        for i=1:hfb.NACTHFB(min(length(hfb.NACTHFB),iPer))
            % Do nothing because no parameters are used in mfLab: hfb.NACHTFB(iPer).Pname(i);
        end
    end
end

fclose(fid);
