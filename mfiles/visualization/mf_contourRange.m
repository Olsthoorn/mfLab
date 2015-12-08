function range=mf_contourRange(H,dH,Lay,fldname)
%MF_CONTOURRANGE computes suitable range of contour levels
%    gets a range min:dH:Max for the data in struct H to contour
%    H is obtained from readDat or readMT3D yielding H(i).values
%    dH is range interval. The number of significant didgets is also obtained
%    fldname is used or if unspecified filed 'values' is used
%    CountourRange(B,dPsi,'','Psi') will get contours from Budget file after
%    running mf_Psi to add the stream function first
%    from dH
%
% USAGE:
%    range=mf_contourRange(H,dH,Lay,fldname)
%    range=mf_contourRange(H,dH,Lay)
%    range=mf_contourRange(H,dH)
%
%    TO 100529 110319

%    Example:
%       A=peaks;
%       c=contourrange(0.5,A);
%       contour(peaks,c); colorbar
%
%    See also READDAT, READBUD, READMT3D, CONTOUR.

%    Copyright 1999-2013 Theo Olsthoorn, TUDelft, Waternet 
%    $Revision: 1 $  $Date: 2007/11/13 00:10:21 $


m=Inf; M=-Inf;

if nargin<2, error('Insufficient arguments'); end
if nargin<4, fldname='values'; end

if isstruct(H)
    if nargin<3 || ~exist('Lay','var') || isempty(Lay)
        m=+Inf; M=-Inf; dH=abs(dH);
        for i=1:length(H)
            m=min([m;H(i).(fldname)(:)]);
            M=max([M;H(i).(fldname)(:)]);
        end
    else
        for i=1:length(H)
            m=min(m,min(min(min(H(i).(fldname)(:,:,Lay)))));
            M=max(M,max(max(max(H(i).(fldname)(:,:,Lay)))));
        end
    end
else
    if nargin==2
        m=min(H(:));
        M=max(H(:));
    else
        m=min(min(min(H(:,:,Lay))));
        M=max(max(max(H(:,:,Lay))));
    end
end

N=(M-m)/dH+2;

if m<0, n=fix(m/dH)-1; else n=fix(m/dH); end

NMax=100; NMin=4;
range=dH*(n:n+N);
redfac=2;

fprintf('%s: length of range will be %d.\n',mfilename,length(range));
if length(range)>2500
    error('mfLab:mf_contourRrange:range',...
        'Your need more than 2500 contour lines, check your range %s',fldname);
end

if length(range)>100
    fprintf(['     However, it must be < %d, check your input! ',...
          'Lowest=%g step=%g Highest=%g yields %d contour lines\n'],NMax,m,dH,M,length(range));
    dH=redfac*dH;
    fprintf('     I will increase the contour stepsize by factor %g to %g to give about %.0f lines.\n',redfac,dH,length(range)/redfac);
    switch nargin
        case 2, range=mf_contourRange(H,dH);
        case 3, range=mf_contourRange(H,dH,Lay);
        case 4, range=mf_contourRange(H,dH,Lay,fldname);
    end
elseif length(range)<2
    fprintf(['     However, it should be > %d, check your input! ',...
          'Lowest=%g step=%g Highest=%g yields %d contour lines\n'],NMin,m,dH,M,length(range));
    dH=dH/redfac;
    fprintf('     I will reduce the contour stepsize by factor %g to %g to give about %.0f lines.\n',redfac,dH,length(range)*redfac);
    switch nargin
        case 2, range=mf_contourRange(H,dH);
        case 3, range=mf_contourRange(H,dH,Lay);
        case 4, range=mf_contourRange(H,dH,Lay,fldname);
    end
end

