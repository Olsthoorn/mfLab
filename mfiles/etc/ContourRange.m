function range=ContourRange(varargin)
%CONTOURRANGE computes a suitable range of values for contouring.
%
% USAGE:
%   range=ContourRange(H,NMax,Lay,fldname1,index) 
%
%   The function tries to prevent unregular numbers.
%   You will not have exactly NMax values but a number between 4 and NMax.
%   [NMin NMax] is allowed instead of NMax to limit the range of contours.
%
%   H may be a struct with a field holding 3D values, such as that
%   obtained from readData or readMT3D.
%   The fldname must be specified if it is not the default: 'values' such as
%   it is the case with the results of readBud and readBud6.
%
%   The index must be specified if the fldname1 refers to a cell array holding
%   more than a single 3D array, such as is the case for the field 'term'
%   obtained from readBud and readBud6.
%
%   H may be a cell array. In that case the seconc struct array or value will
%   be subtracted from the first struct array before the range is computed.
%
%   H may also be a numeric 3D array. In that case not fldname1 and index should
%   be specified.
%
% USAGE:
%  range=ContourRange(B,'Psi');
%  range=ContourRange(B,N,'Psi');
%  range=ContourRange(B,N,iLay,'Psi');
%         where N hte number of desired
%         contour values, iLay the layer number and 'Psi' the field name, in
%         this specific case 'Psi'.
%  range=ContourRange(B,'term')
%  range=ContourRange(B,'term',index)
%  range=ContourRange(B,N,'term',index)
%  range=ContourRange(B,N,iLay,'term',index)
%         same as before, but the budget struct B has many term cells each being
%         a complete 3D aarray. In thata case add the index as last argument so
%         that the ContourRange will be based on B(:).term{index}(:,:,:).
%  range=ContourRange(H)
%  range=ContourRange(H,N)
%  range=ContourRange(H,N,Lay)
%  range=ContourRange({H,h0},N); -- first subtract h0 from H
%  range=ContourRange({H,h0},N,iLay); -- first subtract h0 from H
%
%  Example:
%       A=peaks;
%       crange=ContourRange(A,50);
%       contour(peaks,crange); colorbar
%
%  See also READDAT, READBUD, READMT3D, CONTOUR.
%
%  TO 100529 110319 130221 130322

%    Copyright 1999-2013 Theo Olsthoorn, TUDelft, Waternet 
%    $Revision: 1 $  $Date: 2007/11/13 00:10:21 $

%% default number of contours
NMin=10;

%% Input arguments
if isempty(varargin)
    error('ContourRange: Not enough input data\n');
end

if ~isempty(varargin) && (iscell(varargin{1}) || isstruct(varargin{1}) || isnumeric(varargin{1})) || islogical(varargin{1})
    H=varargin{1};
    varargin(1)=[];
else
    error(['%s: first arg must be a logical array, a numeric array, a strucure read from readDat, readMT3D or readBud\n',...
        '  or a cell array of these. The current is of illegal class <<%s>>'],...
        mfilename,class(varargin{1}));
end

%% NMax
if ~isempty(varargin) && isnumeric(varargin{1})
        NMax = varargin{1};
        varargin(1)=[];
else
    NMax=40;
end

NMin=min(NMin,round(NMax(  1)));
NMax=max(NMin,round(NMax(end)));

%% Lay
if ~isempty(varargin) && (isnumeric(varargin{1}) || isempty(varargin{1}))
    Lay = varargin{1};
    varargin(1)=[];
else
    Lay = [];
end

%% fldname1
if ~isempty(varargin)
    if ischar(varargin{1})
        fldname1  = varargin{1};
        fldname2 = varargin{1};
        varargin(1)=[];
    elseif iscell(varargin{1})
        fldname1  = varargin{1}{  1};
        fldname2 = varargin{1}{end};
        varargin(1)=[];
    else
        fldname1  = 'values';
        fldname2 = 'values';
    end
else
    fldname1  = 'values';
    fldname2 = 'values';
end

%%
if ~isempty(varargin) && isnumeric(varargin{1})
    index = varargin{1};
    varargin = [];
else
    index=1;
end

if ~isempty(varargin)
    display(varargin);
    msgId = 'mfLab:ContourRange:varargin';
    warning(msgId,'on');
    warning(msgId,'%s: some elements of varargin not used\n',mfilename);
    warning(msgId,'off');
end

%% Use all layers or only specified ones
if ~exist('Lay','var') || isempty(Lay)
    if iscell(H)
        if isstruct(H{1})
            Lay = 1:size(H{1}(1).(fldname1),3);
        else
            Lay = 1:size(H{1},3);
        end
    else
        if isstruct(H)
            Lay = 1:size(H(1).(fldname1),3);
        else
            Lay = 1:size(H,3);
        end
    end
end

%% If difference is desired:
% in case H = {Var1,Var2} Var2 will be subtracted from Var1 where H is a
% struct with field values like from readMT3D and Var2 an idential struct
% or an array what will be subtracted.
if iscell(H)
    if isstruct(H{1})
        if isstruct(H{2})
            try
                Nt1=length(H{1});
                Nt2=length(H{2});
                for it=1:length(H{1})
                    H{1}(it).(fldname1) = H{1}(max(it,Nt1)).(fldname1) - H{2}(max(it,Nt2)).(fldname2);
                end
            catch ME
                fprintf(1,ME.message); fprintf(1,'\n');
                error('%s: both structs in first argument have compatible range');
            end
        elseif isnumeric(H{2})
            try
            for it=1:length(H{1})
                H{1}(it).(fldname1) = H{1}(it).(fldname1) - H{2};
            end
            catch ME
                fprintf(1,ME.message); fprintf(1,'\n');
                error('%s: both struct and array in first argument have compatible range');
            end
        else
            error(['%s: The first argument is a cell, it should have two compatible\n',...
                    'structs and/or 3D numerical arrays.'],mfilename);
        end        
    else
        if isstruct(H{2})
            try
                Nt2=length(H{2});
                for it=1:length(H{1})
                    H{2}(it).(fldname2) = H - H{2}(max(it,Nt2)).(fldname2);
                end
            catch ME
                fprintf(1,ME.message); fprintf(1,'\n');
                error('%s: aray and struct in first argument have compatible range');
            end
        elseif isnumeric(H{2})
            try
                H{1} = H{1} - H{2};
            catch ME
                fprintf(1,ME.message); fprintf(1,'\n');
                error('%s: both arrays in first argument have compatible range');
            end
        else
            error(['%s: The first argument is a cell, it should have two compatible\n',...
                    'structs and/or 3D numerical arrays.'],mfilename);
        end        
    end
    H = H{1};
end

%% Determine type of input
if isstruct(H)
    m=+Inf; M=-Inf;
    for i=1:length(H)
        if iscell(H(i).(fldname1))
            m=min(m,min(min(min(H(i).(fldname1){index}(:,:,Lay)))));
            M=max(M,max(max(max(H(i).(fldname1){index}(:,:,Lay)))));
        else
            m=min(m,min(min(min(H(i).(fldname1)(:,:,Lay)))));
            M=max(M,max(max(max(H(i).(fldname1)(:,:,Lay)))));
        end
    end
else
    m=min(min(min(H(:,:,Lay))));
    M=max(max(max(H(:,:,Lay))));
end

mm=m; MM=M;

%% Compute range
if m<0,
    m =-10^round(log10(-m*3));
elseif m>0
    m =+10^round(log10(+m/3));
else
    m=0;
end

if M<0,
    M =-10^round(log10(-M/3));
elseif M>0
    M =+10^round(log10(+M*3));
else
    M=0;
end

%dH    = 10^round(log10((M-m)/NMax));

f=floor((M-m)/(MM-mm)); % always > 1
 
if     f>=5  , f=5;
elseif f>=2.5, f=2.5;
elseif f>=2  , f=2;
else   f=1;
end

dH=(M-m)/NMax/f;

range=m:dH:M;   range=range(range>=mm & range<=MM);

if isempty(range)
    warning(['%s: Your values range for field <<%s>> is empty.\n',...
    'Hence, there are no contours to plot!\n',...
    'Most probably your computed data are uniform. Please check and resolve this.\n',...
    'You may more easily track the error after switching on\n',...
    'Debug>Stopif Errors/Warnings>always stop if error\n',...
    'from the menu bar of the editor, and then run again.'],mfilename,fldname1);
end

