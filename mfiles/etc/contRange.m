function range=contRange(varargin)
%CONTRANGE computes a suitable range of values for contouring.
%
% USAGE:
%   range=ContourRange(H[,[N[[,Lay][,fldname1[,index]]]]]) 
%
%   The function tries to prevent unregular numbers.
%   You will not have exactly nMax values but a number between 4 and nMax.
%   [nMin nMax] is allowed instead of nMax to limit the range of contours.
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
nMin=10;

%% Input arguments
if isempty(varargin)
    error('%s: Not enough input data\n',mfilename);
end

[H    ,varargin] = getNext(varargin,{'cell','struct','double'},[]);
[nMax ,varargin] = getNext(varargin,'double',50);
[Lay  ,varargin] = getNext(varargin,'double',[]); 

[field, varargin] = getNext(varargin,{'cell','char'},'values');
if ~isempty(field)
    if iscell(field)
        fldname1 = field{1};
        fldname2 = field{end};
        index= getNext(varargin,'double',1);
    end
else
    fldname1 = field;
    fldname2 = field;
end


%%

nMin=min(nMin,round(nMax(  1)));
nMax=max(nMin,round(nMax(end)));

%%

if ~isempty(varargin)
    display(varargin);
    msgId = 'mfLab:ContourRange:varargin';
    warning(msgId,'on');
    warning(msgId,'%s: some elements of varargin not used\n',mfilename);
    warning(msgId,'off');
end

%% Use all layers or only specified ones
if isempty(Lay)
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
    H = H(:,:,Lay);
    
    m = min(H(:));
    M = max(H(:));
end


%% Compute range

% compute decent contour line range
dH = (M-m)/nMax;
d = 10^floor(log10(dH));
f = dH/d;
if     f>=7.50,            d=10 *d;
elseif f>=4.50 && f<7.50,  d=5  *d;
elseif f>=3.25 && f<4.50,  d=4  *d;
elseif f>=2.25 && f<3.25,  d=2.5*d;
elseif f>=1.50 && f<2.25,  d=2  *d;
else   % skip d=d;
end

M=d*(floor(M/d)+1); m = d*(floor(m/d));

range=m:d:M;

if isempty(range)
    warning(['%s: Your values range for field <<%s>> is empty.\n',...
    'Hence, there are no contours to plot!\n',...
    'Most probably your computed data are uniform. Please check and resolve this.\n',...
    'You may more easily track the error after switching on\n',...
    'Debug>Stopif Errors/Warnings>always stop if error\n',...
    'from the menu bar of the editor, and then run again.'],mfilename,fldname1);
end

