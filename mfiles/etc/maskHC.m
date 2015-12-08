function H=maskHC(H,MASK,values,fieldName)
%MASKHC masks H(i).(fieldName) if  isstruct(H) or H(i).(fieldName({j} if isstruct(H(i).fieldName)
%
% USAGE:
%     H = maskCH(H,MASK,values,fieldname)
% 
%    H can be struct with fiel values or H is is a numeric 3D array.
%
%    Default fieldname is "values"
%
%    MASK may be 2 values or a MASK array where ~=0 is false.
%
%    works for structs obtained from readDat, readMT3D and readBud, for
%    arrays and struct arrays in general by specigying fieldName.
%
% Example:
%    H = maskHC(H,MASK,value,fieldName)
%        any values in H where MASK==0 are set to value or NaN if value is omitted
%    H = maskHC(H,[V1,V1],[w1 w2])
%        any H <min(V1,V2) --> w1; any H>max(V1,V2) --> w2
%    H = maskHC(H,[V1,V2],w1)
%        any H>min(V1,V2)&& H<max(V1,V2) -->w1
%    H = maskHC(H,[V1,V2])
%        identical to previous with w1=V1 and w2=V2;
%    H = maskHC(H,V1,w1)
%        any H ~=V1 -> w1  V1 full array of size H or H(i).values
%        this works to set concentrations from MT3DMS and SEAWAT to zero where
%        they are slightly below zero due to solver constraints.
%
% CONCRETE examples
%    H = maskHC(H,[-Inf 1e4],[NaN NaN]) removes H(i).values 10^10 and 10^30 typically
%        appllied in MODFLOW to indicate dry and inactive cells.
%    C = maskHC(C,[0,Inf],[0,Inf]); % sets any C(i).values<0 to zero. Which is
%        useful to remove concentractions that are a little below zero caused
%        by solver limitations.
%    B = maskHC(B,[-Inf 1e6],[-Inf NaN],'term'); would set any values in any
%        of the flowterms B(i).term{j}(:,:,:) larger than 1e6 to NaN
%
% See also maskBud readDat readMT3D readBud
%
% TO 090315 100420 100718 100911 120409 120509

% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<2, error('mflab:maskHC:insufficientInputArguments',...
    ['maskHK: 2 or 3 inpu arguments required:\n',...
        '%s(H,MASK,value)\n'],mfilename);
end

if nargin<3, values=[]; end
if nargin<4, fieldName='values'; end

if ~(isnumeric(MASK) || islogical(MASK))
    error('maskCH: if second argument is used, it must be scalar, two scalars or a full array of ize H');
end
if ~isnumeric(values)
    error('maskHC: if the third argument is used, it must be scalar, two sclars or a full array size of H');
end


if isnumeric(H) || islogical(H)
    H=dispatch(H,MASK,values);
elseif isstruct(H)
    for i=1:length(H)
        if iscell(H(i).(fieldName))  % this also works for Budget B(i).term{j}
            for ic=1:length(H(i).(fieldName))
                H(i).(fieldName){ic}=dispatch(H(i).(fieldName){ic},MASK,values);
            end
        else
                H(i).(fieldName)=dispatch(H(i).(fieldName),MASK,values);
        end
    end
else
    error('maskCH: H must be a numerical array or a struct as from readDat or readMT3D');
end


function H=dispatch(H,MASK,values)
    
    if isempty(values)
        if isscalar(MASK)
            H(H~=MASK)=NaN;
        elseif numel(MASK)==2
            MASK=sort(MASK);
            H(H<MASK(1))=MASK(1);
            H(H>MASK(2))=MASK(2);
        else
            if all(size(H)==size(MASK))
                H(~MASK)=NaN;
            else
                return;
            end
        end
    else
        if isscalar(MASK)
            if isscalar(values)
                H(H~=MASK)=values;
            elseif numel(values==2)
                H(H<min(values) && H>max(values))=MASK;
            else
                H(H~=MASK) = values(H~=MASK);
            end
        elseif numel(MASK)==2
            if isscalar(values)
                H(H<MASK(1))=values;
                H(H>MASK(2))=values;
            elseif numel(values==2)
                H(H<MASK(1))=values(1);
                H(H>MASK(2))=values(2);
            else
                H(H<MASK(1) & H>MASK(2)) = C(H<MASK(1) & H>MASK(2));
            end
        else
            if isscalar(values)
                H(~MASK)=values;
            elseif numel(values==2)
                H(H<MASK)=min(values(1));
                H(H>MASK)=max(values(2));
            else
                if all(size(H)==size(MASK)) && all(size(H)==size(values))
                    H(~MASK)=values(~MASK);
                else
                    return;
                end
            end
        end
    end

