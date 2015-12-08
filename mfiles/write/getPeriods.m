function [PERnams,PERvals,NPER]=getPeriods(XLSF,varargin)
%GETPERIODS reads stress period information from workhsheet and expands when necessary
%
% Example / Usage
%       [PERnams,PERvals,NPER] = getPeriods(XLSF[,sheetNm[,Hdrs]])
%       [PERnams,PERvals,NPER] = getPeriods(XLSF)  --- default sheeNm 'PER'
%       [PERnams,PERvals,NPER] = getPeriods(XLSF,sheetNm)
%       [PERnams,PERvals,NPER] = getPeriods(XLSF,sheetNm,Hdrs)
%             only columns with header names given in cell array Hdrs
%       NPER                   = getPeriods(XLSF)  --- only get NPER
%
% What it does:
%    Fetches the column labels PERnams and data PERvals from the worksheet PER
%    in the parameter workbook XLSF
%    Then completes the stress period list from the top down
%    where the highest encountered period number is NPER by definition.
%
% TO 091219 100120 130309

[fill, varargin ] = getProp(varargin,'fill',true);

[sheetNm,varargin] = getNext(varargin,'char','PER');
[Hdrs,   varargin] = getNext(varargin,'cell',{});  %#ok

[PERnams,pervals]=getExcelData(XLSF,sheetNm,'Horizontal','fill',fill);

% make sure that the list of stress periods is consistent, apply rules
% for flexibility as defined in the manual
% TO 091211 091218

JP=strmatchi('IPER',PERnams);  % column which holds stress period number
pervals(:,JP)=round(pervals(:,JP)); % in case period numbers are computed in Excel

% Only the stress periods that have iper>0 are valid
pervals=pervals(pervals(:,JP)>0,:);
pervals=sortrows(pervals,JP);

NPER=pervals(end,JP);

%% if nargout==1, only NPER is requested !!
if nargout==1,
    PERnams = NPER;
    return;
end

% Allow gaps by way of short-hand input

PERvals=NaN(NPER,size(pervals,2)); % Allocate memory for all stress periods 
for iP=size(pervals,1):-1:1         % fill specified layers in, backward
    perNr=min(NPER,pervals(iP,JP));
    PERvals(perNr,:)=pervals(iP,:); % insert given stress periods
end
% filling in backward implies that if there are duplicates, the first value wins

% fill in the gaps
for iPer=NPER-1:-1:1
    if isnan(PERvals(iPer,1))
        PERvals(iPer,:)=PERvals(iPer+1,:);
    end
end

for iPer=1:NPER
    PERvals(iPer,JP) = iPer;
end

% if we only require some of the field
if nargin>2 % i.e. exist('Hdrs','var')
    if ~(iscell(Hdrs) || ischar(Hdrs))
        error('%s: 3rd argument (Hdrs) must be string or cell array of strings',...
            mfilename);
    end
    I = strmatchi(Hdrs,PERnams);
    if I(1)
        PERnams  = PERnams(I);
        PERvals  = PERvals(:,I);
    else
        % skip, we use all headers
        % error('%s: none of the headers { %s } found in worksheet <<%s>> of workbook <<%s>>',...
        %    mfilename,sprintfs(' %s',Hdrs),XLSF);
    end
end
