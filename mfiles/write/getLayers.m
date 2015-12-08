function [LAYnams,LAYvals]=getLayers(varargin)
%GETLAYERS reads layer information from LAY worksheet, expands where necessary
%
% Example:
%    [LAYnams,layVals]=getLayers(XLSF,gr|'Nlay',Nlay[,'LAYCBD',LAYCBD])
%
% INPUT:
%    XLSF      workbook file name for this model with sheet LAY
%    gr        gridObj given, which yields the number of layers
%    'Nlay',Nlay     propName,propValue pair to specify number of layyrs
%    'LAYCBD',LAYCBD propName,propValye pair to specify LAYCBD, whose length
%                  is the number of layers.
% OUTPUT
%    LAYnams headers of columns in LAY worksheet
%    LAYvals values in LAY worksheet
%
% if the second argument is a gridObj, then LAYCBD will be obtained from
% the gridObj and not from the workhseet, also the third argument will be
% ignored.
%
% TO 091212 091224 120407 130307 130312

[XLSF,varargin] = getNext(varargin,'char','');
if isempty(XLSF)
    error('%s: first argument must be the basename of the Excel workbook',mfilename);
end

% The second argument in varargin must be either gr,Nlay or a LAYCBD vector
% Nlay can be specified using 'Nz',Nlay or 'Naly,Nlay propName, propValue
% pair and LAYCBD as 'LAYCBD',LAYCBD propName,propValue pair


% layVals comes from getExelData
% LAYvals is the full layer array to ouput
[LAYnams,layVals]=getExcelData(XLSF,'LAY','Horizontal');

layNrCol = strmatchi('LAYER',LAYnams); % column with layer number

[gr,varargin] = getType(varargin,'gridObj',[]);
if ~isempty(gr)
    Nlay   = gr.Nlay;
    LAYCBD = gr.LAYCBD;
else
    [Nlay,varargin] = getProp(varargin,{'Nz','Nlay'},[]);
    if ~isempty(Nlay)
        LAYCBD = 0;
    else
        [LAYCBD,varargin] = getProp(varargin,'LAYCBD',[]);
        if isempty(LAYCBD)
            Nlay = max(layVals(strmatchi('Nr',LAYnams)));
            LAYCBD = zeros(Nlay,1);
        else
            Nlay = numel(LAYCBD);
            LAYCBD = LAYCBD(:);
        end
    end
end


if isempty(Nlay)
    Nlay = size(layVals,1);
end

layVals(layVals(:,layNrCol)> Nlay,:)=[]; % make sure numbering is compatible with Z
layVals(layVals(:,layNrCol)<= 0,:)=[]; % only layers with Layer>0 are valid

if isempty(layVals)
    error('mf_setup:getLayers:noLayers',...
        ['%s: No layers, check LAY sheet to see that any layers>0 and <Nlay of model\n',...
        'At least one layer with number 1 must be present in LAY worksheet.\n'],mfilename);
end

% at least the first and the last layer are required
layVals=sortrows(layVals,layNrCol);
layVals(end,layNrCol)=Nlay;             % ensure to number the last layer

LAYvals = NaN(Nlay,numel(LAYnams));

il = 1;
for iLay= 1:Nlay
    if iLay<=layVals(il,layNrCol)
        LAYvals(iLay,:) = layVals(il,:);
        LAYvals(iLay,layNrCol) = iLay;
        if iLay==layVals(il,layNrCol), il=il+1; end
    end
end

% Check and correct LAYCBD
 if ~exist('LAYCBD','var'),
     iCBD=strmatchi('LAYCBD',LAYnams);
     layVals(end,iCBD)=0;
     LAYCBD=layVals(:,iCBD);
 
    if ~all(layVals(:,iCBD)==LAYCBD)
        msgId ='mflab:LAYCBD:sheetVSgrid';
        warning('on',msgId);
        warning(msgID,'LAYCBD in worksheet LAY differs from that in grid');
        warning('off',msgId);
    end
 end
 
if ~isempty(varargin)
    warning('mfLab:getLayers:vararginNotCompletelyUsed',...
        'not all values of varargin have been used');
end