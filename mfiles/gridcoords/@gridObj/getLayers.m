function [LAYnams,LAYvals]=getLayers(XLSF,gr)
% [LAYnams,LAYvals]=gr.getLayers(XLSF,gr)
% XLSF is workbook for this model with sheet LAY
%
% TO 091212 091224 120407 121122

if nargin<2,
    error('Not enough input arguments in getLayers');
end

[LAYnams,layvals]=getExcelData(XLSF,'LAY','Horizontal');

iNr=strmatchi('LAYER',LAYnams); % column with layer number

layvals(layvals(:,iNr)> gr.Nlay,:)=[]; % make sure numbering is compatible with Z
layvals(layvals(:,iNr)<= 0,:)=[]; % only layers with Layer>0 are valid

if isempty(layvals)
    error('mf_setup:getLayers:noLayers',...
        ['No layers, check LAY sheet to see that any layers>0 and <NLay of model\n',...
        'At least one layer with number 1 must be present in LAY worksheet.\n']);
end

% at least the first and the last layer are required
layvals=sortrows(layvals,iNr);
if size(layvals,1)==1
    layvals=[layvals;layvals];
    layvals(end,1)=gr.Nlay;
end
layvals(end,iNr)=gr.Nlay;             % sure to number first and last layer
layvals(  1,iNr)=1;                % correctly


% Fill specified layers in the overall layer matrix, gaps are permitted
LAYvals=NaN(gr.Nlay,size(layvals,2)); % Allocate memory to hold all layers
    
iL=NaN(size(layvals,1),2);
iL( :     ,2) = layvals(:,iNr); iL(:,1)=1;
if size(layvals,1)>1, iL(2:end,1)=iL(1:end-1,2)+1; end

for iP=1:size(iL,1)       % fill in the specified layers, backward
    for i=iL(iP,1):iL(iP,2)
            LAYvals(i,:)=layvals(iP,:); % insert given layers
    end
end

% Check and correct LAYCBD
iCBD=strmatchi('LAYCBD',LAYnams);
% if ~exist('LAYCBD','var'),
%     LAYvals(end,iCBD)=0;
%     LAYCBD=LAYvals(:,iCBD);
% end
% [~,LAYCBD]=isAquifer(length(LAYCBD)+sum(LAYCBD>0),LAYCBD);

LAYvals(:,iCBD)=gr.LAYCBD;



