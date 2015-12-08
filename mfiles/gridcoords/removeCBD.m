function [Model,grNew] = removeCBD(Model)
%REMOVECBD turn all confining beds of model into model layers
%
% Example;
%    [Model,grNew] = removeCBD(Model)
%
% Model is an array of class modelObj, grNew a gridObj
%
% Narrative:
%    Removes all confining beds in modelObj Model turning them into computed
%    model layers, thereby replacing VKCB by HK and VK resp.
%    Model is an array of modelObj, one of which must be the grid.
% 
%    The new model layers get HK and VK parameters respectively
%
%    These function are superseded by methods of the gridObj
%
% SEE ALSO: gridObj modelObj cutoutBCN cutoutXTA RefineBCN RefineGrid JoinBCN JoinLayers removeCBD
%
%  TO 100601 100610

grOld = Model{strmatchi('gr',Model(:,1)),2};

%% To remove all CBD, set LAYCBDnew equal to all zeros accoring to:
LAYCBDnew = zeros(length(grOld.LAYCBD)+sum(grOld.LAYCBD>0),1);

grNew = gridObj(grOld.xGr,grOld.yGr,grOld.Z,LAYCBDnew,grOld.MINDZ,grOld.AXIAL);

Intermediate = NaN(grOld.Ny,grOld.Nx,grOld.Nlay+grOld.Ncbd);

[L2N, LC2N, CC2N, N2L, N2C] = convertLAYCBD(grOld.LAYCBD,LAYCBDnew);


for i=1:size(Model,1)
    % switch Model Variable Type
    switch Model{i,3}
        case 'zlist'
            Inter(L2N( :,2),1) = Model{i,2}(L2N( :,1)); %#ok<AGROW>
            Inter(LC2N(:,2),1) = Model{i,2}(LC2N(:,1)); %#ok<AGROW>
            Model{i,2} = Inter(N2L(:,1));
        case '3Dtime',
            % Do nothing;
        case '3Dcbd'
            % Do nothgin;
        case '3Dlay'
            Intermediate(:,:,L2N(:,2))=Model{i,2}(:,:,L2N(:,1));
            switch Model{i,1}
                case {'HK','VK'}
                    ivkcb = strmatchi('VKCB',Model(:,1));
                    Intermediate(:,:,CC2N(:,2))=Model{ivkcb,2}(:,:,CC2N(:,1));
                    Model{i,    2} = Intermediate(:,:,N2L(:,1));
                    Model{ivkcb,2} = Intermediate(:,:,N2C(:,1));
                case {'TRAN','HY'}
                    ivkcb = strmatchi('VCONT',Model(:,1));
                    Intermediate(:,:,CC2N(:,2))=Model{ivkcb,2}(:,:,CC2N(:,1));
                    Model{i,    2} = Intermediate(:,:,N2L(:,1));
                    Model{ivkcb,2} = Intermediate(:,:,N2C(:,1));
                otherwise
                Intermediate(:,:,LC2N(:,2))=Model{i,2}(:,:,LC2N(:,1));
            end
            Model{i,2} = Intermediate(:,:,N2L(:,1));            
        case 'struct'
            if strmatchi('values',fieldnames(Model{i,2}))
                for iper=1:length(Model{i,2})
                    Intermediate(:,:,L2N( :,2))=Model{i,2}(iper).values(:,:,L2N( :,1));
                    Intermediate(:,:,LC2N(:,2))=Model{i,2}(iper).values(:,:,LC2N(:,1));
                    Model{i,2}(iper).values = Intermediate(:,:,N2L(:,1));
                end
            elseif strmatchi('term',fieldnames(Model{i,2}))
                for iper=1:length(Model{i,2})
                    for it=1:length(Model{i,2}(iper).term)
                        Intermediate(:,:,L2N( :,2))=Model{i,2}(iper).term{it}(:,:,L2N( :,1));
                        Intermediate(:,:,LC2N(:,2))=Model{i,2}(iper).term{it}(:,:,LC2N(:,1));
                        Model{i,2}(iper).term{it} = Intermediate(:,:,N2L(:,1));
                    end
                end
            else
                error('%s: Struct with unknown field <<%s>>.',mfilename,Model{i,3});
            end
        case 'stress'
            %% This works only if LAYCBDnew is all zeros (as is the case for this function)
            Model{i,2}(:,2) = N2L(L2N(Model{i,2}(:,2),2),1);
        case 'wellObj'
            well = Model{i,2};
            for iw = 1:length(well)
                well(iw).iLay     = N2L(L2N(well(iw).iLay,2),1);
                well(iw).LRC(:,1) = well(iw).iLay;
                well(iw).idx      = cellIndex(well(iw).ix,well(iw).iy,well(iw).iLay,grOld.size);
            end
            Model{i,2}=well;
        case 'wellSeriesObj'
            % Do Nothing
        case 'gridObj'
            Model{i,2} = grNew;
        otherwise
            error('mfLab:removeCBD:unknownVariableType',...
                '%s: Unknown type field <<%s>>',mfilename,Model{i,3});
    end
end

for i=size(Model,1):-1:1
    if strcmpi(Model{i,1},'VKCB') || strcmpi(Model{i,1},'VCONT');
        Model(i,:)=[];
    end
end

%% Consider VKA as vertical anisotropy if it exists

iVKA= strmatchi('VKA',Model(:,1));

if iVKA
    VKA = Model{iVKA,2};

    % Replace VKA in Model by VK
    Model{iVKA,1} = 'VK';
    % Replace by HK as vk is not given
    Model{iVKA,2} = Model{strmatchi('HK',Model(:,1)),2};
    % Set its new type
    Model{iVKA,3} = '3Dlay';

    % Consider VKA vertical anisotropy
    iVK=iVKA;
    for i=1:length(VKA)
        % Multiply VK with this anisotropy VKA
        Model{iVK,2}(:,:,i)=Model{iVK,2}(:,:,i)*VKA(i);
    end
end
