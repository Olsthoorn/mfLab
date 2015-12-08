classdef cObsObj < observationObj
    % cObsObj -- definition of class cObsObj by calling
    % obs = cObsObj(basename,obsSheetNm) where
    %  basename is the workbook holding obsSheetNm
    %  obsSheetNm is the worksheet with a table defining the observations. See
    %  constructor by typing
    % help cObsObj.cObsObj.
    % cObsObj defines observations for use in groundwater models. The objects hold
    % property values and share methods that facilitate using observations.
    % Some methods need a grid object. But a grid object is not necessary
    % to create obs objects. This is because upon creation, obs objects
    % only contain their pertinent data. Later on, when a grid object is
    % available they may be put into the grid, upon which additional grid
    % values are added to the observations, so that each obs then also holds its
    % position in the grid. See the methods of cObsObj and of gridObj for
    % more information.
    %
    % SEE ALSO gridObj wellObj
    %
    % TO 120807 130318
    properties
        C
        NCOMP
        compound
    end
    methods
        function o=cObsObj(varargin)
            % generate observation objects
            % USAGE:
            %       obs = cObsObj(xyz[,gr[,compoundNames,['use',criterion]]])
            %       obs = cObsObj(nr,x,y,z[,gr[,compountNames,['use',criterion]]])
            %       obs = cObjObj(basename,sheetNm[,gr[,tracerNames[,'use',criterion')
            %
            %         basename   = name of Excel workbook
            %         obsSheetNm = name of sheet holding the table of observations
            %         gr         = grid object holding current finite difference grid
            %         compoments are the names in a cell array of the
            %                    compoments to be read. 1:ncompoments are read
            %                    from files MT3D001.UCN .. MT3D00i.UCN
            %         use        = obs numbers or regexp for obs names to
            %                      select the ones to include
            %
            % SEE ALSO: gridObj wellObj
            %
            % TO 120807 130318
            o = o@observationObj(varargin{:}); % call observationObj constructor

            if nargin<1, return; end
            
            % remove numeric arguments
            I = cellfun(@isnumeric,varargin);
            if I(1), varargin(I)=[]; end
            
%             % remove grid
%             gr = [];
%             for i=numel(varargin):-1:1,
%                 if strcmp(class(varargin{i}),'gridObj');
%                     gr = varargin(i);
%                     varargin(i)=[];
%                     break;
%                 end
%             end
            
            % remove use
            for i=numel(varargin):-1:1
                if ischar(varargin{i}) && strcmpi(varargin{i},'use'), varargin(i:i+1)=[]; end
            end
            
            % get compound names and make sure they are in a cell array
            if numel(varargin)>2
                if ischar(varargin{3})
                    compoundNames = varargin(3);
                else
                    compoundNames = varargin{3};
                end
            elseif numel(varargin)>0
                if ischar(varargin{1})
                    compoundNames = varargin(1);
                else
                    compoundNames = varargin{1};
                end
            else
                compoundNames = 'concentration';
            end
            
            if ~isempty(gr) % then o contains gridInfo
                o = o.getConc(compoundNames);
            end
        end
                
        function plot(o,iComp,varargin)
            % cObsObj.plotCin(iComop,plotOptions); plots input
            % concentration of component iComp and according to plotOptions
            % 'color','k' etc.
            if nargin<2 || ~exist('iComp','var') || isempty(iComp)
                iComp=1;
            end
            for iw=1:length(o)
                if isempty(varargin)
                    plot(o(iw).t,o(iw).C(iComp,:),mf_color(iw));
                else
                    plot(o(iw).t,o(iw).C(iComp,:),varargin{:});
                end
            end
        end
    end
end
