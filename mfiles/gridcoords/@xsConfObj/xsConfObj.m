classdef xsConfObj
% MF_CONF: Gets grid values of the parameter parnam from specified zones
%   the zones are specified in the
%   worksheets "config" and "material" of Excel file [basename '.xls']
%   Used for convenient and efficient generation of cross sections.
%   See the examples that use it, e.g.
%        mflab>examples>swt_v4>freshkeeper>
%
% USAGE:
%   xsConfObj=xsConfObj(basename,configSheetNm,materialSheetNm)
%   see example VechtDredging under examples/mf2005/DutchTop/VechtDredging
%   for usage and its html file in the htlm directory
%   The data below originate from the VechtDredging example.
%
% Layout of config worksheet:
%  line 1: label label followed by arbitrary zone names as off column 3
%  line 2: xZone,         followed by xL of left-most zone. Ends with xR of all zones
% or
%  line 2: width|breedte, followed by left-most coordinate and the width of each zone
%  line 3: label, blank, ztop of all zones
%  line 4: label, blank, prescribed head of all zones, NaN is empty cell,
%    meaning head is not prescribed and will be calculated as is the case under dikes
%    adjacent to rivers and canals.
%  line 5: layer#, layerBotElev, character code id for all zones (codes may have more than one letter)
%  line 6: same
%  line 7: same etc
%  See workbook for exact alignment of columns.
%  NaN stands for empty cells in the spreadheet.
%
%Example:
% Layer   Basis(m) Polder2W DikeARKW ARK  DikeARKE PolderW Polder1W DitchW DikeW VechtW Dredged VechtE DikeE DitchE PolderE Polder2E Polder3E
% Width   -4297    2500     24       100  24       1500    100      2      12    20     30      20     12    2     800   1500  2500
% Head       NaN   -2.1     NaN      -0.4 NaN      -2.02   -2.02   -2.02   NaN -0.4    -0.4   -0.4     NaN  -1.86  -1.86 -1.71 -1.71
% Top (m)    NaN   -1.4     1.2      -0.4 1.2      -1.2    -1.2    -2.02   0.6  -4.25  -4.25   -4.25   0.5  -1.86  -1.3  -1.3  -1.7
% 1    -2      K    K    W    K    K    K    W    K    W    W    W    K    W    K    K    W
% 2    -2.5    K    K    W    K    K    K    K    K    W    W    W    K    K    V    V    W
% 3    -3      K    K    W    K    Z    Z    Z    K    W    W    W    K    K    V    V    W
% 4    -3.5    K    K    W    K    Z    Z    Z    K    W    W    W    K    K    V    V    W
% 5    -4      K    K    W    K    Z    Z    Z    K    W    W    W    K    K    V    V    W
% 6    -4.25   K    K    W    K    K    K    Z    K    W    W    W    K    K    V    V    W
% 7    -4.5    K    K    W    K    K    K    Z    Z    S    W    S    K    K    V    V    W
% 8    -5      K    K    K    K    K    K    Z    Z    Z    S    Z    K    Z    V    V    W
% 9    -6      K    K    K    K    K    K    K    K    Z    Z    Z    Z    Z    V    V    W
% 10   -6.5    V    V    V    V    V    V    V    V    Z    Z    Z    Z    Z    V    V    W
% 11   -10     P    P    P    P    P    P    P    P    P    P    P    P    P    P    P    B
% 12   -15     P    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P
% 13   -20     P    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P
% 14   -50     P    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P
% 15   -64     P    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P
% 16   -65     P    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P
% 17   -100    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P
% 18   -150    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P
% 19   -175    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P    P%
%
%
% The layout of the material worksheet
% line 1: headers of material properties
% line 2: dimensions of the materials properties
% line 3: values of material nr 1
% line 4: same for material nr 2
% line 5: etc.
%
% headers in line 1 are arbitrary except for 'material' and 'code'.
% material is the full name of the material and code its short name which
% must correspond to the mateials in the configuration layout specified
% above. This short name may be more than 1 character long.
%
% Example:
% Material                    Code    kh    kv    Red  Green  Blue  Porosity  rhoSolid  rhoDry   rhoWet
% dimension                    dim    [m/d] [m/d] [-]  [-]    [-]   [-]       [kg/m3]   [kg/m3]  [kg/m3]
% Pleistocene sand               P    30    10    1    1      0.5   0.35      2640      1716     2066
% Holocene sand                  Z    1     1     1    0.75   0     0.35      2640      1716     2066
% Sludge                         S    0.02  0.02  0.2  0.8    0     0.5       1400      700      1200
% Clay                           K    0.01  0.01  0    1      1     0.6       2100      840      1440
% Holland peat                   V    0.04  0.04  0.5  0.5    0     0.6       1400      560      1160
% Water                          W    1e6   1e6   0    0      0.8   1         1000      0        1000
% Water bottom                   B    0.01  0.01  0.75 0.5    0     0.6       2100      840      1440
% Clay intersected by ditches    D    0.1   0.1   0    1      1     0.6       2000      800      1400
%
% To prevent missing some layers, make sure that the top and bottom
% elevations of the zones match with vertical grid lines. You can do that
% by merging the two sets of elevations when defining zGr. Just put the two
% elevation sets in one vector and use that in gr=gridObj(xGr,yGr,zGr)
% gridObj guarantees correct merging.
%
% SEE ALSO: mf_conf, mf_zone mf_plotConf as non-object versions of xsConfObj.
%           and especially mfLab/examples/mf2005/DutchTop/VechtDredging/html
%           see the methods when viewing the object to see its options.
%
% TO 110319

    properties
        matProps      % material property headers
        matPropValues % material property matrix (numeric)
        matNames      % long names of materials
        matDims       % dimension of material properties
        matCodes      % short names of materials (= config codes)
        
        zoneNames  % name of zones
        zoneLeft;  % x-coord of left  end of zone
        zoneRight; % x-coord of right end of zone
        zoneHeads, % given zone head
        zoneTops,  % zone top elevation
        zoneBot,   % zone bottom elevation
        zoneZ      % zone Z values
        kD,        % total zone transmissivty
        c,         % total zone resistance
        dzMin = 0.1; % default minimum zone thickness (will be maintained)
        zoneCodes, % charater code of zones in config list
        zoneImat,  % index of each char zone in material list        
        Nz         % nr of layes per zone (= length zoneBot)
        NZone,     % nr of zones
		UserData,
    end
    properties (Dependent = true)
        xL, xR, xGr, xm, dx, zoneWidths,

    end
    methods
        % Constructor:
        function o=xsConfObj(basename,configSheetnm,materialSheetnm)
            % xsConf = xsConfObj(basename,configSheetnm,materialSheetnm) -
            % genererates xsConfObj (xs Configuration Objects)
            % basename = as ususal (name of workbook)
            % configSheetnm = sheetName in workbook holding the configuration
            % materialSheetnm = sheetName in workbook holding the Materials definitions table
            % See help of xsConfObj for the layout of these tables.
            % TO 1200731
            
            if nargin==0,
                return;
            end

            %% Get material
            try
                warning('off'); %#ok
                [~,txt] = xlsread(basename,materialSheetnm,'','basic');
                warning('on'); %#ok
            catch ME
                error('%s: %s\nCan''t find sheet <<%s>> in workbook <<%s>>.',...
                    mfilename,ME.message,materialSheetnm,basename);
            end
            
            txtHdr     = txt(1,1:2);
            o.matProps = txt(1,3:end);
            o.matDims  = txt(2,3:end);

            [~,o.matPropValues,~,txtValues]=getExcelData(basename,materialSheetnm,'Hor','fill',false);

            o.matNames = txtValues(:,strmatchi('Mat' ,txtHdr));
            o.matCodes = txtValues(:,strmatchi('Code',txtHdr));

            %% Get zones
            try
             [~,~,Raw]=xlsread(basename,configSheetnm,'','basic');
            catch ME
                error('%s: %s\nCan''t find sheet <<%s>> in workbook <<%s>>.',...
                    mfilename,ME.message,configSheetnm,basename);
            end
             I = ~(cellfun(@isempty,Raw) | cellfun(@ismynan,Raw));

             %% Remove empty line (shift data to left top corner)
             datarow = any(I,2);
             datacol = any(I,1);
             Raw = Raw(datarow,datacol);
            
            %% There is a bit of hard wiring of positions here
            if any(strmatchi({'xZone','yZone'},Raw(:,1)))
                o.zoneLeft  = [Raw{2,2:end-1}]; % left  zone coords
                o.zoneRight = [Raw{2,3:end  }];   % right zone coords
            elseif strmatchi('breedte',Raw(:,1),'exact') || strmatchi('width',Raw(:,1),'exaxt')        
                o.zoneLeft  = Raw{2,2} + [0 cumsum([Raw{2,3:end-1}])];
                o.zoneRight = Raw{2,2} + cumsum([Raw{2,3:end}  ]);
            else
                error('%s: Cell A2 in sheet %s must either be ''%s'' or ''%s''',mfilename,configSheetnm,'width','xzone');
            end            
            
            o.zoneHeads = [Raw{3,3:end}];
            o.zoneTops  = [Raw{4,3:end}];
            o.zoneBot   = [Raw{5:end,2}]';
            o.zoneZ     = [o.zoneTops ; o.zoneBot * ones(size(o.zoneTops))];
            o.NZone     = length(o.zoneTops);
            o.Nz        = length(o.zoneBot);
            
            %% make sure all zone bottoms are lower than zone tops
            for iL=2:size(o.zoneZ,1)
                o.zoneZ(iL,:)=min(o.zoneZ(iL,:),o.zoneZ(iL-1,:)-o.dzMin);
            end

            %% getting the character zoneCodes separately:
            o.zoneNames = Raw(1,3:end);
            
            if numel(unique(o.zoneNames))~=numel(o.zoneNames)
                error('%s: zoneNames must be unique in sheet <<%s>> workbook <<%s>>\n',...
                    mfilename,configSheetnm,basename);
            end
            
            o.zoneCodes = Raw(5:end,3:end);

            %% replace zoneCodes by direct index into material list
            o.zoneImat = NaN(size(o.zoneCodes));

            for ic=1:numel(o.matCodes)
                I = strmatchi(o.matCodes{ic},o.zoneCodes,'exact');
                if I(1)
                    o.zoneImat(I) = ic;
                end
            end
            if any(o.zoneImat==0)
                fprintf('Indices in Material list of zone codes (must be <>0):\n');
                display(o.zoneImat);
                error('Not all zone codes in the configuration are defined in the material list;\nSee index table just printed.');
            end
            
            kh_  = o.matPropValues(:,strmatchi('kh',o.matProps));
            kv_  = o.matPropValues(:,strmatchi('kv',o.matProps));

            o.kD = reshape(kh_(o.zoneImat),size(o.zoneImat)).*abs(diff(o.zoneZ));
            o.c  = abs(diff(o.zoneZ)) ./ reshape(kv_(o.zoneImat),size(o.zoneImat));

        end
        
        function A = zonePropArray(o,PropName)
            % A = ConfObj.confArray(PropName) - get an array of Properties for the zones of the configuration
            % PropName is the name of one of the material properties from
            % the materials list (under header Code (see xsConfObj.matProps))            
            iProp = strmatchi(PropName,o.matProps,'exact');
            if ~iProp,
                error('xsConfObj:confArray:unknownPropname',...
                    'can''t find propName <<%s>>, use one of <<%s>>',o.matProps);
            end
            A = reshape(o.matPropValues(o.zoneImat,iProp),size(o.zoneImat));
        end
        
        function xL  = get.xL(o)
            % xL = xsConf.xL() - gets left coordinate of zones
            xL = o.zoneLeft;
        end
        
        function xR  = get.xR(o)
            % x = xConf.xR()  - gets right coordinate of zones
            xR = o.zoneRight;
        end
        
        function xGr = get.xGr(o)
            % xGr = xsConf.xGr() - gets the grid lines of the zones
            xGr = [o.zoneLeft(1) o.zoneRight];
        end
        
        function  xm = get.xm(o)
            % xm = xsConf.xm() - gets the centers of the zones
            xm = 0.5*(o.xGr(1:end-1)+o.xGr(2:end));
        end
        
        function zoneWidths = get.zoneWidths(o)
            % zoneWidths = xsConf.zoneWidths -- get widths of zones
            zoneWidths = o.zoneRight - o.zoneLeft;
        end

        function  dx = get.dx(o) % abbreviation for convenience
            % dx = xsConf.dx() - gets widths of zones
            dx = o.zoneWidths;
        end
        
        function A = array2D(o,propName,xGr,Z,code)
            % A = xsConf.array2D(propName,xGr [Z [code]])
            % - generates a full model array for property propName
            % propName ==>
            %         'Z'     ==> generates Z array implied by xsConf.zoneTop, xsConf.zoneBot and xGr
            %         'head'  ==> top head extended to full array
            %         any propName in Conf.matNames
            % xGr ==> is the xGr of the model array not of the zones
            % Z   ==> the Z   of the model array (Nz+1,Nx), 2D because
            %         xsConf is 2D for cross sections with Z vertical downward
            %         (1st dimension) and X horizontal (2nd dimension)
            % if Z is omitted or empty, then Z is assumed to be implied by xsConf.zoneTop and xsConf.zoneBot
            %         i.e. that the layers in xsConfObj are the same as those of
            %         the model.
            % code ==> can be 'geometric', 'harmonic', see help gridsTransfer()
            %         if code is omitted 'geometric' is used.
            % A has dimensions (Nz,Nx) use A=XS(A) to make it 3D.
            %
            % TO 120803
            
            dim = 1;
            Nx=length(xGr)-1;
            xm = 0.5*(xGr(1:end-1)+xGr(2:end));
            
            if ~exist('code','var'), code='geometric'; end

            Idx = floor(interp1(o.xGr,1:o.NZone+1,xm));
            
            for iZone=o.NZone:-1:1
                Ix{iZone}=find(xm>=o.xL(iZone) & xm<=o.xR(iZone)); % ix indices model in this conf zone
            end

            if ~exist('Z','var') || isempty(Z) % use Z implied by config
                Nz_= o.Nz;
                Z  = NaN(o.Nz+1,Nx);
                for iZone=1:o.NZone
                    Z(:,Ix{iZone}) = o.zoneZ(:,iZone) * ones(size(Ix{iZone}));
                end
            else  % Z is given in its entirety by the call

                Nz_ = size(Z,1)-1;
            end
            
            for iZone=o.NZone:-1:1
                layersareuniform(iZone) = all(abs(mean(Z(:,Ix{iZone}),2)-Z(:,Ix{iZone}(1)))<1e-1);
            end

            lpropName = lower(propName);
            switch lpropName
                case 'z'
                    A = Z(:,Idx);
                case 'top' % full array of heads, using o.zoneHeads interpolating when isnan(zoneHeads)
                    A = blockInterp(o.zoneTops,o.xGr,xm);
                case 'tophead' % full array of heads, using o.zoneHeads interpolating when isnan(zoneHeads)
                    A = NaN(1,Nx);
                    for iZone= 1:o.NZone
                        A(xm>o.xL(iZone) & xm<o.xR(iZone)) = o.zoneHeads(iZone);
                    end
                case 'head' % full array of heads, using o.zoneHeads interpolating when isnan(zoneHeads)
                    A = ones(Nz_,1)*blockInterp(o.zoneHeads,o.xGr,xm);
                case 'ibound'
                    A = ones(Nz_,Nx);
                case {'matindex','mat'}
                    A = NaN(Nz_,Nx);
                    for iZone = 1:o.NZone
                        if layersareuniform(iZone)
                            A(:,Ix{iZone})   = gridsTransfer(o.zoneZ(:,iZone),o.zoneImat(:,iZone),Z(:,Ix{iZone}(1)),code,dim)*ones(1,length(Ix{iZone}));
                        else
                            for i=1:length(Ix{iZone})
                                A(:,Ix{iZone(i)}) = gridsTransfer(o.zoneZ(:,iZone),o.zoneImat(:,iZone),Z(:,Ix{iZone}(i)),code,dim);
                            end
                        end
                    end
                otherwise
                    iProp = strmatchi(propName,o.matProps,'exact');
                    if ~iProp,
                        error('xsConfObj:confArray:unknownPropname',...
                            ['can''t find propName <<%s>>, use one of <<', ...
                            repmat(' %s',[1, length(o.matProps)]), '>>'],propName,o.matProps{:});
                    end
                    A = NaN(Nz_,Nx);
                    for iZone=1:o.NZone
                        if layersareuniform
                            A(:,Ix{iZone})   = gridsTransfer(o.zoneZ(:,iZone),o.matPropValues(o.zoneImat(:,iZone),iProp),Z(:,Ix{iZone}(1)),code,dim)*ones(size(Ix{iZone}));
                        else
                            for i=1:length(Ix)
                                A(:,Ix{iZone}(i))   = gridsTransfer(o.zoneZ(:,iZone),o.matPropValues(o.zoneImat(:,iZone),iProp),Z(:,Ix{iZone}(i)),code,dim);
                            end
                        end
                    end
            end
        end
        
        function A = array3D(o,varargin)
            % xsConfOjb.model3D(gr,propName,code,dir)
            % gr is grid object
            % propName is the name of one of the material properties (under
            %    heading "Code" in the materials worksheet
            % Spcial property names are 'head'  'top' and 'bottom' to get
            %    the full 3D array of the heads, top and bottom of the model
            %    cells as derived from the zones. Note that only the top head
            %    of the zones is defined in the configuration. Hence all
            %    layers will have this head.
            % code is a code determining the way of averaging when the the
            %    zone and the model grid are intersected:
            % The codes:
            %   'k' or 'geometric' meaning geometric mean is used
            %   'c' or 'harmonic'  meaning harmonic  mean is used
            % Example:
            %   xsConfObj.getMdlArray(gr,'kh','k');
            %   xsConfObj.getMdlArray(gr,'kv','c');
            % where 'kv' and 'kh' are material properties and 'k' and 'c'
            % Note: any value outside the Conf layout will be a NaN
            % The output model is of size (1,Nx,Nz) i.e. 3D vertical sheet
            % TO 120731
            
            [gr,varargin] = getType(varargin,'gridObj',[]);
            if isempty(gr)
                error('gridObj missing in input');
            end
            
            [propName, varargin ] =getNext(varargin,'char','');
            if isempty(propName)
                error('%s: Second arguments must be a material code',mfilename);
                if ~strmatchi(propName,o.matNames)
                    error('%s: no such material code <<%s>> in xsConfObj',mfilename,propName);
                end
            end
            
            [code,varargin ] = getNext(varargin,'char','geometric');
            
            [dir, ~        ] = getNext(varargin,'char','x');
            
            if dir == 'x'
                A = XS( o.array2D(propName,gr.xGr,XS(gr.Z),code) );
            else
                A = o.array2D(propName,gr.yGr',YS(gr.Z),code);
                A = permute(A,[2 3 1]);
            end
        end
        
        function ax =plotOverview(o,varargin)
            % xsConfObj/plotOverview -- plots overview of cross section
            % methods of xsConf
            % USAGE:  xsConfObj.plotOverview(gr,H,B,ttl,leg,varargin)
            %   H(end).values(1,:,:) are contoured and mf_Psi(B(end))
            %   H is obtained from readDat and B from readBud
            %   ttl is struct{3,1} of titles for the three graphs
            %   leg is struct{3,1} of legends for the three graphs
            %   varargin is optional axis properties.
            %
            % TO 120101, 141209 (not yet finished)

            [gr,varargin]  = getType(varargin,'gridObj',[]);
            [H ,varargin]  = getNext(varargin,'struct' ,{});
            [B ,varargin]  = getNext(varargin,'struct' ,{});
            [ttl,varargin] = getNext(varargin,'cell', '');
            [leg,varargin]  = getNext(varargin,'cell', '');
            
            if isempty(gr), error('grid not specified'); end
            if isempty(H) , error('specify heads as H(i) from readDat'); end
            if isempty(B) , error('specify flow  as B(i) from readBud'); end
            if isempty(ttl) || numel(ttl)~=3, error('specify ttl as cell{3,1} of strings'); end
            if isempty(leg) error('specify legends in cell array of strings'); end
            
            % Compute net inflow of nodes
            iR = strmatchi('FLOWR',B(end).label);
            if iR, Qx = [0 sum(XS(B(end).term{iR}(end,:,:)),1)]; end
            
            Q  = diff(Qx,1,2);

            if isempty(varargin), varargin = {'visible','on'}; end

            args ={'nextplot','add','xgrid','on','ygrid','on'};
            args = [args varargin];
            xlbl = 'x [m]';
            ylbl = {'head [m]';'discharge [m2/d]';'seepage mm/d'};

            figure;

            for i=3:-1:1
                ax(i) = subplot(3,1,i,args{:});
                xlabel(xlbl); ylabel(ylbl{i}); title(ttl{i});
                switch i
                    case 1
                        plot(gr.xm,o.array2D('top',gr.xGr),'k','linewidth',0.5);
                        plot(gr.xm,XS(H(end).values(1,:,:)));                        
                    case 2
                        plot(gr.xGr,Qx);
                    case 3
                        plot(gr.xm,Q*1000./gr.dx); % mm/d inflow from top
                end
                if ~isempty(leg)
                    legend(leg{:});
                end
            end
            % make sure axes stay the same (especially for interactive zooming)
            linkaxes(ax,'x');
        end
        
        function ax = plotContours(o,varargin)
            % xsConfObj/plotContours --- contours cross section with both
            % the heads H(end).values and the stream function Psi obtained
            % from B(end), where H is obtained from readDat and B from
            % readBud. varargin is other options accepted by contour.
            %
            % USAGE ax = xsConfObj.plotContours(gr,H,B,axisOptions)
            %
            % TO 141209
            
            [gr,varargin] = getType(varargin,'gridObj',[]);
            [H ,varargin] = getNext(varargin,'struct',[]);
            [B ,varargin] = getNext(varargin,'struct',[]);
            
            if isempty(gr), error('grid argument missing'); end
            if isempty(H ), error('missing heads (use H(i) from readDat'); end
            if isempty(B ), error('missing flows (use B(i) from readBud'); end
                
            B = mf_Psi(B(end));
            %% Plot stream-line patterns with heads and flows

            args = {'nextplot','add','xgrid','on','ygrid','on'};
            args = [args varargin];


            ax(1) = o.paint(args{:}); % plot background

            ax(2)= axes('position',get(ax(1),'position'),'color','none',...
                        'xlim',get(ax(1),'xlim'),'ylim',get(ax(1),'ylim'),args{:});

            xlabel('x [m]');  ylabel('z [m NAP]');

            linkaxes(ax);

            hrange =ContourRange(H.values,50); dPhi = diff(hrange(1:2));
            prange =ContourRange(B.Psi   ,50); dPsi = diff(hrange(1:2));

            title(sprintf('Head and stream lines, dPhi=%.2f m, dPsi=%.2f m2/d',dPhi,dPsi));

            contour(ax(2),gr.xc,gr.zc,XS(H(end).values),hrange,'g');
            contour(ax(2),gr.xp,gr.zp,B(end).Psi       ,prange,'b');

            plot(ax(2),gr.xm,o.array2D('top',gr.xGr),'k'); % top of model
            plot(ax(2),gr.xm,H(end).values(1,:,1),'b'); % water table

            legend('head lines','stream lines','ground surface','prescribed head')
        end
         
        function ax = paint(o,varargin)
            % xsConfObj.paint() - paints the zones each iwth unique color
            % if materialsTable contains columns 'Red', 'Green' and 'Blue'
            % these will be used, else mf_color is used with its default
            % sequential standard colors (try mf_color(1:10));
            %
            % TO 120731

            figure;
            args = [{'nextplot','add','xgrid','on','ygrid','on'},varargin];
            
            ax = axes(args{:});
            
            iR = strmatchi('Red'  ,o.matProps);
            iG = strmatchi('Green',o.matProps);
            iB = strmatchi('Blue' ,o.matProps);
            if all([iR iG iB]>0)
                propColor = [o.matPropValues(:,strmatchi('Red'  ,o.matProps)), ...
                             o.matPropValues(:,strmatchi('Green',o.matProps)), ...
                             o.matPropValues(:,strmatchi('Blue' ,o.matProps))];
            else
                propColor = mf_color(1:max(o.zoneImat(:)));
            end
            for i=1:size(o.zoneImat,2)
                zGr = [max(o.zoneTops(i),o.zoneHeads(i)); o.zoneBot];
                for j=1:size(o.zoneImat,1)
                    fill(o.xGr([i i+1 i+1 i  ]),...
                           zGr([j j   j+1 j+1]),...
                           propColor(o.zoneImat(j,i),:),'edgecolor','none','parent',ax);
                end
            end
        end
        
        function show(o,what)
            if ~exist('what','var') || isempty(what), what='c'; end
            
            w=lower(what(1));
            switch w
                case {'m','materials'}
                    fprintf('Overview materials in this model:\n');
                    fprintf('Material properties:      ');
                    for im=1:length(o.matProps), fprintf('%12s',o.matProps{im}); end; fprintf('\n');
                    fprintf('Property dimensions:      ');
                    for im=1:length(o.matDims) , fprintf('%12s',o.matDims{im});  end; fprintf('\n');
                    for i=1:length(o.matNames)
                        fprintf('%20s%6s',o.matNames{i},o.matCodes{i});
                        for j=1:size(o.matPropValues,2)
                            fprintf('%12g',o.matPropValues(i,j));
                        end
                        fprintf('\n');
                    end
                    fprintf('\n');
                case {'c','config'}
                    fprintf('\nConfiguration:\n');
                    fprintf('zone xL [m]:                          ');  fprintf(' %10.1f',o.xL); fprintf('\n');
                    fprintf('zone xR [m]:                          ');  fprintf(' %10.1f',o.xR); fprintf('\n');
                    fprintf('zone width [m]:                       ');  fprintf(' %10.1f',o.zoneWidths); fprintf('\n');
                    fprintf('top (ground surface [m NAP])          ');  fprintf(' %10.2f',o.zoneTops);   fprintf('\n');
                    fprintf('prescribed head peil: [m NAP]:        ');  fprintf(' %10.2f',o.zoneHeads);  fprintf('\n');
                    fprintf('Totaal transmissivity: [m2/d]:        ');  fprintf(' %10.4g',sum(o.kD,1));  fprintf('\n');
                    fprintf('Totale vert. hydraul. resistance [d]: ');  fprintf(' %10.4g',sum(o.c ,1));  fprintf('\n');

                    for iL=1:o.Nz
                       fprintf('Ground layer %2d [m NAP]: %10.0f:',iL,o.zoneBot(iL));
                       fprintf(' %3s',o.zoneCodes{iL,:}); fprintf('\n');
                    end
                    fprintf('\n');

                otherwise
                    error('xsConfObj:display:unknownOption',...
                        '%s: unknown option <<%s>>',mfilename,what)
            end
        end
end

end

