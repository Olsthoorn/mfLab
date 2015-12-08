%KMLFOLDEROBJ get the coordinates from a kml folder file with paths
% and at the same time compute its length and area. You may also directly
% compute the depth of the object interpolating from start to end given the
% depths at start and end. Finally you can print the object such that it
% can directly be copied into Excel.
% Not the kml file can contain a kmlpath or an entire directly with
% kmlpaths as obtained by copying a folder with such paths from Google
% Earth. The kmlfile thus obtained contains the paths of all the paths
% that were in the folder.
%
% SEE ALSO: kmlPath, kmlFolder, wgs2utm, utm2wgs gridObj
%
% USAGE:
%     p = kmlPathsObj(kmlPathsFileName[,gr])
%     p = kmlPathsObj(kmlPathsFileName,gr,depthToe,depthHeel);
%     p = kmlFilderObj(kmlFilderFileName,gr,depthToe,depthHeel,cellArray);
%     where cellArray defines the toe and heel depths below ground surface
%     for named khettaras or other buried objects:
%     cellArray = {name1 toeDepth1 heelDepth1; ...
%                  name2 toeDepth2 heelDepth2; ...
%                  ...
%                  };
%     It will apply these depths to specifally named khettaras if present
%     otherwise it uses the default
%
% TO 13106

classdef kmlPathsObj
    properties
        name
        E
        N
        X
        Y
        zDem
        z
        UserData
    end
    properties (Dependent=true)
        length
        perimeter
        area
        isClosed
    end
    methods
        function o = kmlPathsObj(varargin)
            %KMLPATHSOBJ constructor of a kmlPathsObj
            %
            % USAGE obj  = kmlFolderObj(kmlFolderName.kml)
            %
            %  kmlFolderName is name of kmlfile that contains a set of
            %  paths as obtained by copying a folder of paths from Google
            %  Earth to another directory.
            %
            % SEE ALSO kmlPathsObj kmlPathsObj
            %
            % TO 131010

            % get data from kml
            if nargin==0
                return;
            end
            
            [kmlPathsFileName,varargin] = getNext(varargin,'char',[]);
            if isempty(kmlPathsFileName)
                error('First argument must be the kmlFolderFleName');
            end
            [gr,varargin] = getType(varargin,'gridObj',[]);
                        
            LLstruct = kmlFolder(kmlPathsFileName);
            
%            o(numel(LLstruct,1))=kmlPathsObj();
            for io=numel(LLstruct):-1:1
                o(io).name = LLstruct(io).name;
                o(io).E = LLstruct(io).E;
                o(io).N = LLstruct(io).N;
                o(io).X = LLstruct(io).X;
                o(io).Y = LLstruct(io).Y;
                if ~isempty(gr)
                    o(io).getDem(gr);
                end
            end
            if ~isempty(varargin)
                o = o.setElevation(varargin);
            end
        end
        function o = setDem(o,gr)
            for io=numel(o):-1:1
                o(io).zDem = interp2(gr.xc,gr.yc,gr.Z(:,:,1),o(io).X,o(io).Y);
                I = isnan(o(io).zDem);
                if all(I)
                    o(io)=[];
                else
                    o(io).E(I)=[];
                    o(io).N(I)=[];
                    o(io).X(I)=[];
                    o(io).Y(I)=[];
                    o(io).zDem(I)=[];
                end
            end
        end
        function o = setElevation(o,varargin)
            
            defaultDepthToe      = 17;
            defaultDepthHeel     = 1.5;
            
            [depthToe ,varargin] = getNext(varargin,'double',defaultDepthToe);
            [depthHeel,varargin] = getNext(varargin,'double',defaultDepthHeel);
            [depthNames,      ~] = getNext(varargin,'cell',[]);
            
            names = depthNames(:,2);
            for io=numel(o):-1:1

                name_ = o(io).name;
                i = regexp(o(io).name,' ','once');
                if ~isempty(i)
                    name_ = name_(i+1:end);
                    name_(find(name_==' ',1,'first'):end)=[];
                end
                
                i = strmatchi(name_,names);
                if i
                    depthToe = depthNames{i,3};
                    depthHeel= depthNames{i,4};
                end

                s = cumsum([0; sqrt(diff(o(io).X).^2+diff(o(io).Y).^2)]);
                zHeel = o(io).zDem(end)-depthHeel;
                zToe  = min(o(io).zDem(1),...
                            max(o(io).zDem(  1)-depthToe,...
                                zHeel + s(end)/1000)...
                            );
                if numel(s)>1
                    o(io).z = interp1(s([1 end]),[zToe zHeel],s);
                else
                    o(io).z = zToe;
                end
            end
        end
        function plot(o,varargin)
            %KMLPATHSOBJ/PLOT -- plots the objects
            %
            % USAGE:  kmlPathsObj.plot(plotOptions)
            % SEE ALSO: lineObj area2Obj
            %
            % TO 131015
            
            [ax,varargin] = getType(varargin,'axis',gca);
            for io=numel(o):-1:1
                plot(ax,o(io).X,o(io).Y,varargin{:});
            end
        end
        function print(o)
            %PRINT: prints out object such that it can be copied directly
            %into excel
            % TO 1301006
            for io=1 % use first obj for settings
                fn = fieldnames(o(io));
                fn = fn(ismember(fn,{'E','N','X','Y','zDem','z'}));
                for i=numel(fn):-1:1
                    if isempty(o(io).(fn{i}))
                        fn(i)=[];
                    end
                end
            end
            fmt = repmat('%s\t',[1 numel(fn)]); fmt(end-1:end)='\n';
            fprintf('name\t'); fprintf(fmt,fn{:});
            for io=1:numel(o)
                for i=1:numel(o(io).E)
                    if i == 1
                        fprintf('%s',o(io).name);
                    else
                       % skip
                    end
                    for j=1:numel(fn)
                        fprintf('\t%.9g',o(io).(fn{j})(i));
                    end
                    fprintf('\n');
                end
                fprintf('\n');
            end
        end
        
        function isClosed = get.isClosed(o)
            n=6; % 6 significant digits
            isClosed = roundn(o.X(1),n)==roundn(o.X(end),n) && ...
                       roundn(o.Y(1),n)==roundn(o.Y(end),n);
        end
        function length = get.length(o)
                length = sum(sqrt(diff(o.X).^2+diff(o.Y).^2));
        end           
        function perimeter = get.perimeter(o)
                    perimeter = o.length + ...
                    sqrt(diff(o.X([1 end]))^2+diff(o.Y([1 end]))^2);
        end
        function area = get.area(o) 
            %AREA: computes surface area of a polyline considered as polygon
            % computation based on outer products of vectors spanned
            % between point in polygond and each line segment.

            xv=o.X(:); yv=o.Y(:);
            if xv(end)~=xv(1) || yv(end)~=yv(1)
                xv = [xv; xv(1)];
                yv = [yv; yv(2)];
            end

            xp=mean(xv);
            yp=mean(yv);

            dx = xv-xp;
            dy = yv-yp;

            dA = NaN(3,numel(xv)-1);

            for i = size(dA,2):-1:1
                dA(:,i) =  cross([dx(i); dy(i); 0], [dx(i+1); dy(i+1); 0]);
            end
            
            area = abs(sum(dA(3,:))/2);
        end
        function print2(o)
            %%PRINT2 -- prints out object so that it can be copied to Excel
            % this version is more general than print
            %
            % USAGE: kmlPathsObj.print2()
            %
            % TO 140226 Deyang
            
            for io=1:numel(o)
                N_ = numel(o(io).X);
                hdrs = o(io).fieldnames();
                for iN=numel(hdrs):-1:1
                    if numel(o(io).(hdrs{iN}))<N_
                        hdrs(iN)=[];
                    end
                end
                
                if io==1
                    fprintf('%s','name');
                    for iN=1:numel(hdrs)
                        fprintf('\t%s',hdrs{iN});
                    end
                fprintf('\n');
                end
                
                for iLine=1:N_
                    fprintf('%s',o(io).name);
                    for iN=1:numel(hdrs)
                        fprintf('\t%.8g',o(io).(hdrs{iN})(iLine));
                    end
                    fprintf('\n');
                end
                fprintf('\n');
            end
        end
    end
end

