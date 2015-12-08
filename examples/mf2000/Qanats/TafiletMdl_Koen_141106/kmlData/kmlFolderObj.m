%KMLFOLDEROBJ get the coordinates from a kml folder file with paths
%
% USAGE:
%     p = kmlFolderObj(kmlFolderFileName[,gr])
%     p = kmlFolderObj(kmlFolderFileName,gr,depthToe,depthHeel);
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
classdef kmlFolderObj
    properties
        name
        E
        N
        X
        Y
        zDem
        z
    end
    properties (Dependent=true)
        length
        perimeter
        area
        isClosed
    end
    methods
        function o = kmlFolderObj(varargin)
            %KMLFOLDEROBJ constructer of a kmlFolderObj
            %
            % USAGE obj  = kmlFolderObj(kmlFolderName)
            %
            %  kmlFolderName is name of kmlfile that contains a set of
            %  paths as obtained by copying a folder of paths from Google
            %  Earth to another directory.
            %
            % SEE ALSO kmlPathsObj kmlPathsObj
            %
            % TO 131010
            
            if nargin==0
                return;
            end
            
            [kmlFolderFileName,varargin] = getNext(varargin,'char',[]);
            if isempty(kmlFolderFileName)
                error('First argument must be the kmlFolderFleName');
            end
            [gr,varargin] = getType(varargin,'gridObj',[]);
                        
            LLstruct = kmlFolder(kmlFolderFileName);
            
%            o(numel(LLstruct,1))=kmlFolderObj();
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
        function o = getDem(o,gr)
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
            %SETELEVATION -- set the elevation of line objects such as
            % drainage tunnels (qanats, khettaras) given the beginning (toe)
            % and end (heel) depth below ground surface.
            %
            % USAGE  obj.setElevation(toeHeels);
            %
            % toeHeels is a cell aray in which each line has the name of
            % the object, its toeDepth and its heelDepth.
            % for instance:
            %   toeHeels = {'Sayed' 16 1.5
            %               'Fougania' 15.5 1.5
            %               'Lakdima' 17.2 1.5
            %               'Bouia Lakadima' 17.1 1.5
            %               'default'  17,1.5
            %              };
            %
            % Any object found with a name matching the first column of
            % toeHeels will use the given toe and heel depts. Any other
            % objects uses defaults. The final elevation is matched with a
            % desired inclination of 1/1000.

            
            [depthToe ,varargin] = getNext(varargin,'double',[]);
            [depthHeel,varargin] = getNext(varargin,'double',[]);
            [depthNames,      ~] = getNext(varargin,'cell',[]);
            
            if isempty(depthNames)
            
                for io=numel(o):-1:1
                    s = cumsum([0; sqrt(diff(o(io).X).^2+diff(o(io).Y).^2)]);
                    zToe  = o(io).zDem(  1)-depthToe;
                    zHeel = o(io).zDem(end)-depthHeel;
                    zToe  = min(o(io).zDem(1),...
                                max(zToe,zHeel + s(end)/5000)...
                                );
                    if numel(s)>1
                        o(io).z = interp1(s([1 end]),[zToe zHeel],s);
                    else
                        o(io).z = zToe;
                    end
                end
            else
                msgId = 'kmlFolder:cantfindname';
                for i=1:size(depthNames,1)
                    io = strmatchi(depthNames{i},{o.name});
                    depthToe = depthNames{i,2};
                    depthHeel= depthNames{i,3};
                    if ~io
                        warning('on',msgId);
                        warning(msgId,'%s: Can''t fine name <<%s>> in profile names',mfilename,depthNames{i});
                        warning('off',msgId);
                        continue;
                    end
                    s = cumsum([0; sqrt(diff(o(io).X).^2+diff(o(io).Y).^2)]);
                    
                    zToe  = o(io).zDem(1  )-depthToe;
                    zHeel = o(io).zDem(end)-depthHeel;
                    inclination = max(2500,s(end)/(zToe-zHeel));
                    zToe  = min(o(io).zDem(1),...
                                max(o(io).zDem(  1)-depthToe,...
                                    zHeel + s(end)/inclination)...
                                );
                    o(io).z = interp1(s([1 end]),[zToe zHeel],s);
                end
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
                    fprintf('%s',o(io).name);
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
            
            area = sum(dA(3,:))/2;
        end
    end
end

