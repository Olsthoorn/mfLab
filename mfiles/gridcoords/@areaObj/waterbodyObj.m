classdef waterbodyObj <hydObj
    properties
        havg
        c_entry
        z_bot
        izone  % zone numer in mesh for this waterbody
        izbot  % cell z-index of bottom of this waterbody 
        I      % indices of waterbody this waterbody in mesh
        Cond   % Conductance for cells I = dx*dy/ct
        LRC    % Lay Row Col of this waterbody in mesh (redundant)
		UserData % for any purpose
    end
    methods
        function obj=waterbodyObj(S,type)
            switch nargin
                case 0
                    return
                case 2
                    obj(numel(S),1)=waterbodyObj();
                    for i=1:numel(S)
                        obj(i).type     = type;
                        obj(i).name     = S(i).NAAM;
                        obj(i).shortnm  = S(i).NAAM;
                        obj(i).code     = S(i).CODE;
                        obj(i).x        = S(i).x;
                        obj(i).y        = S(i).y;
                        switch type
                            case {'canal'}
                                obj(i).nr    = S(i).KANAAL_AWD;
                                obj(i).tstart= datenum(S(i).START,1,1);
                                obj(i).tend  = datenum(S(i).EINDE,1,1);
                                obj(i).gauge = S(i).PEILSCHAAL;
                            case {'geul'}
                                obj(i).nr    = S(i).GEULEN_AWD;
                                obj(i).cbot  = S(i).C_BODEM;
                                obj(i).havg  = S(i).GEM_PEIL_7;
                                obj(i).gauge = S(i).GEULEN_AWD;
                        end
                    end
                otherwise
                    error('%s constructor wrong number of elements(%d)',class(obj),nargin);
            end
        end
        function obj=merge(obj,grid,izone)
            ZONE=zeros(grid.Ny,grid.Nx,grid.Nz);
            obj.izone=izone;
            obj.izbot = xyzindex(obj.z_bot,grid.zGr);
            obj.I=find(inpolygon(grid.XM(:,:,1),grid.YM(:,:,1),obj.x,obj.y));
            %figure; hold on
            %plot(obj.x,obj.y,'b');                     % debug
            %plot(grid.XM(:,:,1),grid.YM(:,:,1),'ro');  % debug
            if ~isempty(obj.I)
                for izb=1:obj.izbot
                    ZONE(obj.I+grid.Ny*grid.Nx*(izb-1))=obj.izone;
                end
                obj.I=find(ZONE);
                obj.LRC=cellIndices(obj.I,size(ZONE),'LRC');
                fprintf('waterbody %s nr(%d) in mesh\n',obj.shortnm,obj.izone);
            end
        end
        function obj=setCond(obj,grid)
            if ~isempty(obj.I)
                obj.Cond=(grid.DY(obj.I).*grid.DX(obj.I))/obj.c_entry;
            end
        end
        function h=head(obj,Piez,time)
            k=strmatchi(obj.gauge,{Piez.name},'exact');
            if k(1)>0
                h=interp1([-Inf; Piez(k).t(:); Inf] , Piez(k).h([1 1:end end] ),time);
            else
                error('Can''t find Piezometer % for waterbody %s,Piez(k)',name,obj.gauge);
            end
        end
    end
end