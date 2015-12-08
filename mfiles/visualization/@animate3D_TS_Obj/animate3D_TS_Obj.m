classdef animate3D_TS_Obj
%% animate 3D concentration and temperature
% TO 120420 120901
    properties
        basename
        titleStr             ='My title string';
        numberOfHeadContours = 50;
        numberOfTempContours = 50;
        numberOfConcContours = 50;
        numberOfStreamLines  = 50;
        concThreshold,
        xLim
        yLim
        framerate            = 3;
        quality              = 75;
        
        H, % heads
        T, % temperatures
        C, % concentrations
        B, % budget
        well, % wells
        ax,, % axes
        TSTART, % start temperature
        CSTART, % start concentration
    end
    properties (Dependent=true)
        hrange, % contour range of heads
        trange, % contour range for temperatures
        crange, % contour range of concentrations
        prange, % contour range of stream function
        dtrange, % contour range of temp relative to initial values
        dcrange, % contour range of conc relative to initial values
        cMax    % max concentration
        cMin    % min concentration
        cLim
        Tmax    % max temperature
        Tmin    % min temperature
        Tlim
        hMax
        hMin
        hLim
        DT,
        DC,

    end
    methods
        function o=animate3D_TS_Obj(basename,well,STCONC)
            o.basename = basename;
            o.H=readDat([basename '.HDS']); o.H=maskHC(o.H,[-1000,1000],[NaN,NaN]);
            o.T=readMT3D('MT3D001.UCN');    o.T=maskHC(o.T,[o.Tmin o.Tmax]);
            o.C=readMT3D('MT3D002.UCN');    o.C=maskHC(o.C,[o.cMin o.cMax]);
            o.B=readBud([basename,'.bgt']);  % get only flow rightface
            
            o.B= mf_Psi(o.B);
            
            if exist('well','var');
                o.well=well.setCout(o.C,1);
            end
            
            if nargin>2
                o.TSTART = STCONC{1};
                o.CSTART = STCONC{2};
            end
        end
        function hrange = get.hrange(o)
            hrange = ContourRange(o.H,o.numberOfHeadContours);
        end
        function crange = get.crange(o)
            crange = [0 ContourRange(o.C,o.numberOfConcContours)];
        end
        function trange = get.trange(o)
            trange = ContourRange(o.T,o.numberOfTempContours);
        end
        function prange = get.prange(o)
            prange = ContourRange(o.B,o.numberOfStreamLines,[],'Psi');
        end
        function dtrange = get.dtrange(o)
            dtrange = ContourRange(o.B,o.numberOfTempContours,[],'Psi');
        end
        function dcrange = get.dcrange(o)
            dcrange = ContourRange(o.B,o.numberOfConcContours,[],'Psi');
        end
        
        function cMax = get.cMax(o), cMax = o.crange(end);      end
        function cMin = get.cMin(o), cMin = max(0,o.crange(1)); end
        function cLim = get.cLim(o), cLim = [o.cMin o.cMax]; end
        function Tmax = get.Tmax(o), Tmax = max(o.trange);   end
        function Tmin = get.Tmin(o), Tmin = min(o.trange);   end
        function Tlim = get.Tlim(o), Tlim = [o.Tmin o.Tmax]; end
        function hMax = get.hMax(o), hMax = max(o.hrange);   end
        function hMin = get.hMin(o), hMin = min(o.hrange);   end
        function hLim = get.hLim(o), hLim = [o.hMin o.hMax]; end
        function DT   = get.DT(o)
            DT = o.T;
            for it=1:length(DT),
                DT(it).values = o.T(it).values - o.TSTART;
            end
        end
        function DC   = get.DC(o)
            DC = o.C;
            for it=1:length(DC),
                DC(it).values = o.C(it).values - o.CSTART;
            end
        end
        
    end
end

        