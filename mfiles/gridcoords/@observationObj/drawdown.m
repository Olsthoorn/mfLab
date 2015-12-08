function o=drawdown(o,t0,tau,varargin)
%DRAWDOWN -- puts diff of data into observation object (suitable for drawdowns
% etc. Interpolates on time scale
%
% The drawdown is relative to the head at t0
% The time of the object is set to o.data.head.time-t0
% while all times and data less than t0 are removed.
% tau is the times vector on which the data can be interpolated
%
% USAGE:
%   DD = observationObj.drawdown(t0,tau)
%   DD = observationObj.drawdown(t0,tau,'base',otherObservationObj)
%
% EXAMPLE:
%   piezom = observationObj(basename,sheetName,gr,'head',H)
%      where
%           basename the workbookName
%           sheetNm  the worksheet with the piezometer specifications
%           gr is a gridObj
%           H as obtained from readDat
%
%  BR = piezomSR,SR2BR(t0,tau);
% or
%  BR = piezomSR.SR2BR(t0,tau,'base',piezomBase);
%    in this second case, 'base',piezomBase uses another piezometer to
%    subtract from the heads of the objects instead of the head at t0.
%    Both piezometer objects need not have the same time base but should overlap.
%

delta = 1e-4;

[base,~] = getProp(varargin,'base',[]);
for io=numel(o):-1:1
    t = o(io).data.head(1,:);
    h = o(io).data.head(2,:);
    h0= interp1(t,h,t0);

    if isempty(base)
        tb = t;
        hb = h0*ones(size(h));
    else
        tb = base(min(numel(base),io)).data.head(1,:);
        hb = base(min(numel(base),io)).data.head(2,:);
    end

    h  = h( t >=t0);
    hb = hb(tb>=t0);
    
    t  = t( t >=t0);
    tb = tb(tb>=t0);
   
    t  = t - t0;
    tb = tb- t0;
    
    tau = [0 tau(tau>0 & tau<=t(end))];
    if isempty(tau)
        error('%s: tau is empty, it must overlap time',mfilename);
    end
    
    if any(diff(tb)>delta)
        hb = interp1(tb,hb,t);
    end

    d = interp1(t,h-hb,tau);

    o(io).data.drawdown = [tau; d];
end
