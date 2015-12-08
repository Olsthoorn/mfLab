function o=stepResponse(o,t0,tau,varargin)
%STEPRESPONSE -- computes step response
%
% It is assumed that the input in head contains a step response starting at t=t0.
% that can be computed by subtracting the head for t>t0 from that in t0 or
% fromm the head in another piezometer, which can optionally be specified
% using parameter value pair 'base',observationObj.
% The step reponse is assumed to start in t0 and we want the result on the
% time series tau which starts at t0 with the value zero. Hence tau should
% start with 0 and represents t-t0.
%
% USAGE:
%   SR = observationObj.stepResponse(t0,tau)
%   SR = observationObj.stepResponse(t0,tau,'base',otherObservationObj)
%
% EXAMPLE:
%   run a model in which something starts at a given moment
%   while everything else is kept constant (e.g. rechage is increased by 
%   1 mm/d or a well is swithed on with a fixed Q or a river is suddenly
%   set to a new stage, which is constant thereafter.
%   piezom = observationObj(basename,sheetName,gr,'head',H)
%      where
%           basename the workbookName
%           sheetNm  the worksheet with the piezometer specifications
%           gr is a gridObj
%           H as obtained from readDat
%
%   make sure one observationObj contains the data of the base model
%  BR = piezomSR,SR2BR(t0,tau);
% or
%  BR = piezomSR.SR2BR(t0,tau,'base',piezomBase);
%
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
        hb = interp1(tb,hb,t,'cubic');
    end

    sr= interp1(t,h-hb,tau,'cubic');  if isnan(sr(1)), sr(1)=0; end
        
    o(io).data.stepResponse = [tau(abs(sr)>delta); sr(abs(sr)>delta)];
end
