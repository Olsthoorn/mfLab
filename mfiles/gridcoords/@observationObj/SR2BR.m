function o=SR2BR(o,t0,tau,varargin)
%SR2BR -- converts Step Response to Block Response implied by tau
%
% It is assumed that the input is a step response starting at t=t0.
% The time of the o is set to o.data.head.time-t0
% while all times and data less than t0 are removed.
% tau is the times vector for the block response. It should have
% fixed step size and start at 0. If not 0 is prefixed.
% The step response is relative to the data at t=t0.
% It is possilble to compute the step response realtive to another
% observation object, that represents the results from a base run. To do
% that use 'base',observationobj as a parameter name value pair in the
% input.
%
% USAGE:
%   BR = observationObj.SR2BR(t0,tau)
%   BR = observationObj.SR2BR(t0,tau,'base',otherObservationObj)
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
        hb = interp1(tb,hb,t);
    end

    sr= interp1(t,h-hb,tau);  if isnan(sr(1)), sr(1)=0; end
    
    br= [0 diff(sr)];
    
    o(io).data.blockResponse = [tau(abs(br)>delta); br(abs(br)>delta)];
end
