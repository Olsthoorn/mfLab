function o=blockResponse(o,t0,tau,varargin)
%BLOCKRESPONSE -- computes block response
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
%   BR = observationObj.blockResponse(t0,tau)
%   BR = observationObj.blockResponse(t0,tau,'base',otherObservationObj)
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

o = o.stepResponse(t0,tau,varargin{:});

for io=numel(o):-1:1

    tau =         o(io).data.stepResponse(1,:);
    br  = [0 diff(o(io).data.stepResponse(2,:))];

    o(io).data.blockResponse = [tau(abs(br)>delta); br(abs(br)>delta)];
end
