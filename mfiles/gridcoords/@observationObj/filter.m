function o=filter(o,tp,F)
%FILTER -- filters o.data.blockResponse
% simulates the input time series by with input series tp
% tp is a time series of the "driving force" having the same time step as
% the block response that is contained in o(io).data.blockReposponse(1,:);
%
% USAGE:
%   observationObj = observationObj.filter(tp[,F])
%
%   tp is a time series with first column time and second column the value
%   of the driving force. F is a multiplyer for the driving force, which
%   may be 1. This can be used to scale the block response.
%
% EXAMPLE:
%   piezom = observationObj(basename,sheetName,gr,'head',H)
%      where
%           basename the workbookName
%           sheetNm  the worksheet with the piezometer specifications
%           gr is a gridObj
%           H as obtained from readDat
%   Then get a block response from a time series using
%   observationObj.blockResponse. For that you need a time series of heads which
%   contains a step response starting at some t0 and optionally a time
%   series of heads that are not affected.
%   BR = piezom.blockResponse(t0,tau); or
%   BR = piezom.blockResponse(t0,tau,'base',piezom(j));
%   Then use filter
%   h = BR.filter(tp);
%   where tp is a time series of rain or well pumping
%

delta = 1e-3;

if nargin<3, F=1; end

for io=numel(o):-1:1
    tau = o(io).data.blockResponse(1,:);
    br  = o(io).data.blockResponse(2,:);

    t   = tp(:,1)';
    p   = tp(:,2)';

    if abs(mean(diff(tau))-mean(diff(t)))>delta
        error('%s: time stepsize must equal that of block response',mfilename);
    end
    
    h   = filter(br*F,1,p);
    
    o(io).data.filtered = [t; h];
end
