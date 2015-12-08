function y = modelWrapper(p)
%%MODELWRAPPER -- computes measurements - modelOutcome
% all p are assumed to be the ln of the original parameters
% so exp(p) will be passed to the context model function below.
%
% USAGE:
%      e = modelWrapper(p)
% p is a vector with parameter values
%     the order of p must be in the order of unused parameters in the
%     the vector used as described below.
% e is measurement - modelOutcome
%
% other parameters are passed through global
%   global par  -- full set of parametervalues
%   used        -- set in par that is used
%       this is a vector or logicals or integers where 0 means false
%       the length of use is the full set or parameter. The total number
%       of falses must equal the length of p. The order of the false must
%       be the same as the order of p.
%       Whereas the parameters passed by p are assumed the ln of the
%       original parameters, the parameters passed by are are not. They are
%       use as passed.
%   measurements -- measurements a struct that allows the model
%      to know what info to extract from the model output
%
global Par defaults meas iter % all initial parameters

%% Assemble the full set of parameters

[~,~,~,pars] = Par.Par2p(p);

y  = model('x',meas.x,pars{:},defaults{:}) - meas.y;

iter = iter+1;
switch rem(iter,length(p));
    case 0,
        fprintf(' %10.4g',sqrt(y'*y));
        fprintf('\n');
    case 1, 
        fprintf('Iter%5d',iter); 
        fprintf(' %10.4g',sqrt(y'*y));
    otherwise
        fprintf(' %10.4g',sqrt(y'*y));
end

