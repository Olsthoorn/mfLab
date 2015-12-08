function [wv,v,n,y,yM,sigma,weights] = tsFun_a(par,tm,h,t,P,E)
    % tsFun_a_seep -- sime seris model standard fixed drainage level
    % [wv,v,n,y,sigma] = tsFunMdl(par,tm,h,IN)
    % Continuous time series model with irregular time steps.
    % input : par parameter vector
    % output: wv(tm)  weighted innovations,, tm indicates measurement times
    %         v(tm)   innovations
    %         n(tm)   residuals (h-y)
    %         y(t)    estimate of deterministic model
    %         yM(tm)  y(t) interpolated on tm
    %         sigma^2 variance due to noise process after unit time (one day)
    %         
    % TO 141010

    %global tm h t P E
    
    if nargin==0
        % return default inital parameters in as
        % [wv =pDefaults, v=pLog n=pNames y=pDim]
        [wv,v,n,y] = getDefaults();
        return;
    end
    
    %% deterministic model
    y = simulate(par);
    
    if isempty(tm)
        wv = []; v=[]; n=[]; yM=[]; sigma=[]; weights=[];
        return
    end

    %% noise model, residuals and innovations
    yM       = interp1(t,y,tm); % intepolate simulated on measurement times
    
    dtm      = diff(tm); % measurement time steps
    v        = zeros(size(tm)); % innovations (v(tm)=n(tm)-n(tm-dtm)e(-dtm/T)
    T        = exp(par(end)); % decay factor or noise model
    n        = h-yM; % residual at measurement time
    v(1)     = n(1);
    v(2:end) = n(2:end) - n(1:end-1).*exp(-dtm/T); % dtm ??
    correct  = 1-exp(-2*[dtm(1); dtm]/T);

    %% Estimate sigma; sigma^2 is the variance increase rate of the noise process (see theory)
    sigma    = sqrt( 2 * mean( v.^2 ./ correct) / T );

    %% Solve either by LSE or MLE
    LSE = false;
    if LSE % Least squares estimation
        wv      = v; %#ok
        weights = ones(size(v));
    else % Maximum Likelihood estimation
        N       = numel(v);
        F       = prod(correct.^(1/N));
        weights = sqrt(F./correct);
        wv      = weights .* v;
    end

    function y = simulate(par)
        % Simulate specific deterministic model

        d     =     par(1);  % top    drainage level
        S     = exp(par(2)); % specific yield (storage coefficient)
        c     = exp(par(3)); % top    drainage resistance

        y  = zeros(size(t)); % deterministic model estimate
        dt = diff(t);         % simulation  time steps

        for it=1:numel(t)    % simulate over all simulation times t
            if it==1         % initialize model
                y(it) = d;  % start model at drainage base
            else             % compute y
                e = exp(-dt(it-1)/S/c);
                y(it) = d + (P(it)-E(it))*c *(1-e) + (y(it-1)-d)*e;
            end
        end    
    end

end
function [pDef,pLog,pNames,pDim] = getDefaults()
    % Default initial parameters (Suitable for Terwisscha)
    d      = -1.0;
    S      = 0.23;
    c      = 450;
    T      = S*c/3;
    pDef   = [d   log([S c    T])]'; 
    pLog   = [false true true true]';
    pNames = {'d' 'S' 'c' 'T'}';
    pDim   = {'d' '-' 'd' 'd'};
end

