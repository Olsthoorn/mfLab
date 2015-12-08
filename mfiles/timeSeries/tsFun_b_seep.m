function [wv,v,n,y,yM,sigma,weights] = tsFun_b_seep(par,varargin)
    % tsFun_b_seep -- tsa simulator non lin, no drainage below drainage base
    % [wv,v,n,y,sigma] = tsFunMdl(par) -- non linear
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

    global t tm h

    if nargin==0
        % return default inital parameters in as
        % [wv =pDefaults, v=pLog n=pNames y=pDim]
        [wv,v,n,y] = getDefaults();
        return;
    end

    [t,y] = myModel(par);

    if nargin>1 % return only y as first argument
        wv=y;
        v=[]; n=[]; sigma=[]; weights=[];
        return;
    end

    %% Noise, part residuals and innovations

    yM       = interp1(t,y,tm); % intepolate simulated on measurement times
    
    dtm      = diff(tm); % measurement time steps
    v        = zeros(size(tm)); % innovations (v(tm)=n(tm)-n(tm-dtm)e(- dtm/T)
    T        = exp(par(end)); % decay factor or noise model
    n        = h-yM; % residual at measurement time
    v(1)     = n(1);
    v(2:end) = n(2:end) - n(1:end-1).*exp(-dtm/T); % dtm ??
    correct  = 1-exp(-2*[dtm(1); dtm]/T);

    %% Estimate sigma (see theory).
    % sigma^2 is the variance increase rate of the noise process
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

end

function [pDef,pLog,pNames,pDim] = getDefaults()
    % Default initial parameters (suitable for Terwisscha)
    d      = -1.0;
    q      = -1e-3;
    S      = 0.23;
    c      =  450;
    T      = S*c/3;
    pDef   = [d     q   log([S  c    T])]'; 
    pLog   = [false false true true true]';
    pNames = {'d' 'q'   'S' 'c' 'T'}';
    pDim   = {'d' 'm/d' '-' 'd' 'd'};
end

function [t_ y] = myModel(par,samplingDt,IN)
    % Simulate specific deterministic model

    global t P E
    
    if nargin>1
        [t_,P_,E_] = IN.prepareData(samplingDt);
    else
        t_=t;
        P_=P;
        E_=E;
    end
    
    d     =     par(1);  % top    drainage level
    q     =     par(2);  % seepage upward positive
    S     = exp(par(3)); % specific yield (storage coefficient)
    c     = exp(par(4)); % top    drainage resistance

    y  = zeros(size(t_)); % deterministic model estimate
    dt = diff(t_);        % simulation  time steps

    for it=1:numel(t_)    % simulate over all simulation times t
        if it==1          % initialize model
            y(it) = d;    % start model at drainage base
        else              % compute y
            if y(it-1)>d
                e = exp(-dt(it-1)/S/c);
                y(it) = d       + (P_(it)-E_(it)+q)*c *(1-e) + (y(it-1)-d)*e;
            else
                y(it) = y(it-1) + (P_(it)-E_(it)+q)/S*dt(it-1);
            end
        end
    end    
end