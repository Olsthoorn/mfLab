function SR = menyanthesSR(A,a,n,b,t,varargin)
%% menyanthesSR -- compute IR, SR and BR for Menyanthes simulations
%
% USAGE:
%    SR = menyanthesSR(A,a,n,b,t);       % step response
%    SR = menyanthesSR(A,a,n,b,t,'SR');  % step response
%    BR = menyanthesSR(A,a,n,b,t,'BR');  % block response
%         menyanthesSR();  % selftest
%    IR = menyanthesSR(A,a,n,b,t,'IR');  % impulse response
%
% TO141021
    if nargin==0
        selftest;
        SR = [];
        return;
    end


    [ir,varargin] = getWord(varargin,'IR');
    [sr,varargin] = getWord(varargin,'SR');
    [br,  ~     ] = getWord(varargin,'BR');

    if ~any([ir sr br]), sr=true; end
    
    if ir
        SR = A * a^(n) * (t-b).^(n-1).*exp(-a*(t-b))/gamma(n);
        SR(t<=b)=0;
        return;
    elseif br
        SR = menyanthesSR(A,a,n,b,t);
        SR = SR-[0 SR(1:end-1)];
    elseif sr
        TOL =1e-4; maxIter = 20; SR1=Inf;
        tauM = maxTau(A,a,n,b,t,TOL);
        for i=1:maxIter;
            dtau = tauM/(i^2*1000);
            tau = 0:dtau:tauM; tau(end)=tauM;
            IR  = menyanthesSR(A,a,n,b,tau,'IR');
            SR  = [0 cumsum(0.5*(IR(1:end-1)+IR(2:end)).*diff(tau))];
          %  fprintf('%2d %12.2f %12.6f %12.6f %12.6e\n',i,tauM,SR(end),SR1(end),SR(end)-SR1(end));
            if abs(SR(end)-SR1(end))<TOL
                break;
            else
                SR1=SR;
            end
        end
        if i==maxIter, error('no convergence reached in computing SR'); end
        SR =interp1(tau,SR,t);
    else
        error('see help of this function');
    end
end

function tau = maxTau(A,a,n,b,t,TOL)
% find maximum tau for numeric integration of Menyanthes' IR
%
% USAGE tau= maxTau(A,a,n,b,t)
%
% find t for which  t * theta(t) < TOL
% where theta(t) is menyanthes IR:
%
% theta(t) = A*a^n*(t-b)^(n-1)*exp(-a*(t-b))/gamma(n)
%
% with y=t-b and taking the log on both sides
% F(y)    = n * log(y) - a*y - ln(TOL/A/a^n*gamma(n))
% Facc(y) = n*(1-log(y))-a;
% y(k) = y(k-1) - F(y)/Facc(y); %newton raphson
%
% TO 141022

    if nargin<5, TOL=1e-6; end
    tau = max(t(:));  % initial estimate

    y = tau-b;

    for k=1:100
        F   = n*log(y)-a*y-log(gamma(n)/A/a^n*TOL);
        Facc= n/y-a;
        dy = F/Facc;
        if abs(dy)<TOL
            tau = y+b;
            if tau>t(end)  % no need to go beyond t(end)
                tau=t(end);
            end
            return
        else
            y=y-dy;
        end
    end
    error('root not found after %d iterations',k);
end
        
function selftest
    n = [1 1.3 1.7 2.3]; % n>=1, n<1 gives trouble and is unnecessary
    A0= 100;   % gain A>0 is this example A0*n(i)
    a = 0.01;  % a>0
    b = 100;   % delay
    Dt= 7;     % width of block response
    t = 0:Dt:300;
    
    figure('pos',screenPos(1,0.5)); hold on;
    
    %% Impulse response (See Von Asmuth (2012, p76)
    ax1 = subplot(1,3,1,'nextPlot','add','xGrid','on','yGrid','on','ylim',[0 1]);
    leg=cell(1,numel(n));
     xlabel(ax1,'time [d]'); ylabel(ax1,'response (-)');
    title(ax1,'Menyanthes Impuls Responses');
    for i=1:numel(n)
        plot(ax1,t,menyanthesSR(A0*n(i),a,n(i),b,t,'IR'),mf_color(i));
        leg{i} = sprintf('n = %.1f',n(i));
    end
    legend(ax1,leg);

    %% Step response
    ax2 = subplot(1,3,2,'nextPlot','add','xGrid','on','yGrid','on');
    leg=cell(1,numel(n));
    xlabel(ax2,'time [d]'); ylabel(ax2,'response (-)');
    title('Menyanthes Step Responses');
    for i=1:numel(n)
        plot(ax2,t,menyanthesSR(A0*n(i),a,n(i),b,t,'SR'),mf_color(i));
        leg{i} = sprintf('n = %.1f',n(i));
    end
    legend(ax2,leg);
    
    %% Block response
    ax3 = subplot(1,3,3,'nextPlot','add','xGrid','on','yGrid','on');
    leg=cell(1,numel(n));
    xlabel(ax3,'time [d]'); ylabel(ax3,'response (-)');
    title(ax3,sprintf('Menyanthes Block Responses, Dt=%.2d',t(2)));
    for i=1:numel(n)
        plot(ax3,t,menyanthesSR(A0*n(i),a,n(i),b,t,'BR'),mf_color(i));
        leg{i} = sprintf('n = %.1f',n(i));
    end
    legend(ax3,leg);
end
