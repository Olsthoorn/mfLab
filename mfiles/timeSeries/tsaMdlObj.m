classdef tsaMdlObj % --- time series analysis model object
    properties
        name      % [-] name of this model
        date      % [d] date of creation
        xcoord    % [m] coordinate of well
        ycoord    % [m] coordinate of well
        fun       % [-] pointer to simulation function
        Hname     % [-] always only one
        INname    % [-] plural mostly
        INtype    % [-] same number as INname
        pIn       % [?] initial par vector
        p         % [?] final par vector
        pNames    % [-] parmeter names
        pLog      % [-] whether or not to take log
        pDim      % [?] dimension of parameters
        tm        % [d] measurement times
        in        % [-] boolean, true where tm>=t(1) & tm<=t(end)
        t         % [d] simulation times
        y         % [m] y(t)  simulated deterministic model, result of calibration
        yM        % [m] yM(tm, y(t) interpolated on tm, measurement times
        h         % [m] h(tm) head measurements
        n         % [m] h(tm) residuals
        v         % [m] h(tm) innovations
        sigma     % [m/root(d)] root of varaince increase rate of noise process
        weights   % [-] MLE weights
        J         % [?]from lsqnonlin
        resnorm   % [m^2] from lsqnonlin
        residual  % [m] from lsqnonlin
        exitflag  % [-] from lsqnonlin
        R         % [m] distance (e.g. from pumping station)
    end
    properties (Dependent = true)
        dt        % [d] simulation time step series
        dtm       % [d] measurement time step series
        Cov       % [?] parameter covariance matrix
        Std       % [?] parameter standard deviations
        Cor       % [-] parameter correlation matrix
        Rn        % [-] autocorr(o(io).n);
        Rv        % [-] autocorr(o(io).v);
        Ry        % [-] autocorr(o(io).y);
        Rh        % [-] autocorr(o(io).h);
        Rvy       % [-] crosscorr(o(io).v,o(io).yM
        Rvn       % [-] crosscorr(o(io).v,nOtilde
        explVar   % [-] explained variance
        d         % [m] drainage depth
        q         % [m/d] upward seepage
        S         % [-] specific yield
        c         % [d] drainage resistance
        T         % [d] characteristic time of noise process
        yMean     % [m] mean of simulated values

    end
    methods
        function o=tsaMdlObj(fun,pIn,H,IN)
            % o=tsaMdlObj -- model for time series analysis
            %
            % USAGE:    o = tsaMdlObj(fun,pIn,H,IN)
            %
            % TO 141028
            
            if nargin==0; return; end
            
            for io=numel(H):-1:1
                o(io).name    = H(io).name;          % default name
                o(io).tm      = H(io).values(:,1);   % measurment times
                o(io).h       = H(io).values(:,end); % measurements
                o(io).xcoord  = H(io).xcoord;
                o(io).ycoord  = H(io).ycoord;
                o(io).R       = H(io).R;    % distance to PS or whatever, obtained from H 
                o(io).date    = now();      % stamp
                o(io).fun     = fun;        % funcion used to calibrate
                o(io).pIn     = pIn;        % initial parameters
                
                % parameter log, names, dims
                [~,o(io).pLog,o(io).pNames,o(io).pDim] = fun();
                
                o(io).Hname   = H(io).name; % only name
                o(io).INname  = {IN.name};  % only name
                o(io).INtype  = {IN.type};  % only name
            end
        end
        function showPar(o)
            % show the parameter values
            mode = {'(lin)' '(log)'};
            fprintf('\nParameters:\n');
            warning('off','all');
            for io=1:numel(o)
                fprintf('%s\n',o(io).name);
                for i=1:numel(o(io).p)
                    if o(io).pLog(i)
                        fprintf('%5s   %s [%3s] %10g std: %10.0f%%\n',...
                            o(io).pNames{i},...       % parameter name
                            mode{o(io).pLog(i)+1},... %'(lin)' or '(log)'
                            o(io).pDim{i},...         % [ dim ]
                            exp(o(io).p(i)),...            % value
                            round(100*(exp(o(io).Std(i))-1)));            % std dev
                    else
                        fprintf('%5s   %s [%3s] %10g std: %g\n',...
                            o(io).pNames{i},...       % parameter name
                            mode{o(io).pLog(i)+1},... %'(lin)' or '(log)'
                            o(io).pDim{i},...         % [ dim ]
                            o(io).p(i),...            % value
                            o(io).Std(i));            % std dev
                    end
                end
                fprintf('%5s %s [%3s] %7.1f\n','sigma^2','(lin)','cm^2/d',(o(io).sigma*100)^2);
                fprintf('%5s %s [%3s] %10.0f%%\n','explVar','(lin)',' - ',o(io).explVar*100);
                fprintf('\n');                
            end
            warning('on','all');
        end
         function par(o)
            % show the parameter values
            Par = NaN(numel(o),numel(o(1).pNames) + 2);
            N = numel(o(1).p);
            for i=N:-1:1
                hdr{2*i-1} = [o(1).pNames{i},'[',o(1).pDim{i},']'];
                hdr{2*i  } = ['std[',o(1).pDim{i},'])'];
            end
            hdr{2*N+1} = '[sigma^2 [m^2/d]]';
            hdr{2*N+2} = '[explVar [fraction]]';            
            
            fprintf(' %s\t',hdr{1:end-1});
            fprintf(' %s\n',hdr{  end  });

            warning('off','all');
            for io=1:numel(o)
                for i=1:numel(o(io).p)
                    if o(io).pLog(i)
                        Par(io,2*i-1) = exp(o(io).p(i));
                        Par(io,2*i  ) = exp(o(io).Std(i))-1;
                    else
                        Par(io,2*i-1) = o(io).p(i);
                        Par(io,2*i  ) = o(io).Std(i);            % std dev
                    end
                end
                Par(io,2*i+1) = o(io).sigma^2;
                Par(io,2*i+2) = o(io).explVar;
            end
            warning('on','all');
            fmt = repmat('%.3g\t',[1,size(Par,2)]); fmt([end-1 end])='\n';
            fprintf(fmt,Par);
         end
        function o = calibrate(o,t,P,E)
            % tsaMdlObj/calibrate -- calibrate tsa model
            %
            % USAGE: mdl = mdl.calibrate(t,P,E)
            %  mdl = tsaModlObj array
            %  t  = time
            %  P  = preciitation
            %  E  = evapo-transpiration
            %
            % TO 141029           
                                 
            Nfiles= numel(o);
            
            fprintf('Optimizing %d time series:\n',Nfiles);
            
            for io=Nfiles:-1:1
                                
                o(io).in = o(io).tm>=t(1) & o(io).tm<=t(end);
                o(io).t  = t;
                                
                fun_ = @(p)  o(io).fun(p,o(io).tm(o(io).in),o(io).h(o(io).in),t,P,E);                
                
                % Make use of options (may be omitted, just to show how)
                lb   = -Inf * ones(size(o(io).pIn));
                ub   =  Inf * ones(size(o(io).pIn));
                options = optimset('TolX',1e-6,'TolFun',1e-6,'Display','off');  % default TolX = 1e-6

                % Optimize parameters
                [o(io).p,o(io).resnorm,o(io).residual,o(io).exitflag,~,~,o(io).J] = ...
                    lsqnonlin(fun_,o(io).pIn,lb,ub,options);

                % Simulate including noise model
                [~,o(io).v,o(io).n,o(io).y,o(io).yM,o(io).sigma,o(io).weights] = ...
                    o(io).fun(o(io).p,o(io).tm(o(io).in),o(io).h(o(io).in),t,P,E); % Simulate with final   parameters
                
                fprintf('.'); if rem(Nfiles-io+1,50)==0, fprintf('%d\n',Nfiles-io+1); end
            end
            if rem(Nfiles-io+1,50)~=0, fprintf('%d\n',Nfiles-io+1); end
            fprintf(' done !\n');  
        end
        function showStat(o)
            % Show statistics of calibration
            for io=1:numel(o)
                fprintf('%s, parameter statistics\n',o(io).name);
                fprintf('Cov ='); display(o(io).Cov);
                fprintf('Std ='); display(o(io).Std);
                fprintf('Cor =');display(o(io).Cor);
            end
        end
        function plotXCor(o)
            fsz = 14; defaults = {'nextPlot','add','fontSize',fsz,'xGrid','on','yGrid','on'};            
            for io=1:numel(o)
                figure('name',o(io).name,'pos',screenPos(0.8,0.5));
                ax1 = subplot(1,2,1,defaults{:});
                xlabel(ax1,'v'); ylabel(ax1,'y'); title(ax1,'Rvy');
                plot(ax1,o(io).v,o(io).yM,'.');
                
                ax2 = subplot(1,2,2,defaults{:});
                xlabel(ax2,'v'); ylabel(ax2,'n'); title(ax2,'Rvn');

                nTilde = [o(io).n(1); o(io).n(1:end-1).*(1-exp(-diff( o(io).tm( o(io).in ) ) / o(io).T))];
                
                plot(ax2,o(io).v,nTilde,'.');
            end
        end
        function plotCor(o)
            fsz = 14;
            for io=1:numel(o)                
                dflt = {'nextPlot','add','fontsize',fsz,'xGrid','on','yGrid','on'};
                
                ierr = 0;
                if isempty(o(io).Rn),  fprintf('Rn not defined in model %s\n',o(io).name);  ierr=ierr+1; end
                if isempty(o(io).Rv),  fprintf('Rv not defined in model %s\n',o(io).name);  ierr=ierr+1; end
                if isempty(o(io).Ry),  fprintf('Ry not defined in model %s\n',o(io).name);  ierr=ierr+1; end
                if isempty(o(io).Rvy), fprintf('Rvy not defined in model %s\n',o(io).name); ierr=ierr+1; end
                if isempty(o(io).Rvn), fprintf('Rvn not defined in model %s\n',o(io).name); ierr=ierr+1; end
                if ierr>0, return; end
                
                figure('name',o(io).name,'pos',screenPos(0.8));
                subplot(5,1,1,dflt{:}); plot(1:numel(o(io).Rn) ,o(io).Rn,'bo-'); ylabel('Rn'); title('autocorrelation');
                subplot(5,1,2,dflt{:}); plot(1:numel(o(io).Rv) ,o(io).Rv,'bo-'); ylabel('Rv');
                subplot(5,1,3,dflt{:}); plot(1:numel(o(io).Ry) ,o(io).Ry,'bo-'); ylabel('Ry');
                subplot(5,1,4,dflt{:}); plot(1:numel(o(io).Rvy),o(io).Rvy,'bo-');ylabel('Rvy'); title('cross correlation');
                subplot(5,1,5,dflt{:}); plot(1:numel(o(io).Rvn),o(io).Rvn,'bo-');ylabel('Rvn');
            end
        end
        function ax = showResults(o,varargin)
            %% Show the optimized results
            %
            % USAGE    ax = mdl.showResults       % plots and returns axes
            %             = mdl.showResults(ax);  % uses ax to plot on
            % TO 141029
            
            [ax   , ~      ] = getType(varargin,'axis',[]); 
            newFig = isempty(ax);
            fsz      = 14;
            
            defaults = {'nextPlot','add','fontsize',fsz,'xGrid','on','yGrid','on'};
            for io=1:numel(o)
                if newFig
                    figure('name',o(io).name,'pos',screenPos(0.8));
                    ax(1) = subplot(3,1,1,defaults{:},'xlim',o(io).tm([1 end]));
                        plot(    ax(1),o(io).tm,o(io).h,'r.');
                        plot(    ax(1),o(io).t ,o(io).y,'b');
                        legend(  ax(1),'head','model'); %,'initialModel');
                        title(   ax(1),sprintf('measured (h) and modeled (y) Model: %s',o(io).name));
                        datetick(ax(1),'x','yyyy','keeplimits');
                    ax(2) = subplot(3,1,2,defaults{:},'xlim',o(io).tm([1 end]));
                        plot(    ax(2),o(io).tm,o(io).n,'g');
                        legend(  ax(2),'n');
                        title(   ax(2),'residuals');
                        datetick(ax(2),'x','yyyy','keeplimits');
                    ax(3) = subplot(3,1,3,defaults{:},'xlim',o(io).tm([1 end]));
                        plot(    ax(3),o(io).tm,o(io).v,'k');
                        legend(  ax(3),'v');
                        title(   ax(3),'innovations');
                        datetick(ax(3),'x','yyyy','keeplimits');
                else
                    plot(    ax(min(1,numel(ax))),o(io).tm,o(io).h,'r.');
                    plot(    ax(min(1,numel(ax))),o(io).tm,o(io).yM,'b');

                    plot(    ax(min(2,numel(ax))),o(io).tm,o(io).n,'g');
                    plot(    ax(min(3,numel(ax))),o(io).tm,o(io).v,'k');
                end
            end
        end        
        function o = sim(o,t,P,E)
            % tsMdl/sim -- simulates the time series using the model
            %
            % USAGE tsaMdlObj = tsaMdlObj.sim(IN)
            %
            % TO 141028
            
            % [wv,v,n,y,yM,sigma,weights] =  ...
            %             fun(pIn,tm,h,t,P,E);  % Simulate with initial parameters
            for io=1:numel(o)
                o(io).t     = t;
                [~,o(io).v,o(io).n,o(io).y,o(io).yM,o(io).sigma,~] = o(io).fun(o(io).p, o(io).tm, o(io).h, t,P,E);
            end
        end
        function showSimulated(o)
            % show simulation results
            figure('name','showSimulated','pos',screenPos(0.8,0.6));
            ax = axes('nextPlot','add','xGrid','on','yGrid','on');
            xlabel(ax,'t [d]'); ylabel(ax,'y [m]'); title(ax,'simulated heads');
            for io=1:numel(o)
                plot(ax,o(io).t,o(io).ySim,mf_color(io));
            end
            datetick(ax,'x','yyyy','keeplimits');
            legend(ax,{o.name});
        end
        function plot(o,varargin)
            [ax,varargin] = getType(varargin,'axis',gca);
            if isempty(varargin), end % dummy, skip
            
            ierr=0;
            if isempty(o.y), fprintf('y not defined in mdl %d\n',o.name); ierr=ierr+1; end
            if isempty(o.h), fprintf('h not defined in mdl %d\n',o.name); ierr=ierr+1; end
            if isempty(o.n), fprintf('n not defined in mdl %d\n',o.name); ierr=ierr+1; end
            if isempty(o.v), fprintf('v not defined in mdl %d\n',o.name); ierr=ierr+1; end
            if ierr>0, return; end
                                    
            for i=1:numel(o)
                plot(ax,o.t ,o.y,'b'); % varargin{:});
                plot(ax,o.tm,o.h,'r.');% varargin{:});
                plot(ax,o.tm,o.n,'g'); % varargin{:});
                plot(ax,o.tm,o.v,'m'); % varargin{:});
            end
        end
        function [R,o] = distance(o,x,y)
            % compute distance from well location to point x,y
            % add distance to fields in o
            % TO 141105
            R = sqrt(([o.xcoord]-x).^2 + (([o.ycoord]-y).^2));
            for io=numel(o):-1:1
                o(io).R = R(io);
            end
        end
        %% Compute and show autocorrelations and cross correlations
        function d = get.d(o)
            j = strmatchi('d',o.pNames,'exact');
            if o.pLog(j)
                d = exp(o.p(j));
            else
                d = o.p(j);
            end
        end
        function S = get.S(o)
            j = strmatchi('S',o.pNames,'exact');
            if o.pLog(j)
                S = exp(o.p(j));
            else
                S = o.p(j);
            end
        end
        function c = get.c(o)
            j = strmatchi('c',o.pNames,'exact');
            if o.pLog(j)
                c = exp(o.p(j));
            else
                c = o.p(j);
            end
        end
        function T = get.T(o)
            j = strmatchi('T',o.pNames,'exact');
            if o.pLog(j)
                T = exp(o.p(j));
            else
                T = o.p(j);
            end
        end
        function q = get.q(o)
            j = strmatchi('q',o.pNames,'exact');
            if ~j
                q = NaN;
            else
                if o.pLog(j)
                    q = exp(o.p(j));
                else
                    q = o.p(j);
                end
            end
        end
        function Rn = get.Rn(o), if isempty(o.n), Rn=[]; else Rn = autocorr(o.n); end, end
        function Rv = get.Rv(o), if isempty(o.v), Rv=[]; else Rv = autocorr(o.v); end, end
        function Ry = get.Ry(o), if isempty(o.y), Ry=[]; else Ry = autocorr(o.y); end, end
        function Rh = get.Rh(o), if isempty(o.h), Rh=[]; else Rh = autocorr(o.h); end, end
        function Rvy = get.Rvy(o)
            if isempty(o.v) || isempty(o.y)
                Rvy=[];
            else
                Rvy = crosscorr(o.v,o.yM);
            end
        end
        function Rvn = get.Rvn(o)
            if isempty(o.n)
                Rvn=[];
            else
                % unbiased between innovation and estimate of n, ntilde
                nTilde = [o.n(1); o.n(1:end-1).*(1-exp(-diff(o.tm)/o.T))];
                Rvn = crosscorr(o.v,nTilde);  % n is dependent on v
            end
        end
        function Cov = get.Cov(o)
            % covariance matrix
            if isempty(o.v) || isempty(o.J)
                Cov=[];
            else
                Cov = var(o.v).*((full(o.J)'*full(o.J))^(-1));
            end
        end
        function Std = get.Std(o)
            % parameters standard deviation
            if isempty(o.Cov)
                Std = [];
            else
                Std = sqrt(diag(o.Cov));
            end
        end
        function Cor = get.Cor(o)
            % correlation matrix
            if isempty(o.Std)
                Cor=[];
            else
                Cor = o.Cov./(o.Std*o.Std');
            end
        end
        function explVar = get.explVar(o)
            % explained variable
            if isempty(o.h) || isempty(o.y)
                explVar=[];
            else
                explVar = 1-var(o.n)/var(o.h);
            end
        end
        function yMean = get.yMean(o), yMean = mean(o.y); end
    end
end