classdef parObj
%%PAROBJ -- parObj holds parameters for calibration purposes
%
% It has fields name, value, LB, UB, logFlag useFlag
% where LB is lower boundary, UB upper boundary logFlag a flag indicating to apply
% ln transformation when calibrating and useFlag a flag indicating to incclude
% this parameter in the caliration.
%
% This parameter object is useful in conjunction with calibration using the
% non linear least squares method of Matlab lsqnonlin.
% It has several methods as outlined below.
% TO 130621
    properties
        name    = '';
        value   =[];
        LB      = -Inf;
        UB      = +Inf;
        logFlag = true;
        useFlag = true;
    end
    methods
        function o = parObj(varargin)
            %%PAROBJ/PAROBJ -- constructor for parObj objects
            % parObj objects hold parameters, their name, value, lower bound,
            % upper bound, logFlag and useFlag (see classdef above).
            %
            % USAGE
            %    par = parObj(inCell)
            %     par = parObj
            %     inCell = {parNm1,parVal1,LB,UB,LogFlag,UseFlag;...
            %               parNm2,parVal1,....};
            %     parNm,ParVal and UseFlag are obligatory
            %     interpresation of inCell depends on size(inCell,2)
            %     case 3   parNm,parValue,useFlag
            %     case 4   parNm,parValue,-Inf Inf,logFlag,useFlag
            %     case 5   parNm,parValue,-Inf,upperBound,logFlag,useFlag
            %     case 6   parNm,parValue,lowerBound upperBound,logFlag,useFlag
            % TO 130621
            
            if nargin<1
                return;
            end
            
            [inCell,varargin] = getNext(varargin,'cell',[]);
            if isempty(inCell)
                [basename,varargin] = getNext(varargin,'char',[]);
                [sheetNm , ~      ] = getNext(varargin,'char',[]);
                if isempty(basename) || isempty(sheetNm)
                    error('%s: need basename and sheetNm',mfilename);
                end
            end
            
            if ~isempty(inCell)
            
                w = size(inCell,2);

                err = {'parNm must be of type char',...
                   'value must be a numeric scalar',...
                   'useFlag must be a numeric scalar or a logical',...
                   'logFlag must be a numeric scalar or a logical',...
                   'UB  must be a numeric scalar',...
                   'LB  must be a numeric scalar'};

                if ~iscell(inCell)
                    error('input must be a cell array see help parObj');
                end
                if size(inCell,2)<2
                    error('input cell array must have at least 2 columns see help parObj');
                end

                nPar = size(inCell,1);
                o(nPar).name = 'dummy';

                if ~all(cellfun(@ischar,inCell(:,1))), error(err{1}); end
                [o.name] = deal(inCell{:,1});

                if ~all(cellfun(@isscalar,inCell(:,2))), error(err{2}); end
                for io=1:numel(o)
                    o(io).value = inCell{io,2};
                end
                if w>2
                    if ~(all(cellfun(@islogical,inCell(:,end))) || all(cellfun(@isscalar,inCell(:,end))))
                        error(err{3});
                    end
                    for io=1:numel(o)
                       o(io).useFlag = inCell{io,end};
                    end
                end
                if w>3
                    if ~(all(cellfun(@islogical,inCell(:,end))) || all(cellfun(@isscalar,inCell(:,end))))
                        error(err{4});
                    end
                    for io=1:numel(o)
                        o(io).logFlag = inCell{io,end-1};
                    end
                end
                if w>4
                    if ~all(cellfun(@isscalar,inCell(:,end)))
                        error(err{5});
                    end
                    for io=1:numel(o)
                        o(io).UB  = inCell{io,end-2};
                    end                
                end
                if w>5
                    if ~all(cellfun(@isscalar,inCell(:,end)))
                        error(o.err{6});
                    end
                    for io=1:numel(o)
                        o(io).LB = inCell{io,end-3};
                        if o(io).LB>o(io).UB
                            p=o(io).LB; o(io).LB=o(io).UB; o(io).UB=p;
                        end
                    end
                end
            else       
                [initHdr,initVals,initTxtHdr,initTxt] = getExcelData(basename,sheetNm,'Hor');
                for io = size(initVals,1):-1:1
                    o(io).name    = initTxt{ io,strmatchi({'name','par'}   ,initTxtHdr)};
                    o(io).useFlag = initVals(io,strmatchi('useFlag',initHdr));
                    o(io).logFlag = initVals(io,strmatchi('logFlag',initHdr));
                    o(io).value   = initVals(io,strmatchi('value'  ,initHdr));
                    o(io).LB      = initVals(io,strmatchi('LB'     ,initHdr));
                    o(io).UB      = initVals(io,strmatchi('UB'     ,initHdr));
                end
            end
        end
        function [P,lb,ub,pars,o] = Par2p(o,p)
            %%PAROBJ/PAR2P -- get P to be calibrated together with LB,UB and pars.
            % optionally also the parNew can be computed but is not
            % required in any calibration as Par is not changed during
            % the calibration, only p is.
            %
            % This function prepares for calibration a model by non-linear optimization
            % routines of Matlalb, such as p = lsqnonlin(@fun,p0);
            % Make sure that Par is global in fun and the workspace, so that all
            % parameters are known to fun. In fun reassemble the set using also the
            % parameters passed by lsqnonlin.

            % USAGE:
            %   [p,lb,ub,pars,newPar] = Par.toCalib(p) % with p
            %   [p,lb,ub,pars,newPar] = Par.toCalib()  % without p
            %    P = list of ln transformed parameters to be calibrated
            %    lb= lower bound of transformed parameter
            %    ub= upper bound of transformed parameter
            %    pars = {parNm1,parVal1,parNm2,parVal2,...} all parameters
            %    parNew is of class parObj, the p-converted parameters.
            %    Par is the original parameter set.
            %
            %    P,lb,ub and p all pertain to the parameters that are to be
            %    calibrated
            %    pars, parNew and Par to all parameters that are not in te
            %    defaults set. I.e. also those with Par.useFlag = false.
            %    p is vector of parameters to calibrate

            % get parameters
            
            if  nargin<2
                P = NaN(size(o(:)));  lb= NaN(size(P));  ub= NaN(size(P));

                if nargin<2
                    p=NaN(size(o));
                    j=0;
                    for io=1:numel(o)
                        if o(io).useFlag
                            j=j+1;
                            if o(io).logFlag
                                p(j)=0;
                            else
                                p(j)=o(io).value;
                            end
                        end
                    end
                    p= p(~isnan(p));
                end
                if size(p(:)) ~= size(P)
                    error('size(p) does not match size(find([Par.useFlag]))');
                end
            end
            
            pars{1,2*numel(o)}='';
            
            j=0;
            for io=1:numel(o)
                pars{2*io-1} = o(io).name;   % allocate
                if o(io).useFlag % if parameter needs to be calibrated
                    j=j+1;
                    if o(io).logFlag
                        % multiplyer, logFlag par, initially 1
                        P( io) = p(j);     % log(multiplier) i.e. logFlag(Par{ip,2}/Par{ip,2}); 
                        lb(io) = log(exp(p(j))*o(io).LB/o(io).value);
                        ub(io) = log(exp(p(j))*o(io).UB/o(io).value);
                        o(io).value = exp(p(j))*o(io).value;
                    else % linear parameter
                        lb(io)      = o(io).LB;
                        ub(io)      = o(io).UB;
                        P( io)      = p(j);      % adds to parameter
                        o(io).value = p(j);
                    end
                end
                pars{2*io} = o(io).value;
            end
            P( isnan( P))=[];
            lb(isnan(lb))=[];
            ub(isnan(ub))=[];
        end
        function o = changeBy(o,changeBy)
            %%PAROBJ/CHANGEBY -- multiply the paramters with the values in changeBy for as far
            %              as Par.useFlag is true or nonzero
            %
            % This function changes the parameter values in the parObj by using
            % the values in parObj changeBy.
            % parameter values with useFlag=0 or false are not changed
            % parameter values with logFlag=1 or true  are multiplied by the
            % values in changeBy
            % parameter values with logFlag=0 or false get their old value
            % polus the value in changeBy.
            %
            % USAGE:
            %   par = par.changeBy(changeBy)
            %   changeBy is also a parObj, its values are useFlag to change
            %   those in par depending on the useFlag and the logFlag in par
            %   as well as on the LB and UB in par.
            %   The useFlag and logFlag in changeBy are ignored.
            %   Only its pararmeter names of the parmates in changeBy will be
            %   checked against that of par.
            %        and its value will be used to change Par if useFlag is on
            %        if ln is on  Par is multiplied by changeBy
            %        if ln is off changeBy is added to Par
            %
            %    Par should be pase
            % TO 130620

            chNames = {changeBy.name};

            for io=1:numel(o)
                if o(io).useFlag
                    i = strmatchi(o(io).name,chNames); i=i(1);
                    if ~i
                        error('parName %s in par not in changeBy <<%s>>',...
                            o(io).parName,sprinfs(chNames));
                    end
                    if o(io).logFlag
                        o(io).value = o(io).value * changeBy(i).value;
                    else
                        o(io).value = o(io).value + changeBy(i).value;
                    end
                    if o(io).value<o(io).LB, o(io).value = o(io).LB; end
                    if o(io).value>o(io).UB, o(io).value = o(io).UB; end
                end
            end
        end
        function L = equal(o,Par)
            %%PAROBJ/EQUAL -- see if the two Par objects are equal
            %
            % USAGE:
            %   logical = Par.equal(Par2);
            %
            % TO130625
            
            tol = 1e-6; 
            L = numel(o) == numel(Par) ...
                && all(strcmpi({o.name},{Par.name})) ...
                && all(abs([o.value] - [Par.value])<tol) ...
                && all([o.logFlag] == [Par.logFlag]) ...
                && all([o.useFlag] == [Par.useFlag]);
        end
        function L = checkEquality(o,Par)
            %%PAROBJ/CHECKEQUALITY -- checks for equality of the two parObj
            %
            % USAGE:
            %     logical = Par1.checkEquality(Par2);
            %
            % TO 130525
            tol = 1e-6; 
            if  numel(o) ~= numel(Par)
                error('%s: number of parameters differs',mfilename);
            elseif ~all( strcmpi({o.name},{Par.name}) )
                error('%s: parameter names differ',mfilename);
            elseif any(abs([o.value] - [Par.value])>tol)
                error('%s: values differ',mfilename);
            elseif any(o.logFlag ~= Par.logFlag)
                error('%s: logFlags differ',mfilename)
            elseif any(o.useFlag ~= Par.useFlag)
                error('%s: useFlags differ',mfilename);
            else
                L = true;
            end
        end
        
        function display(o)
            %%PAROBJ/DISPLAY -- shows paameter object
            %
            % USAGE:
            %   Par.display();
            %
            % TO 130625
            fprintf('%12s','Parameter','value','LB','UB','logFlag','useFlag');
            fprintf('\n');            
            for io = 1:numel(o)
                fprintf('%12s%12.4g%12.4g%12.4g%12d%12d',...
                    o(io).name,...
                    o(io).value,...
                    o(io).LB,...
                    o(io).UB,...
                    o(io).logFlag,...
                    o(io).useFlag);
                fprintf('\n');
            end            
        end
        
        function show(o,cases,Pars,stdP)
            %%PAROBJ/SHOW -- compars parametes Par with ParNew
            %
            % Puts next to eachother parameter Par and ParNew
            % also adds stdPar and uncertainty of parameters
            %
            % USAGE:
            %   Par.show(cases,{Par1,Par2,Par3,...},stdP);
            %    Par_i is parObj of same size and same order as Par
            %    useFlag and logFlag of Par will be used and of ParNew will
            %    be ignored. A warning will be issued if they differ.
            %    stdP is a vector of paramter standar deviations of length
            %    P, the calibrated parameter set. This length and order
            %    must be as that of [Par.useFlag]==true;
            %    uncert = 100*exp(stdP)/exp(P) if logFlag is true
            %    uncert = stdP if logFlag is false
            %
            % TO 130620

            msgId='parObj:show:logFlagOrUseFlag';
            warning('on',msgId);
            
            fprintf('%12s%8s%8s','Parameter','useFlag','logFlag');
            fprintfs('%12s',cases);
            fprintf('%12s','stdP');
            fprintf('\n');

            for ip = numel(Pars)
                if ~all(cellfun(@strcmpi,{o.name},{Pars{ip}.name}))
                    error('names in Par must match those in newPar exactly');
                end
                if ~all([o.useFlag]==[Pars{ip}.useFlag])
                    warning(msgId,'not all Par.useFlag==newPar.useFlag');
                end
                if ~all([o.logFlag] ==[Pars{ip}.logFlag])
                    warning(mgsId,'not all [Par.logFlag]==[newPar.logFlag]');
                end
            end
            
            j=0;
            for io = 1:numel(o)
                fprintf('%12s%8d%8d',o(io).name,o(io).useFlag,o(io).logFlag);
                fprintf('%12.4g',      o(io).value);
                for ip=1:numel(Pars)
                    fprintf('%12.4g', Pars{ip}(io).value);
                end
                if o(io).useFlag
                    j=j+1;
                    fprintf('%12.3g',stdP(j));
                end
                fprintf('\n');
            end            
            warning('off',msgId);
        end
        
        
        function meas = syntheticMeas(o,NP,stdMeas)
            %%PAROBJ/SYNTHETICMEAS -- saves new synthetic measurements as
            % a vector of length NP.
            % Also adds random noise with std = stdMeas
            % o must bue the true parameter set, i.e. the set for which the
            % model is supposed to provide the synthetic true values, i.e. the
            % measurements
            %
            % USAGE:
            %   meas = Par.meas(NP);
            %   measv= Par.meas(relPos)
            %
            % NP is number of points
            % relPos is vector of relative points locations [0..1]
            %
            % TO 130621

            global defaults

            NP = NP(:)';
            if ~isnumeric(NP)
                error('NP must be a numeric scalar');
            end

            L = getProp(defaults,'L',[]);

            if numel(NP)==1
                x  = L * rand(1,NP);
            else
                x  = L * NP;
            end

            [~,~,~,pars] = o.Par2p();

            meas.y = model('x',x,pars{:},defaults{:}) + stdMeas * randn(size(x(:)));
            meas.x = x(:);
        end
        
        function o = p2Par(o,p)
            %%PAROBJ/P2PAR -- converts paramter vector to parObj
            %
            % USAGE:
            %     Par2 = Par1.p2Par(p);
            %
            % TO  130625
            j=0;
            for io=1:numel(o)
                if o(io).useFlag
                    j=j+1;
                    if o(io).logFlag
                        o(io).value = exp(p(j)) * o(io).value;
                    else
                        o(io).value = p(j);
                    end
                end
            end
        end
        function sigP = uncert(o,sigP)
            %%PAROBJ/UNCERT -- computes uncertainty of the inidividual parameters
            %
            % USAGE:
            %   unc = ParNew.uncert(sigP);
            %   sigP is standard dev as obtained from lsqnonlin
            %     sigP = Cov./(diag(Cov)*diag(Cov)')
            %     Cov = sigE^2 [J'J]^(-1)
            %  if logFlag == true --> uncert is 100*sigP(j)/ParNew(io).value
            %  else                   uncert = sigP as is.
            %
            % TO 130621
            
            j=0;
            for io=1:numel(o)
                if o(io).useFlag
                    j=j+1;
                    if o(io).logFlag
                        sigP(j) = 0.5*(exp(sigP(j))-1)+(1-exp(-sigP(j)))*o(io).value;
                    end
                end
            end
        end
    end
end