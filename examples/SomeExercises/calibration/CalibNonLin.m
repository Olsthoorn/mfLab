%% Generic calibration using Matlab's lsqnonlin optimizer
%  Calibration can be done in Matlab using one of its optimizers. The most
%  useful one for groundwater models should be the Marquardt Levenberg
%  solver embeddedin the Matlab function lsqnonlin (non-linear least squares).
%  This function has to be called as follows
%
% [pEnd,....J] = lsqnonlin(@modelWrapper,pStart)
%
% The input is a pointer to the cumtom-made function modelWrapper. modelWrapper
% accepts a parameter vector as input. It provides the the model results
% minus the measurments (or vice versa) as output. lsqnonlin will then
% internally iterrate until the residuals are minimal.It finally provides
% the vector with the so-optmized parameters.

%% Pre and after processing is done in a script called
%
%    CalibNonLin
%
% CalibNonLin calls lsqnonlin with the proper parameter set. CalibNonLin
% also presents the results.

%% Description of the overall process
% The actual model takes more than just the parameters p. These other parameters
% and anything else that the model needs except for the parameters that it
% uses in the calibration, is made available through a global declaration
% both in the workspace and in the modelWrapper.
% The modelWrapper can then assemble the actual and
% complete parameter set required and pass this to the actual
% groundwater model. Parameters may also be converted to their log
% value in the modelWrapper. The modelWrapper also transfers the outcomes
% of the regular model into a set of residuals. This is done by selecting the
% required data from the output of the actual model and subtracting the
% measurements. Note that also these measurements are passed to the
% modelWrapper through a global declaration in both the workspace and the
% modelWrapper. Alternatively, one might prefer reading the measurements
% form disk at any call of the modelWrapper, but this will be much much
% slower.

%% Structure of the calibration
% The structure is then as follows
%
% calib
%    makes sure that measurements are present
%    prepares parameter vector p0
%    calls lsqnonlin(...p0)
%      loop
%         lsqnonlin --> p modelWrapper
%             prepares actual parameter set for model
%             actual input -->  model
%             modelwrapper <--- groundwater results <--- model
%             selects results at measurement locations and times
%             computes residuals by subtracting measurements
%         lsqnonlin <-- residual <--- modelwrapper
%         sees if optimum has been reached
%         if not loop once more
%         else finish loop
%      end loop
%    caib <-- final parameter vecter, statistics <-- lsqnonlin()
%    compute various paramter statistics
%    display output
%      Compare oldParametrs with New parameters
%      Shows diverse statistics
%      plots the measurements, initial model and final model results
%  enb calib

%% Parameter transformation, lower and upper limits
%  It is generally preferred to optimize the log of the actual parameters
%  or rather the log of a multiplyer for the actual parameters. This is
%  beneficial for all parameters that have zero as a natural lowest limit.
%  Hence the parameter is then
%        par = exp(p)*par0
%        p   = ln(par/par0)
%  p=0 is then a good starting value for the parameters.
% This transformation works less well for parameters without a natural zero
% such as head or elevation. In that case we may have to transform the
% actual values such that they always remain positive. Depth may thus be a
% better choice than elevation. We can always transform such that the
% multiplyer remains positive. For instance by adding a suitable value.
% However, the value range of such parameters will not be between 0 and Inf
% as is the case with many physical quantities. Therefore we should in that
% case use lower and upper limits in the call to lsqnonlin
%    [pEnd,....] = lsqnonlin(p0,ll,ul)
% For instance if alpha is a parameter that can vary between 0 and 1
% and we are using p = ln(alpha/alpha0) then we should also use
%    ll = ln(lowestAlpha /alpha)
%    ul = ln(highestAlpha/alpha)
% To make this work, we see that we should set lowestAlpha to a small
% value>0 instead of zero. Perhaps we could add 1 to alpha and let the
% actual parameter vary between 0 and 1 or -1 and 1 for instance.
%
% modelWrapper should take care of all the necessary to and back parameter
% conversions and is, therefore, alsways custom made.
% instead of modelWrapper(), we may prefer calling this function residuals()
% as that is what it computes. The results are only seen within lsqnonlin
%%   
% At the level of calib only actual parameters will be of interest
% However, we have to pass the converted paramters, which is a log-transformed
% subset of the total set of parameters to lsqnonlin. The length of p is
% equal to the parameters that are being calibrated only.
% We also have to deal with the final pEnd outcome of lsqnonlin. These have
% to be converted to actual parameter values for only those parameters that
% were optmized in the calibration.

%% Defining and passing parameters
% To allow this we will work with a parameter set contained in a cell array
% that has the following fields defining one parameter per line
% Par =
% { parName parValue lowerBound upperBound use; ...
%   parName ....
% }
% The use value, which is 0 or nozero or true or false, is used to dicern
% which of the parameters will take part in the calibratin and which will
% be passed to the model as given.
% The model will generally need many more parameters than contained in the
% parameter set just outline. Such other parameters will never be calibrated
% will be calibrated and are called defaults. The default parametes are
% communitated with the modelWrapper through a global declaration in both
% the workspace and in modelWrapper.
% Defaults are specifie das parametername,value pairs. They are contained
% in a cell array consisting of one line:
% defaults = {parName1,parValue1,parName1,parvalue2,....}
%%
% It is a user's choice which parameters will be considere defaults and which
% are conveyed through the Par cell array just explained. Generally it is
% not convenient to have a very long Par list with only a few of the
% defined parameters actually calibrated.
%%
% modelWrapper will take care to combine the defaults with the calibration
% subet. It will apply the parameter set p, i.e. the ln values of the
% parameters to be calibrated on the parameters in Par and leaves the
% parameters in Par that will not be calibrated untouched.
% Then it will join the parameters and the defaults as a cell vector
% of parName,parValue to obtain the full set of parameters necesssary
% to run the model and pass this vector to the model.
% The model will interpret this input as parName,parValue pairs, which can
% be specified in any order as long as the set is complete.

%% Results and calibration statistics
% The final parameters will be converted to actual ones in calib.
% The old and new parametrs will be printed together with their esimated
% standard deviation and uncertainty, which is defined as
% 100*stdPar/abs(Par)
% The parameter covarians and correlation matrices will be computed and
% shown and the measurements will be shown together with the initial and
% final model results.
% TO 130619


%% Calibration using non linear parameter optimizationUsing lsqnonlin

global Par defaults meas iter % parameters that must be statically known in modelWrapper

iter = 0;

%% Defaults: initial parameter values not changed during any optimization
defaults = { 'sL',  0, 'sR', 0, 'L', 2000};

%% Initial parameters of the model
% Initial parameters that may be calibrated depending on useFlag
% log transformation according to log flag.
% Remember Par are original parameters, and are not changed during the
% calibration. Only the p vector is updated with which the modelWrapper and
% lsqnonlin is called.

% To change the initial parameters, change the values in the second column
% To change the lower bound of a parameter change the value in column 3
% TO change the upper bound of a parameter change the value in column 4
% To change the logFlag, set value in column 5 to either 0 or 1
%    values with logFlag=1 are transformed to their ln value during the
%    calibration.
% To change the useFlag, set value in column 6 to either 0 or 1
%    values with useFlag=0 are not calibrated but kept constant.
Par = parObj({  % initial parameters
    % NAME VALUE LB  UB  logFlag useFlag
    'N'      0.001  0.0001 0.5 1 0
    'kL'     100    0.1  500  1  1
    'kR'      10    0.1  500  1  1
    'alpha'  0.5    0      1  0 1   
    'zB'     -10   -5  -40  0 0 
    });


%% Generate synthetic measurements if necessary

% The true parameters are given here. They are ued to generate
% measurements. The measurements must thus be consistent with the true
% parameters. The measurements are generated by them and then saved on
% disk. The true parameters never change unless by hand, after which new
% measurements must be computed.
% Some noise is added to the measurements as given by the passed standard
% deviation for the heads.

% to update the true parameters and to generate new measurements that
% comply with the true parameters, delete meas.mat from the directory and
% run this mfile again.
parTrue = parObj({
        'N'       0.001
        'kL'         10
        'kR'         10
        'alpha'     0.5
        'zB'        -10
        });

stdMeas =  0.01; %0.04; % standard error of the measurements
NP      =  15; % number of head measurement locations
NP = (0:100)/100; NP([1 end])=[];
% If loading measures is succuessful, create them
if ~exist('meas.mat','file')
    truePar = parTrue;
    meas = parTrue.syntheticMeas(NP,stdMeas);
    save meas meas truePar;
else
    load meas
    if ~truePar.equal(parTrue) && length(meas.y)~=length(NP)
        truePar = parTrue;
        meas = parTrue.syntheticMeas(NP,stdMeas);
        save meas meas truePar
    end
end

%% Use of the non-linear solver
[p0,lb,ub] = Par.Par2p(); % get initial parameters

[pEnd,resNorm,reSid,exitflag,output,lambda,J] = lsqnonlin(@modelWrapper,p0,lb,ub);

parEnd = Par.p2Par(pEnd);

%% Finalizing statistics derived from the Jacobian (sentivity matrix J)

J=full(J);  % convert sparce J to full J. J is not sparse

Inv    = (J'*J)^(-1); B   = Inv*J';

[U,S,V] = svd(J,0);

display(S);
display(V);

% Both columns should be the same
fprintf('     pEnd      pEnd\n');
display([pEnd B*meas.y]);       % final parameters from lsqnonlin and J

Cov    = resNorm*Inv;           % covariance matrix of the parameters

sigmaP = sqrt(diag(Cov));       % std of the parameters

Cor    = Cov./(sigmaP*sigmaP'); % correlation matrix of the parameters

display(Cov);
display(Cor);

%% Uncertainty reporting
[~,~,~,~,parNew] = Par.Par2p(pEnd);

sigP = Par.uncert(sigmaP);    % uncertainty

Par.show({'oldPar','truePar','newPar'},{parTrue,parNew},sigP);

figure; hold on;
xlabel('x [m]'); ylabel('head [m]'); title('calibration');

fprintf('\n\n');

%% Plot results
plot(meas.x,                    meas.y, 'ro');
plot(meas.x, modelWrapper(p0  )+meas.y, 'xb');
plot(meas.x, modelWrapper(pEnd)+meas.y, 'gs');

legend('meas','init','end');