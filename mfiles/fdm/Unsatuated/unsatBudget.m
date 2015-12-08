%% Check water budget

load unsat1.mat

% Verify the storage coefficient
fprintf(' %11.4g',sum(SS .* (dz*ones(1,size(SS,2))))); fprintf \n
fprintf(' %11.4g',sum(diff(Theta,1,2)./diff(h,1,2)));  fprintf \n

if showBudget == true
    % The budget consists of
    % fixed flows P, E
    % Cauchy flows, drainage, exchange with ditch
    % Fixed head flows
    % Storage, Ss dh and delta theta
    % internal flows qv
    
    % water budget excluding the first day because we don't know the theta
    % at the beginning of the first day
    t1 = TPE(  2, 1);
    t2 = TPE(end, 1);
    I = TPE(:,1)>=t1 & TPE(:,1)<=t2;
    P  = TPE(:,2);
    E  = TPE(:,3);
    % interception ?
    % seepage
    
    fprintf \n\n\n
    fprintf('Water budget %s model unsat1\n', bodemCrop);
    fprintf('====================================================\n');
    % data as a function of time
    
    %% showing the water balance terms, one line per time step
    fprintf(' %11s','P-Intercepted','a*(E-Intercepted)','seepage','Qfix','QCau','QstoSS','QstoTheta','Qnod','QfixedHd');
    fprintf \n
    display([(P - Intercepted), ...
        -sum(bsxfun(@times,Acrop,rootZone),1)'/sum(rootZone) .* (E - Intercepted), ...
        seepage * ones(Nt,1), ... 
        sum(Qfix)', ...
        sum(QCau)', ...
        -sum(QstoSS)', ...
        -sum(QstoTheta)', ...
        sum(Qnod(1:end-1,:))', ...
        Qnod(end,:)']);
    
    %% Water budget, comparing fows with storage one line per time step
fprintf('Total sum over all each of the individual water budget terms\n');
fprintf('Qfix+QCau+QfixedHd, QstoSs, QstoTheta\n');
    display([sum(Qfix)'+sum(QCau)'-sum(QstoSS)'-Qpond', sum(Qnod)', sum(QstoTheta)', Qpond', Pond'  h(1,:)']);

%%   Show for comparison QstoSS and QstoTheta
%     figure; plot(TPE(:,1),[-sum(QstoSS); -sum(QstoTheta)]); datetick; grid on;
%     xlabel('2015'); ylabel('totaal uit berging [m/d]');
%     title('vergelijk bergingsberekening QstoSS met QstoTheta');
%     legend('QstoSS','QstoTheta');
    
    
%    display(reshape([-Qz;K5],[Nz-1,Nt*2])); % when continuous inflow

    fprintf \n\n
end
