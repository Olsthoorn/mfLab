%% Run verification case
tstart=TPE(  1,1);
tend  =TPE(end,1);
I= TPE(:,1)>=tstart & TPE(:,1)<=tend;
str = sprintf('Van %s tot %s, %.0f d: mean P E P-E = %.4f %.4f %.4f, sum P E P-E = %.4f %.4f %.4f\n',...
    datestr(tstart),datestr(tend),tend-tstart,...
    mean(TPE(I,2)),mean(TPE(I,3)),mean(TPE(I,2)-TPE(I,3)),...
    sum( TPE(I,2)),sum( TPE(I,3)),sum( TPE(I,2)-TPE(I,3)));
fprintf(str);

[h,Qfix,Qnod,Qsto,QCau,Qz,z,bodem,crop,SS,Theta,Kh,Acrop,Intercepted]=unsat1(TPE,'ClapHornberger78SoilData',sheetNm);
