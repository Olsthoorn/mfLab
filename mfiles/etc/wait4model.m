function wait4model(mdlname)
%WAIT4MODEL wait until model that runs in background has finished
%
% USAGE:
%    mf_wait(mdlname);
%
% Waits until model finishes on system after call SYSTEM('model
% inputfile').
%
% TO 120421

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

    for i = 1:100
        [~,result]=system(['ps | grep ',mdlname]);
        fprintf('Waiting for model %s to finish...\n',mdlname);
        if length(strfind(result,mdlname))<3
            fprintf('... ok, ready, continuing after %d seconds.\n',round(toc));
            break;
        end
        pause(1);
    end
