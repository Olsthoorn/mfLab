function o=getConc(o,animate)
% USAGE:
%     obs = obs.getConc(animate);
% Dissolved concentrations will be stored as obs.(fieldNm).data
% sorbed concentrations will be in obs.CS(iComp,1:NT);
%
% SEE ALSO: wellObj.setCout
%
% TO 120423 121017 121129

%     obs = obs.getConc([iComp iComp2 ...])

%% Observation values can be obtained from an animateObj
    if ~strcmpi(class(animate),'animateObj')
        error('%s: first argument must be an animateObj with loaed with required data',mfilename);
    end

    for io=1:numel(o)
        if ~isempty(animate.H)
            o(io).data.HDS = getObsData(animate.H,o(io).idx);
        end
        if ~isempty(animate.D)
            o(io).data.DDN = getObsData(animate.D,o(io).idx);
        end
        for iComp = 1:numel(animate.titleStr)
            o(io).data.(animate.titleStr{iComp}) = getObsData(animate.UCN{iComp},o(io).idx);
        end
        %% To DO: Cell by Cell flows
    end
end

function data = getObsData(C,idx)
    data = NaN(2,numel(C));
    try
        data(1,:) = [C.totim];
    catch ME
        %fprintf('%s: %s\ntrying field ''time''\n',mfilename,ME.message);        
        data(1,:) = [C.time];
    end
    for i=1:numel(C)
        data(2,i) = C(i).values(idx);
    end
end