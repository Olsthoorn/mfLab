function o=getConc(o,varargin)
% USAGE:
%     obs = obs.getConc();
%     obs = obs.getConc(STRNG1, STRNG2, STRNG3 ..)
%     obs = obs.getConc(STRING, iComp)
% sets the computed obs concentration. It is the flow-
% averaged concentration at the obs's model cells during the simulation.
% C, C1, C2 ... are structs read by readMT3D.m.
%
% Dissolved concentrations will be stored as obs.C(iComp,1:NT)
% sorbed concentrations will be in obs.CS(iComp,1:NT);
%
% SEE ALSO: wellObj.setCout
%
% TO 120423 121017 121129

%     obs = obs.getConc([iComp iComp2 ...])

if numel(varargin)==0
    compoundNames = o(1).compound;
    if isempty(compoundNames)
        error('%s: no compoundNames specified, see help cObsObj',mfilename);
    end
else
    if ischar(varargin{1})
        compoundNames = varargin(1);
    elseif iscell(varargin{1})
        compoundNames = varargin{1};
    else
        error('%s: specify compound names as a cell array of strings not as <<%s>>',...
            mfilename,class(varargin{1}));
    end
end

for iComp = numel(compoundNames):-1:1
    if ischar(compoundNames{iComp})
        C(iComp).solved = readMT3D(sprintf('MT3D%03d.UCN' ,iComp));
        try
            C(iComp).sorbed = readMT3D(sprintf('MT3D%03dS.UCN',iComp));
        catch  %#ok
        end

        time = [C(iComp).solved.time];
        NT   = numel(time);
    end
end


% all observation points
for iObs = numel(o):-1:1    
    o(iObs).t     = time;
    o(iObs).Dt    = diff(time);
    o(iObs).NCOMP = numel(compoundNames);
    o(iObs).compound = compoundNames;
    for iComp = numel(compoundNames):-1:1
        for it=NT:-1:1
            o(iObs).C(iComp,it)             = C(iComp).solved(it).values(o(iObs).idx(1));
            try
                o(iObs).C(iComp).sorbed(it) = C(iComp).sorbed(it).values(o(iObs).idx(1));
            catch %#ok
            end
        end
    end
end