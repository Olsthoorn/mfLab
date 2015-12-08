function [BCN PNTSRC]=bcnparser(type,zonearray,zones,values)
%% Stump trial TO 120410
% under development
% does not work as yet
%

for iz=1:length(zones)
    Idx=find(zonearray==zones(iz));
    if ~isempty(Idx)
        LRC=cellindics(Idx,size(IBOUND));
        u=ones(size(LRC,1));
        for it=1:NT
            BCN=[BCN;...
                it*u LCR u*[h1 h2] CHDOPTT];
            %  if values are numeric
            % otherwise, they are headers of columns in PER
        end
    end
end