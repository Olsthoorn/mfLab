function well = distributeQ(o,well,HK)
% well = wellSeries(distributeQ(well,HK)
% distribute Q of well series over is wells accoriding to the transmissivity
% of the well screens relative to their sum over the wellSeries.
% TO 1200308

if ~(strcmpi(class(well),'wellObj') || strcmpi(class(well),'MNW1Obj') || strcmpi(class(well),'MNW2Obj'))
    error('%s: First argument must be class <<wellObj>> <<MNW1Obj>>or <<MNW2Obj>> not <<%s>>.',...
        mfilename,class(well));
end

if ~isnumeric(HK)
        error('%s: Second argument must be a full array of the size of the grid with HK or TRAN');
end
    
%% Transmissivity of screen of individual well
wellNrs = [well.nr];    
for iws = numel(o):-1:1
    if isempty(o(iws).children)
        continue;
    end
    IW = find(ismember(wellNrs,o(iws).children));
    TRAN = 0;
    for iw = IW(:)';
        TRAN = TRAN + sum(well(iw).DZ(:)'.*HK(well(iw).idx));
    end
    for iw = IW(:)'
        well(iw).Q = o(iws).Q * sum(well(iw).DZ(:)'.*HK(well(iw).idx))/TRAN;
        well(iw).t = o(iws).t;
        well(iw).Dt= o(iws).Dt;
    end
end
