function well = distributeC(o,well)
% well = wellSeries(distributeQ(well,HK)
% distribute Q of well series over is wells accoriding to the transmissivity
% of the well screens relative to their sum over the wellSeries.
% TO 1200308

if ~(strcmpi(class(well),'wellObj') || strcmpi(class(well),'MNW1Obj') || strcmpi(class(well),'MNW2Obj'))
    error('%s: First argument must be class <<wellObj>> <<MNW1Obj>>or <<MNW2Obj>> not <<%s>>.',...
        mfilename,class(well));
end

wellNrs = [well.nr];

for iws = numel(o):-1:1
    IW = find(ismember(wellNrs,o(iws).children));

    for iw=IW(:)'
        well(iw).C    = o(iws).C;
        if isempty(well(iw).t)
            well(iw).t = o(iws).t;
            well(iw).Dt= o(iws).Dt;
        end
    end
end
