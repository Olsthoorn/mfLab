classdef particleStartingLocationsObj
    % definition of automatiically generated starting locations for Modpath
    properties
        grid=1;
        gridCellRegionOption
        placementOption
        releaseStartTime
        releaseOption
        CHeadOption
        releasePeriodLength
        releaseEventCount
        minLay,maxLay,minRow,maxRow,minCol,maxCol
        mask
        masklayer
        faceCount
        IFace
        particleLayCount
        particleRowCount
        particleColCount
        particleTimeCount
        releaseTimeIncrement
        timePoints
    end
    properties(Dependent=true)
    end
    methods
        function o=particleStartingLocationsObj(options)
            if nargin==0, return; end
        end
        function write(o,fid)
            fprintf(fid,'%d\n',length(o)); % groupCount
            switch options.particleGenerationOption
                case 1
                    for i=1:length(o)
                        fprintf(fid,'%s\n',o(i).groupName);
                        fprintf(fid,'%d %d %d %g %d %d\n',...
                            o(i).grid,o(i).gridCellRegionOption,o(i).placementOption,...
                            o(i).releaseStartTime,o(i).releaseOption,o(i).CHeadOption);
                        if options.releaseOption==2
                            fprintf(fid,'%g %d\n',...
                                o(i).releasePeriodLength,o(i).releaseEventCount);
                        end
                        switch o(i).gridCellRegionOption
                            case 1
                                fprintf(fid,'%d %d %d %d %d %d\n',...
                                    o(i).minLay,o(i).minRow,o(i).minCol,...
                                    o(i).maxLay,o(i).maxRow,o(i).maxCol);
                            case 2
                                warray(fid);
                            case 3
                                fprintf(fid,'%d\n',o(i).maskLayer);
                                warray(fid,layer as given)
                        end
                        switch o(i).placementOption
                            case 1
                                fprintf(fid,'%d\n',faceCount)
                                for ie=1:faceCount
                                    fprintf(fid,'%d %d %d\n',...
                                        o(i).IFace,o(i).particleRowCount,o(i).particleColCount);
                                end
                            case 2
                                fprintf(fid,'%d %d %d\n',...
                                    o(i).particlLayCount,o(i).particleRowCount,o(i)particleColCount);
                        end
                    end
                case 2
                    fprintf(fid,'%s\n',startloc);
            end
        end
    end
end