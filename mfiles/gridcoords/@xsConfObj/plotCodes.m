function h=plotCodes(o,varargin)
% ax=xsConf.plotCodes(ax,materialNames [ RGB|colorspec [,alpha], varargin)
% plotCodes: colors the materials in configuration specified
% in sheets Config and Materials
%
% materialNames can be specified as material strings
% either a character string or a cell array of strings.
% RGB is an array with one or more lines or RGB (as in colormap)
% alternatiely, use colorspec like 'w' 'g' or 'wgbc'
% if not specified, red green and blue columns in sheet Materials will be
% used.
% If there are more materials than colors, the last will be used for the
% remaining materials.
% Same for alpha (transparency). The default is 0.35
%
% USAGE
%   conf.plotCodes(o,codes[,'color',color][,'alpha',alpha])
%
%   The objects xsConfObj holds info of the conductivity configuration in a cross
%   section. This configuration is a rectangular range with material codes
%   like 'K', 'S' etc. while also the xL, xR, top and bot of the blocks are
%   give.
%   This method patches the material blocks occupied by the specified
%   material codes.
%   'color' and 'alpha' are options that allow choosing a color and
%   transparency.
%
% SEE ALSO: xsConfObj
%
%   TO 110320 120418 130409

[ax    ,varargin] = getNext(varargin,'axis',gca);
[codes ,varargin] = getNext(varargin,'char',[]);
if isempty(codes)
    error('%s: No matarial codes specified for plotting',mfilename);
end

[faceColors ,varargin] = getNext(varargin,{'char','double'},'w');
[alpha  ,varargin] = getNext(varargin,{'double'},0.25);

if ~iscell(codes), codes={codes}; end

for ic=1:length(codes)
    
    code=codes{ic};

    % zoneCodes is a cell matrices with the codes as strings
    % iMat is a numerical array the material numbers corresponding
    % to the materials in the material list.

    idx = find(ismember(o.zoneCodes,code));

    if isempty(idx), continue; end
    
    if ~isempty(faceColors)
        if ischar(faceColors)
            color = faceColors(min(ic,numel(faceColors)));
        else
            if ~isnumeric(faceColors) || size(faceColors,2) ~= 3
                error(['%s: colors specified as second argument is not a legal color(s)\n',...
                       'must be a colorspec or a nx3 array where each line is a legal RGB vector.'],...
                        mfilename);
            end
            m = min(ic,numel(size(faceColors,1)));
            color = faceColors(m,:);
        end
    else
        r = matPropValues(iMat,strmatchi('red'  ,o.matProps));
        g = matPropValues(iMat,strmatchi('green',o.matProps));
        b = matPropValues(iMat,strmatchi('blue' ,o.matProps));
        color = [r g b];
    end


    % which config cells contain this material ?
    RCL = cellIndices(idx,size(o.zoneImat),'RCL');
    
    % index in config cross section
    IY  = RCL(:,1);
    IX  = RCL(:,2);
    

    % patch all rectangles of this material
    for i=1:numel(idx)
        
        h=fill([ o.xL(IX(i)), o.xR(IX(i)), o.xR(IX(i)), o.xL(IX(i))],...
               [ o.zoneZ(IY(i)+1,IX(i)), ...
                 o.zoneZ(IY(i)+1,IX(i)), ...
                 o.zoneZ(IY(i)  ,IX(i)),...
                 o.zoneZ(IY(i)  ,IX(i))],...
                 color,'edgeColor',color,'faceAlpha',alpha,'edgeAlpha',alpha,'parent',ax,varargin{:});
    end
end
