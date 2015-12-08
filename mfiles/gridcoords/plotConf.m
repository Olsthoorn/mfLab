function ax=plotConf(varargin)
%PLOTCONF plots configuration specified in sheets Config and Materials
%
% Example:
%    ax=plotConf([ax,]Conf,codes,color,alpha)
%    (see the examples that use it. (swt_v4>freshkeeper))
%
% Example
%    plotConf(Conf,codes,color,alpha)
%
% The text below is outdated, should be replaced by obj xsConf
%   The struct conf holds info of the conductivity configuration in a cross
%   section. This configuration is a rectangular range with material codes
%   like 'K', 'S' etc. while also the xL, xR, top and bot of the blocks are
%   give. This fucntion allows to patch/fill the material blocks of specified
%   material codes, color and alpha (transparency) value
%   Conf = output of mf_conf
%   codes are material codes for instance {'K' 'S' 'P'}
%   color are colors as string of ['r ' 'b' 'k' etc] or a N times 3 array
%   with RGB values
%   alpha is a ist of transparency values
%   missing color and alpha values are padded, too many are ignored
%
% Used in
%     mflab/examples/Geohydrology2/RemkoNijzink
%     mflab/examples/mf2k/qanat
%     mflab/examples/swt_v4/ATES-WKO/Euronext-Amsterdam
%     mflab/examples/swt_v4/BanglaDesh
%     mflab/examples/swt_v4/BeachDrainage
%     mflab/examples/swt_v4/Freshwater-lens
%
% ToDo: check compatibility with xsConfObj, replace by method
%
% See also:  mf_plotConf xsConfObj
%
%   TO 110320 120418

if ishandle(varargin{1}),
    ax=varargin{1};
    varargin=varargin(2:end);
else
    ax=gca;
end

nargin=size(varargin,2);

if nargin>=1, Conf=varargin{1};  end
if nargin>=2, codes=varargin{2}; end
if nargin>=3, color=varargin{3}; else color='b'; end
if nargin>=4, alpha=varargin{4}; else alpha=1;   end


[Ny,Nx]=size(Conf.mzone);

if isfield(Conf,'xL'),
    xL=Conf.xL; xR=Conf.xR;
else
    xL=Conf.yL; xR=Conf.yR;
end
bot=Conf.bot*ones(1,Nx);
top=ones(Ny,1)*Conf.top; top(2:end,:)=bot(1:end-1,:);

if ~iscell(codes), codes={codes}; end

for i=1:length(codes)
    
    cd=codes{i};
    
    for iy=1:size(Conf.mzone,1)
        IX=strmatchi(cd,Conf.mzone(iy,:),'noerr');
        if any(IX)>0
            for j=1:length(IX)
                h=fill([xL( IX(j)),xR( IX(j)),xR( IX(j)),xL(IX(j))],...
                        [bot(iy   ),bot(iy   ),top(iy   ),top(iy  )],...
                        color(min(i,size(color,1)),:),'parent',ax);
                set(h,'facealpha',alpha(min(i,length(alpha))),...
                    'edgecolor','none');
            end
        end
    end
end
