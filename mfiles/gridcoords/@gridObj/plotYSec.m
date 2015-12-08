function h = plotYSec(o,varargin)
% plot cross section along y-axis. Same as plotXSec but in other direction.
% see plotXSec for options
%
% TO 120531 130321

% complete permure the x and y and then use the same function
% need to revert directions of axes because y runs down and x runs up,
% which must be reversed upon permuting the grid.
gr = gridObj(o.yGr,o.xGr,permute(o.Z(end:-1:1,end:-1:1,:),[2,1,3]),o.LAYCBD,o.MINDZ,o.AXIAL);

i=strmatchi('lay',varargin);
if i
    varargin{i+1} = permute(varargin{i+1}(end:-1:1,end:-1:1,:),[2,1,3]);
end

i = strmatchi('cbd',varargin);
if i
    varargin{i+1} = permute(varargin{i+1}(end:-1:1,end:-1:1,:),[2,1,3]);
end

gr.plotXSec(varargin{:},'dir','y');