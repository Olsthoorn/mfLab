function [hdl]=mf_gcounter(arg1,strs1,vals,formats,strs2,varargin)
%MF_GCOUNTER puts a window with text and values on figure
%   useful to put a counter in a figure
%   [hdl]    = mf_gcounter(xyz,strings,values,formats,strings,'varargin)
%   xyz      = [x y] or [x y z] coordinates to put window
%   strs1    = {str1 str2 str2 ....} header strings
%   vals     = [val1 val2 val3 ....] counter values
%   formats  = {fmt1,fmt2,fmt3,...] format strings (may be just one for all)
%   varargin = {propery1,value,property2,value,....}
%
%EXAMPLE
%  hdl = mf_gcounter([x y z],{'speed1s','pressure','speed2'},[0.2 0.3 0.4],{'%8.1f'},...
%          {'m3/d','Pa','m/s'},'backgroundcolor',[0 1 0],'font','courier'}
%       mf_gcounter(hdl,{hdl,'speed1s','pressure','speed2'},[0.2 0.3 0.4],{'%.3g'},...
%          {'m3/d','Pa','m/s'},'backgroundcolor',[0 1 0],'font','courier'}
%
%  Second example, is to refresh. It takes hdl obtained in first call as its first argument
%
% TO 110422

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if  ~exist('strs2'  ,'var'), strs2   =''; end
if  ~exist('formats','var'), formats =''; end
if   ~exist('vals','var'),   vals    =[]; end

if ~iscell(strs1),   strs1  ={strs1};   end
if ~iscell(formats), formats={formats}; end
if ~iscell(strs2),   strs2  ={strs2};   end

M1=length(strs1); for i=1:length(strs1), M1=max(M1,length(strs1{i})); end
M2=length(strs2); for i=1:length(strs2), M2=max(M2,length(strs2{i})); end

s=cell(size(strs1));
for i=1:length(strs1)
    j=min(length(formats),i);
    k=min(length(strs2),i);
    if ~isempty(vals)
        s{i}=[strs1{i}, blanks(M1-length(strs1{i})), ' ',...
            sprintf(formats{j},vals(i)), ' ',...
            strs2{k}, blanks(M2-length(strs2{k}))];
    else
        s{i}=[strs1{i}, blanks(M1-length(strs1{i})),' ',...
            strs2{k}, blanks(M2-length(strs2{k}))];
    end
end
if nargout==0
    hdl=arg1;
    set(hdl,'string',s);
else
    x=arg1(1); y=arg1(2); if length(arg1)>2, z=arg1(3); else z=0; end
    hdl=text(x,y,z,s);
end

if exist('varargin','var'),
    for i=1:2:length(varargin)-1
        set(hdl,varargin{i},varargin{i+1});
    end
end
