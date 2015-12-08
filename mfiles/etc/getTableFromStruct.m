function str = getTableFromStruct(str,varargin)
%%GETTABLEFROMSTRUCT -- generates a table from a given struct as if it was
% read from an excel sheet. Usind in lineObj constructor (actually in
% pointObj.m) as an alternative to read in the data, desired if the date
% come directly from a shapefile.
%
% USAGE: str = getTableFromStruct(str[,prop,value,prop,value,...])
%
% str is a struct with fields, such as obtained by shaperead.
% prop,value,... are prop,value pairs to be added to the elements of the
% struct, i.e. enhancing the struct.

% Shapefile contains BoundingBox, which we don't want
if isfield(str,'BoundingBox')
    str = rmfield(str,'BoundingBox');
end

% add extra fields to struct read from varargin
% assert that varargin are prop,value pairs
% all str elements get the same fields and values this way
% if each str needs its own values, more code is needed here.
for i=1:2:numel(varargin)
    for iStr=numel(str):-1:1  
        if ischar(varargin{i})
            str(iStr).(varargin{i}) = varargin{i+1};
        end
    end
end
clear varargin

% remaining fields in struct
fldnames = fieldnames(str)';

% get class of data of each of the fields
for i=numel(fldnames):-1:1    
    fldclass{1,i} = class(str(1).(fldnames{i}));
end

% we want either the double and the char fiels, others are ignored
Idbl = strmatchi('double',fldclass);
Istr = strmatchi('char'  ,fldclass);

% double fieds are obligatory (should contain lat lon or X, Y
if isempty(Idbl)
    error('missing numeric data in struct');
end

% generate table (separately with numeric and textual data, exactly as if
% it came from getExcelData callind in the pointObj constructor.
for iStr = numel(str):-1:1
    
    % The headers for double and textual fields
    str(iStr).vertexHdr    = fldnames(Idbl);
    str(iStr).vertexTxtHdr = fldnames(Istr);
    
    % Determine length of table
    L=0;
    for i=Idbl,
        L = max(L,numel(str(iStr).(fldnames{i})));
    end

    %% Get the doubles;
    str(iStr).vertex    = NaN(L,numel(Idbl));
    str(iStr).vertexTxt =cell(L,numel(Istr)); 

    for i=numel(Idbl):-1:1
        values = str(iStr).(fldnames{Idbl(i)})(:);
        str(iStr).vertex(1:numel(values),i)   = values;
        str(iStr).vertex(numel(values):end,i) = values(end);
    end

    %% get the textual values
    if ~isempty(Istr)
        for i=numel(Istr):-1:1
            texts = str(iStr).(fldnames{Istr(i)});
            for iRow=1:L
                str(iStr).vertexTxt{iRow,i}   = texts;
            end
        end
    end    
    
end   

%% remove all fields except the vertex and vertexTxt fields
for i=1:numel(str(1).vertexHdr)
    str = rmfield(str,str(1).vertexHdr{i});
end
for i=1:numel(str(1).vertexTxtHdr)
    str = rmfield(str,str(i).vertexTxtHdr{i});
end





