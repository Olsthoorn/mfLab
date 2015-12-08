function table = getTable(basename,sheetnm,members,direction,use,varargin)
%GETTABLE reads a table into a struct with given regCols and the rest in field UserData
%
% USAGE:
%    table = getTable(basename,sheetnm,members,Vertical|Horizontal)
%
% A very powerful function to readout date from Microsoft Excel Spreadsheets
% It is used in many core function of mfLab to robustly and generically
% extract data from Excel spreadsheets, such as when construcging well
% objects.
%
% The entire table (headers with values below) is red. Any columns that
% correspond to the strinngs (fields) specified in members become fields in
% table. Any other fields in the spreadsheet that are not in members,
% become fields of table.UserData.
% This allows reading in all data from an excel spreadsheet table,
% inlcuding those columns that do not match fixed fields of objects such as
% wellObjects. These fields become fields of table.UserData and, therefore,
% do not get lost. This allows to extend the fields in objects without
% redefining the objects.
%
% varargin can contain 'fill',true|false or 'fill' or 'nofill'
% triggers filling empty spreadsheet cells or setting NaNs for them.
%
% TO 121231

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

[fill, ~ ] = getProp(varargin,'fill',true);

[numHdr,num,txtHdr,txt]=getExcelData(basename,sheetnm,direction,'fill',fill);


%% 3rd argument specifies which rows of the table will be used
% use may be a list of numbers to be find in num column 'nr'
% or it may be a reg expresion to be checked on txt 'name'.

if nargin>5 && ~isempty(use)
    if isnumeric(use)
        iNr= strmatchi('nr',numHdr);
        remove = ~ismember(num(:,iNr),use);
        num(remove,:)=[];
        txt(remove,:)=[];
    elseif ischar(use) % use is considered a reg expression
        expr = use;
        hit = false(size(num(:,1)));

        iName = strmatchi('name',txtHdr);
        for itest=1:numel(hit)
            strt = regexp(txt(itest,iName),expr,'start');
            if ~isempty(strt{:})
                 hit(itest)=true;
            end
        end
        remove = ~hit;
        num(remove,:)=[];
        txt(remove,:)=[];
    end
    
    if isempty(num)
         error('mfLab:getTable:nohits',...
             ['%s: Third argument in call does not provide any hits. Check input,\n',...
             'it must be a list matching column ''Nr'' or a regular expression\n',...
             'that is applicable to colmn ''names''.\n'],...
             mfilename);
    end
end


%% without members, use all
if nargin<3, members = [numHdr txtHdr]; end

try % must be nr in wellObj etc
    numHdr{strmatchi('Nr',numHdr,'exact')}='nr';
catch %#ok
    % ignore
end

nRec = size(num,1);

for i= numel(numHdr):-1:1
    im = strmatchi(numHdr{i},members,'exact');
    if im        
        for j = nRec:-1:1
            table(j).(members{im}) = num(j,strmatchi(numHdr{i},numHdr,'exact'))'; % also gets double hdr like x x y y z z
        end
    else
        for j = nRec:-1:1
            table(j).UserData.(numHdr{i}) = num(j,i);
        end
    end
end

for i = numel(txtHdr):-1:1
    im = strmatchi(txtHdr{i},members);
    if im
        for j = nRec:-1:1
            table(j).(members{im}) = txt{j,i};
        end
    else
        for j = nRec:-1:1
            table(j).UserData.(txtHdr{i}) = txt{j,i};
        end
    end
end

nu= now;

for  j = nRec:-1:1
    table(j).created = nu;
end
