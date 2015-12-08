function [well,WEL,PNTSRC,NPER]=mf_setwells(basename,gr,HK,well)
%MF_SETWELLS puts wells in the grid when they are specifie in the sheet with 
%
% Examples:
%  [well,WEL,PNTSRC,NPER]=mf_setwells(basename,gr,HK,wellsheetname)
%  [well,WEL,PNTSRC,NPER]=mf_setwells(basename,gr,HK,well)
%  [well,WEL,PNTSRC,NPER]=mf_setwells(basename,xGr,yGr,zGr,HK,wellsheetname | well [,welnrs [,LAYCBD]]))
%
%    well is a well object containing all relevant information about the well and
%    its location in the grid.
%    WEL    is the input array for the well package.
%    PNTSRC is the input array for the pointsourcen when running MT3D or Seawat
%    To generate these arrays, it is assume that well flows are specified
%    The well object provides these arrays for each well and all stress
%    periods for as far as well(iw).Q(iPer)~=0.
%    if wells are obtained from another source as wellObj, you may use well
%    instead of wellsheetnm.
%
%    in the PER worksheet as one column per well with the total extraction
%    per stress period. The headers of the well injection
%    columns must be Q_1 Q_2 Q_3 etc. There must be no other columns starting
%    with Q_.
%
%    If PNTSRC is required, the concentratons for each well must be specified
%    in columns of the PER sheet named C1_1 C1_2 C1_3 etc. for concentrations
%    of component 1, C2_1 C2_2 C2_3 for component two etc.
%
%    There must be no other names starting with C1_ C2_ etc in the PER worksheet.
%    HK is required to distribute the well extraction over more than one
%    layer on the basis of the transmissivity.
%
% Used in:
%    mflab/examples/Geohydrology2/Joost
%    mflab/examples/Geohydrology2/RemkoNijzink
%    mflab/examples/swt_v4/BanglaDesh
%
% ToDo: replace by wellObj
%
% See also: wellObj
%
%   TO 110426 120103 120408


if ischar(well)
    
    wellsheetname=well; clear well;

    [welnams,welvals,weltxthdr,weltxt]=getExcelData(basename,wellsheetname,'hor');

    welvals=welvals(welvals(:,1)>0,:); % so we can easily switch wells on and off
     % and even have arbitrary lines in the worksheet welsheetname as long as Nr<0 or Nr=[];

    %% Check existance of header Nr
    if ~strmatchi('Nr',welnams,'exact')
        error('mf_setwells:wellsheetname:NoWellNrs',...
            'Header Nr missing in sheet <<%s>>!\n',wellsheetname);
    end

    NrCol=1;

    %% Remove unwanted wells and possibly empty lines
    welvals=welvals(welvals(:,NrCol)>0,:);

    if isempty(welvals),
        error('mfLab:mf_setwells:noActiveWells',...
              'No active wells in sheet %s\n',wellsheetname);
    end

    %% Create well as wellObj instantiations

    well(size(welvals,1),1)=wellObj(); % memory allocation

    for iw=size(welvals,1):-1:1 % running backward to allow removing while checking

        if ~isempty(weltxthdr)
            j=strmatchi('Name,',weltxthdr');
            if j~=0, well(iw).name=weltxt{iw,j}; end
            j=strmatchi('Remark',weltxthdr');
            if j~=0, well(iw).remark=weltxt{iw,j}; end
        end

        %% see of rw is given
        ii=strmatchi('diam',welnams);
        if ii>0,     rw = welvals(iw,ii)/2;
        else         ii=strmatchi('rw',  welnams);
            if ii>0, rw = welvals(iw,ii);
            else     rw=[];
            end
        end

        well(iw) = wellObj(basename,'wells');welvals(iw,strmatchi('Nr',welnams)),...
                            welvals(iw,strmatchi('x' ,welnams)),...
                            welvals(iw,strmatchi('y' ,welnams)),...
                            welvals(iw,strmatchi('z1' ,welnams)),...
                            welvals(iw,strmatchi('z2' ,welnams)),...
                            rw);

    end
else
    if ~strcmpi(class(well),'wellObj')
        error('mf_setwells:input:wellInput',...
            '%s: input well|wellsheetname  must be either char or wellObj',mfilename);
    end
end

%% put wells in grid
if exist('HK','var') 
    well=well.setWell(gr,HK);        
end

if nargout<2, return; end   % no stress period data required
    
%% Preambule

[pernams,pervals]=getPeriods(basename);
NPER=size(pervals,1);

%%%%% THE COLUMNS IN THE PER SHEET REFERRING TO THE WELLS MUST BE NAMED
%%%%% Q_%d where %d is the well number.

QCOL=NaN(size(well));

for iw=1:length(well)
    qstr=sprintf('Q_%d',well(iw).nr);  % Q columns in PER will be Q_1, Q_2 Q_3 etc.
    QCOL(iw)=strmatchi(qstr,pernams,'exact');
    if QCOL(iw)==0
        error('mf_setwells:QCol:PER',...
            'No Q column <<%s>>in sheet PER for well Nr %d',qstr,well(iw).Nr);
    end
end

%%
Dt=pervals(:,strmatchi('PERLEN',pernams));

for iw=1:length(well)
    well(iw).Q  = pervals(:,QCOL(iw))';
    if any(isnan(well(iw).Q))
        error('mf_setwells:pervals:QhasNaN',...
            '%s.Q has NaN(s), check worksheet PER column <<Q_%d>> for missing Q data!\n',well(iw).name,well(iw).nr);
    end
    well(iw).Dt = Dt';
    well(iw).t  = cumsum(well(iw).Dt,2); % may be overruled by setWelCout(well,C,iComp), using [C.time]
end

if nargout<3,
    WEL=well.WEL;
    return;
end % return, no PNTSRC needed

 
%% if PNTSRC is required (MT3DMS or SEAWAT
    
[MT3nams,MT3vals]=getExcelData(basename,'MT3D','vertical');
NCOMP=MT3vals(strmatchi('NCOMP',MT3nams),1);

%%%%% THE COLUMNS IN THE PER SHEET REFERRING TO CONCENTRATIONS MUST BE
%%%%% NAMED AS FOLLOWS  Cm_n where the first m is the well number and n
%%%%% the species number: i.e.
%%%%% C1_1 C2_1 C3_1 C1_2 C2_2 C3_2 for three species and two wells

% Find the columns in worksheet PER holding the concentrations for the well
for iw=1:length(well)
    well(iw).NCOMP=NCOMP;
    well(iw).C = NaN(NCOMP,NPER);
    
%     well(iw).Cout=NaN(size(well(iw).C));
 
    for ic=1:NCOMP
        cstr=sprintf('C%d_%d',well(iw).nr,ic);
        CCOL=strmatchi(cstr,pernams,'exact');
        if CCOL==0
            error('mf_setwells:CCOL:notfound',...
                ['Can''t find well-concentration column header <<%s>>for component %d well %d in PER worksheet,\n',...
                'Make sure your concentration column headings are "Ca_b" where a=species number and b=well number.\n',...
                'Well number corresponds to the number in the sheet where you specified your wells including\n',...
                'their well numbers. So the well number is specified by you independently of mfLab.\n'],...
                cstr,ic,well(iw).nr);
        end
        well(iw).C(ic,:) =pervals(:,CCOL)'; % plug C into well object
    end
end

[WEL,PNTSRC]=well.PNTSRC();

