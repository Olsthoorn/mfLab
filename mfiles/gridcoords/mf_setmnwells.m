function [MNW,PNTSRC,WEL]=mf_setmnwells(basename,xGr,yGr,Z,QCstr,HK)
%MF_SETMNWELLS puts wells in the model grid.
%
%  These wells are speicfied in a worksheet "wells" in the workbook basename.
%  The function works for both multinode wells and ordinary wells, i.e. where
%  wells are defined in inividual cells.
%
% Example:
%    [MNW, PNTSRC, WEL]=mf_setmnwells(basename,xGr,yGr,Z,dims,QCstr, [,HK])
%
%    MNW is a struct containing mnwObj objects. An mnwobject contains the data
%    necessary to construct the multinode well II input (Konikow, Hornberger,
%    Halford, Hanson 2009, ISBN 978-1-4113-2488-6)
%
%    PNSRC is the point source input as required by MT3DMS and SEAWAT to specify
%    the concentration of the wells upon injeciton into the model.
%
%    WEL is the data required by teh WEl package.
%
%    basename, xGr and yGr as usual. Z is the 3D array specifying top and
%    bottom of layers. If size(Z) == [1 1 Nz] all layers are assumed
%    uniform.
%    dims is the size of the model grid (cells), so dims=size(IBOUND)
%    QCstr are strings of the first uniform letters of the headings of the
%    well flows and concentrations in the PER worksheet. For instance
%    '{Q_,C_,T_}' would corrsespond to columns Q_1, Q_2, Q_3 etc for the
%    flow columns, C_1, C_2 C_3 etc for the concentration columns of the
%    first component (species) and T_ with T_1, T_2, T_3 etc concentration
%    columns of the second species, as for instance temperature. And so on
%    if more species are to be included in a MT3DMS or Seawat simulation
%
% Used in:
%    mflab/examples/Geohydrology2/CoertStrikker
%    mflab/examples/Geohydrology2/Wolfs-DSI-test
%    mflab/examples/mf2k/MNW1_Reilly
%    mflab/examples/swt_v4/ATES-WKO/AXIAL-MONO-Amsterdam
%
% ToDo: replace by MNW1Obj
%
% SEE ALSO: xsConf gridObj
%
%   TO 110426 110806

[welnams,welvals,weltxthdr,weltxt]=getExcelData(basename,'MNW','hor'); % get well info
[pernams,pervals]        =getPeriods(basename);                 % get stress periods info
[laynams,layvals]        =getLayers(basename,gr);

LAYCBD=layvals(:,strmatchi('LAYCBD',laynams));

welvals=welvals(welvals(:,strmatchi('Nr',welnams))>0,:); % only use lines with positive well nrs

% fill the mnwObj objects one by one
NPER=size(pervals,1);

if ~iscell(QCstr), QCstr={QCstr}; end  % make sure cells are used for uniform treatment

NWEL=size(welvals,1);  % number of active wells
MNW=mnwObj(NWEL);      % generate array with NWEL empty mnwObj objects

for iw=1:NWEL
    wellid   = weltxt{iw,strmatchi('WELLID',weltxthdr')};
    welltype = weltxt{iw,strmatchi('TYPE'  ,weltxthdr)};
    
    wellnr = welvals(iw,strmatchi('Nr',welnams));
    
    qstr=sprintf('%s%d',QCstr{1},wellnr); % find stress period Q for this well
    I = strmatchi(qstr,pernams,'exact');
    if I==0, error('mf_setmnwells:qstr',...
            'mf_setmnwells: Can''t find columns %s in worksheet PER !',qstr);
    end
    Q = pervals(:,I);
    
    if length(QCstr)>1     % then continue with concentrations for all components
        NCOMP=length(QCstr)-1;   % the rest of the QCstr indicate component headings
        C=NaN(NPER,NCOMP);       % reserve space
        for ic=1:NCOMP
            cstr=sprintf('%s%d',QCstr{ic+1},wellnr);    % conc for this component of this well
            I = strmatchi(cstr,pernams,'exact');        % search string for headers of per sheet
            if I==0, error('mf_setmnwells:cstr',...
                    'mf_setmnwells: Can''t find column %s in worksheet PER !',cstr);
            end
            C(:,ic)=pervals(:,I); % search string for headers of per sheet
        end
    else
        C=[];
    end
    
    % generate MNW obj for this well
    MNW(iw)=mnwObj(wellid,welltype,...
        xGr,yGr,Z,LAYCBD,...
        welnams,welvals(iw,:),...
        Q,C,HK);

end

% Because the MNW objects know whether they are WEL or MNW
% the can correctly output the PNTSRC data one after the other
if nargout>1
    PNTSRC=[];
    for i=1:length(MNW)
        if i==1
            PNTSRC=MNW(i).pntsrc;
        else
            PNTSRC=[PNTSRC; MNW(i).pntsrc];
        end
    end
    PNTSRC=sortrows(PNTSRC);
end

if nargout>2
    WEL=[];
   
    I=find([MNW.isaWEL]);
    for i=I
        if i==I(1),
            WEL=MNW(i).wel();
        else
            WEL=[WEL; MNW(i).wel()];
        end
    end
    WEL=sortrows(WEL);
end
        