%READMFLOW reads an entire model into Matlab
%
% Example:
%    readModflow(xLL,yLL)
%
% This mfile is still under construction. It is simple to add new packages.
% For some of which readers are already made (see directory
% mflab\mfiles\read). The current one has been used to readin an iMod model
% prepared for the city of Delft in 2010.
%
% Any recognized packages will be read in. Recognation of these pacakges is
% done ussing the files extensions (alternative would be using nam file)
% Unknown packages will issue a warning.
%
% It may (will) be smarter to make use of the nam files. This is a simple
% adaptation that can be made. Just make this script a function with the
% name filename as argument and perhaps the xLL and  yLL to also compute
% the grid line coordinates correctly (assuming the model is unrotated).
% Currently or by default xLL and yLL are assumed zero.
%
% TO 120323

clear; close all

%% LL of GU's model
xLL=0; yLL=0;

%% Get files in local directory and their extensions as extra field
d=dir;
for i=length(d):-1:1
    if d(i).isdir,
        d(i)=[];
    else
        [~,~,ext]=fileparts(d(i).name);
        if length(ext)>1, d(i).ext=lower(ext(2:end)); end
    end
end

packList=cell(size(d)); % list of reconized packages (form file extension)

%% reading

% First get the dis file because the Ny,Nz,Nx are also needed by other
% packages
i=strmatchi('dis',{d.ext});
dis=readDIS(d(i).name);

% Then read the other packages
for i=1:length(d)
    [~,~,ext]=fileparts(d(i).name);
    fprintf('Filename = %s extension=%s\n',d(i).name,ext);
    switch d(i).ext
        case 'bas'
             bas.NLAY=dis.NLAY;
             bas.NROW=dis.NROW;
             bas.NCOL=dis.NCOL;
             bas=readBAS(d(i).name,bas);
             packList{i}='bas';
        case 'ba6'
             ba6.NLAY=dis.NLAY;
             ba6.NROW=dis.NROW;
             ba6.NCOL=dis.NCOL;
             ba6=readBA6(d(i).name,ba6);
             packList{i}='ba6';
        case 'bcf'
             bcf.NLAY=dis.NLAY;
             bcf.NROW=dis.NROW;
             bcf.NCOL=dis.NCOL;
             bcf=readBCF(d(i).name,bcf);
             packList{i}='bcf';
         case 'bc6'
             bc6.isTran=any(dis.isTran);
             bc6.NLAY=dis.NLAY;
             bc6.NROW=dis.NROW;
             bc6.NCOL=dis.NCOL;
             bc6=readBC6(d(i).name,bc6);
             packList{i}='bc6';
        case 'rch'
            rch.NROW=dis.NROW;
            rch.NCOL=dis.NCOL;
            rch.NPER=dis.NPER;
            rch=readRCH(d(i).name,rch);
            packList{i}='rch';
        case 'evt'
            evt.NROW=dis.NROW;
            evt.NCOL=dis.NCOL;
            evt.NPER=dis.NPER;
            evt=readEVT(d(i).name,evt);
            packList{i}='rch';
        case 'ghb'
            ghb.NPER=dis.NPER;
            ghb=readBCN(d(i).name,drn,'GHB');
            packList{i}='ghb';
        case 'drn'
            drn.NPER=dis.NPER;
            drn=readBCN(d(i).name,drn,'DRN');
            packList{i}='drn';
        case 'riv'
            riv.NPER=dis.NPER;
            riv=readBCN(d(i).name,riv,'RIV');
            packList{i}='riv';
        case 'chd'
            chd.NPER=dis.NPER;
            chd=readBCN(d(i).name,chd,'CHD');
            packList{i}='chd';
        case 'wel'
            wel.NPER=dis.NPER;
            wel=readBCN(d(i).name,wel,'WEL');
            packList{i}='wel';
        case 'hfb'
            hfb.NROW=dis.NROW;
            hfb.NCOL=dis.NCOL;
            hfb.NPER=dis.NPER;
            hfb=readHFB(d(i).name,hfb);
            packList{i}='hfb';
        case 'huf' % hydraulic units
            huf.NLAY=dis.NLAY;
            huf.NROW=dis.NROW;
            huf.NCOL=dis.NCOL;
            huf=readHUF(d(i).name,huf);
            packList{i}='huf';
        case 'adv' % advection
            adv.NLAY=dis.NLAY;
            adv.NROW=dis.NROW;
            adv.NCOL=dis.NCOL;
            adv=readADV(d(i).name,adv);
            packList{i}='adv';
        case 'dsp' % dispersion diffusion
            dsp.NLAY=dis.NLAY;
            dsp.NROW=dis.NROW;
            dsp.NCOL=dis.NCOL;
            dsp=readDSP(d(i).name,dsp);
            packList{i}='dsp';
        case 'btn' % basic transport
            btn.NLAY=dis.NLAY;
            btn.NROW=dis.NROW;
            btn.NCOL=dis.NCOL;
            btn.NPER=dis.NPER;
            btn=readBTN(d(i).name,btn);
            packList{i}='btn';
        case 'cfp' % conduit flow package
            cfp.NLAY=dis.NLAY;
            cfp.NROW=dis.NROW;
            cfp.NCOL=dis.NCOL;
            cfp=readCFP(d(i).name,cfp);
            packList{i}='cfp';
        case 'coc' % output control voor CFP
            coc.NLAY=dis.NLAY;
            coc.NROW=dis.NROW;
            coc.NCOL=dis.NCOL;
            coc.NPER=dis.NPER;
            coc=readCOC(d(i).name,coc);
            packList{i}='coc';
        case 'crch' % concentration of recharge
            crch.NLAY=dis.NLAY;
            crch.NROW=dis.NROW;
            crch.NCOL=dis.NCOL;
            crch.NPER=dis.NPER;
            crch=readCRCH(d(i).name,crch);
            packList{i}='crch';
        case 'lvda' % 
            lvda.NLAY=dis.NLAY;
            lvda.NROW=dis.NROW;
            lvda.NCOL=dis.NCOL;
            lvda.NPER=dis.NPER;
            lvda=readLVDA(d(i).name,lvda);
            packList{i}='lvda';
        otherwise
            fprintf('Don''t know package belonging to extension %s and file %s, ignored~\n',...
                d(i).ext,d(i).name);
    end
end

%% Finally save

% s is string with recognized packages
s=' dis'; for i=1:length(packList); s=[s ' ' packList{i}]; end

xGr=xLL+[0, cumsum(dis.DELR)];
yGr=yLL+[0; cumsum(dis.DELC)];

fprintf('Saving model.mat ...\n');
fprintf('save model.mat %s\n',s);

eval(['save model.mat xGr yGr' s]); % dynamic evaluation of packages present in model

fprintf('... model data saved in model.mat\n');
