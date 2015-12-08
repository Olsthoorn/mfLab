%% Generate PEST model input template file
%
% How to generate template, insturciton and control file for PEST ??
% Use this genpestfiles after you adapted it to your model.
%
% PEST has its own parameter values, that can be normally used with their
% default values. I put thise in the accompanying workbook in the worksheet
% PEST. You would normally not need to change these. Note that the
% parameter names in this workhsheet have a comment indicating the meaning
% and often useful values to set it at. Bur for more details look in the
% PEST manual.
%
% However, your parameters, paramter groups and observations and prior
% information rules are deicated to the specific model at hand. Therefore,
% you have to adapt them for each new model. Just replace the cell matrices
% below with the ones that suir your model. Yous their format to make sure
% it works.
%
% TO 100728

fprintf('# generating pestfiles, %s #\n',datestr(now)); 

load name  % loads basename

% The more or less fixe (default) parameters are in the worksheet PEST

%% Definition of the parameter groups
% PEST input variables that define how derivatives are calculated pertain
% to parameter groups rather than individual parameters.
% 1 PARGPNME  name of parameter group (case insensitive <= 8 chars)
% 2 INCTYP    relative|rel_to_max|absolute way of incrementing variable 
% 3 DERINC    increment (interpreted according to INCTYP)
% 4 DERINCLB  absolute lower bound of parameter increments
% 5 FORCEN    always_2|always_3|switch nr of points for derivative
% 6 DERINCMUL multipy with DERINC incase of 3 point derivative calculation
% 7 DERMTHD   parabolic|bestfit|outside_pts

% Adapt to your model
PARGP={ %PARGPNME INCTYP DERINC DERINCLB FORCEN DERINCMUL DERMTHD
         'pgp1','relative', 0.01, 0.00001, 'switch', 2.0, 'parabolic';
         'pgp2','relative', 0.01, 0.00001, 'switch', 2.0, 'parabolic';
         };

%% Definition of the parameters for the calibration
% Each parameter neither tied nor fixed must belong to a parameter group
% i.e. one of the groups defined in PARGP below. Tied or fixed parameters
% may belong to a dummy group none, but this is immaterial.
%
% 1 PARNAME   parameter name, case insensitive <=8 chars
% 2 PARTRANS  none|log|fixe|tied  parameter tranformation
% 3 PARCHGLIM relative|factor     parameter relative or factor limited
% 4 PARVAL1   initial parameter value (value of first iteration)
% 5 PARLBND   parameter lower bound
% 6 PARUBND   parameter upper bound
% 7 PARGP     groupname|none  parameter group to which parameter belongs
% 8 SCALE     value by which parameter if multiplied before offset is applied
% 9 OFFSET    offset applied to parameter after scaling. PEST sees bp=(bm-o)/s
% Note: each parameter must be cited in the template file

% Adapt for each model
PAR={ % PARNAME PARTRANS  CHGLIM PARVAL1 PARLBND PARUBND PARGP SCALE OFFSET
    'phk1'  ,'log'  , 'factor', 1.5, 0.2, 10, 'pgp1', 1.0, 0.0;
    'phk2'  ,'log'  , 'factor', 0.8, 0.2, 10, 'pgp1', 1.0, 0.0;
    'phk3'  ,'fixed', 'factor', 1  , 0.2, 10, 'pgp1', 1.0, 0.0;
    'pvkcb1','log'  , 'factor', 2.0, 0.2, 10, 'pgp1', 1.0, 0.0;
    'pvkcb2','fixed', 'factor', 1,   0.2, 10, 'pgp1', 1.0, 0.0;
    'pdrn'  ,'log'  , 'factor', 0.5, 0.2, 10, 'pgp1', 1.0, 0.0;
    'prch1' ,'fixed', 'factor', 1  , 0.2, 10, 'pgp2', 1.0, 0.0;
    'prch2' ,'fixed', 'factor', 1  , 0.2, 10, 'pgp2', 1.0, 0.0};

%% Definition of parameters second part

% TIEDPAR  name of tied parameter
% TIED2PAR name of parameter tied to (must itself be niether tied nor fixed)
% use as many lines as there are tied parameters or { }; if none

PARTIED={ % TIEDPAR TIEDTOPAR
    };

%% Names of observation groups

% List of observation groupps, each name case insensitive and <=8 chars.

% Adapt to your model
OBGNME={
    'ogrp1';
    'ogrp2';
    'ogrp3';
    'ogrp4'};

%% Definition of the observation data

% 1 OBSNAME unique name of observations, case insensitive, <= 8 chars
% 2 OBSVAL  fied or lab measurement
% 3 WEIGT   weight attached to this observation >=0
% 4 OBGNME  name of observation group to which observation pertains

% Adapt to your model
OBS={ % NAME VAL WEIGTH OBSGPNAME
    'Obs1',  124.15, 1.0, 'ogrp1';
    'Obs2',  151.43, 1.0, 'ogrp1';
    'Obs3',  222.72, 1.0, 'ogrp1';
    'Obs4',  258.77, 1.0, 'ogrp1';
    'Obs5',  222.77, 1.0, 'ogrp1';
    'Obs6',  151.43, 1.0, 'ogrp1';
    'Obs7',  234.43, 1.0, 'ogrp1';
    'Obs8',  258.40, 1.0, 'ogrp1';
    };

%% Model command line, the line PEST excutes to run the model
%
% Assuming PEST is run from the local directory (DOS path ok or full path
% included in command such as in the batfile runpest.bat, which holds
%REM running pest
%"Z:\tolsthoorn On My Mac\GRWMODELS\PEST\PEST12PC\PEST.exe" ex1pest.pst

% assuming run.m exists locally
% run.m should set the paths for mfLab and start mf_setup
% the logfile contains what matlab normally writes to its command window

%PESTCMD='matlab -nosplash -nodesktop -r run.m -logfile run.log';
PESTCMD='run';
%
%alternatively, with this line in run.bat
%PESTCMD=run

%% Template file names paired with input file names produced from templates

% in the mfLab context you generally need only one template file and one
% instruction file. Therefore, you generally don't have to adapt the
% following line. But note that all of thise files have basename equal to
% that of the model (basename) with different extensions

% template and model input file read by mf_adapt.m
TEMPFLE={[basename,'.tpl']}; INFLE ={[basename,'.inp']};
% instruction file and model output file read by PEST
INSFLE ={[basename,'.pif']}; OUTFLE={[basename,'.out']};

%% Definition of prior information / rules

% PROTOCOL
% No more than 300 chars per line, use & + space to continue no next line
% Each item separte by at least one space from its neighbors
% Start with label <8 chars, case insensitive like all other char
% variables.
%
% To the left of ths "=" sign: one or more combinations of
% PIFAC * PARNAME|log(PARNAME) + ...
% All PARNAME use must be adjustable parameters.
% Each parameter can be referenced only once in PRIOR equation.
% The parameter factor must be supplied.
% log(PARNAME) is necessary if PARNAME is log transformed.
% logbase is 10.
%
% To the right side of the "=" sign are two real variables PIVAL and WEIGHT.
% PIVAL is the value that the left side of the the equation aimes to achieve
% WEIGHT is the weight of this prior information rule/article.
% WEIGHT should preferably be inversly proportional to the standard
% deviation of PIVAL, but must always be >=0
%
% No two prior information articles must say the same thing.
%

% Adapt to your model
%PRIOR={ % PILBL PIFAC * PARNME + PIFAC * log(PARNME) .... = PIVAL WIEGHT
%    };
PRIOR={};

%% What follows should not have to be adapted for any model.
% It generates the tamplate file, the initial infile, instruction file
% and the PEST control file for you

%% Generate template file for the parameters

fid=fopen(TEMPFLE{1},'wt');
fprintf(fid,'ptf #\n');
for i=1:size(PAR,1),
    fprintf(fid,'#%-15s#\n',PAR{i,1});
end
fclose(fid);

%% Generate initial model input file

fid=fopen(INFLE{1},'wt');
for i=1:size(PAR,1)
    fprintf(fid,'%15g\n',PAR{i,4});
end
fclose(fid);

%% Generate PEST instruction file

fid=fopen(INSFLE{1},'wt');
fprintf(fid,'pif %%\n');
for i=1:length(OBS)
    fprintf(fid,'l1 !%s!\n',OBS{i,1});
end
fclose(fid);

%% The model output file with the computed observatins is generated by mf_analyze.mm

%% Pest parameter fle basename.PAR

% this file is generated by PEST. Pest uses it to generate the
% model intput files. So the user has nothing to deal with the
% pest par file. Nevertheless the contents in the mfLab context
% is the same as the INFILE obtained from the TMPLFLE, which is
% the one use by mfLab


%% Counters
NPAR    =size(PAR,1);
NOBS    =size(OBS,1);
NPARGP  =size(PARGP,1);
NPRIOR  =size(PRIOR,1);
NOBSGP  =size(OBGNME,1);
NTPLFILE=size(TEMPFLE,1);
NINSFLE =size(INSFLE,1);

%% Generate PEST control file
PESTFILE=[basename,'.pst'];

[parnams,parvals]=getExcelData(basename,'PEST','Vertical');

fid=fopen(PESTFILE,'wt'); 
fprintf(fid,'pcf\n');
fprintf(fid,'* control data\n');
fprintf(fid,'RESTART estimation\n');
fprintf(fid,'%d %d %d %d %d\n',NPAR,NOBS,NPARGP,NPRIOR,NOBSGP);

PRECIS=parvals(strmatchi('PRECIS',   parnams),1);
if PRECIS==1, PRECIS='single'; else PRECIS='double'; end
DPOINT=parvals(strmatchi('DPOINT',   parnams),1);
if DPOINT==0, DPOINT='nopoint'; else DPOINT='point'; end
fprintf(fid,'%d %d %s %s\n',NTPLFILE,NINSFLE,PRECIS,DPOINT);   

fprintf(fid,'%g %g %g %g %d %d LAMFORGIVE\n',...
    parvals(strmatchi('RLAMBDA1', parnams),1),...
    parvals(strmatchi('RLAMFAC',  parnams),1),...
    parvals(strmatchi('PHIRATSUF',parnams),1),...
    parvals(strmatchi('PHIREDLAM',parnams),1),...
    parvals(strmatchi('NUMLAM',   parnams),1),...
    parvals(strmatchi('JACUPDATE',parnams),1)); % JACUPDATE is new version 12
fprintf(fid,'%g %g %g\n',...
    parvals(strmatchi('RELPARMAX',parnams),1),...
    parvals(strmatchi('FACPARMAX',parnams),1),...
    parvals(strmatchi('FACORIG',  parnams),1));
fprintf(fid,'%g\n',...
    parvals(strmatchi('PHIREDSHW',parnams),1));
fprintf(fid,'%d %g %d %d %g %d\n',...
    parvals(strmatchi('NOPTMAX',  parnams),1),...
    parvals(strmatchi('PHIREDSTP',parnams),1),...
    parvals(strmatchi('NPHISTP',  parnams),1),...
    parvals(strmatchi('NPHINORED',parnams),1),...
    parvals(strmatchi('RELPARSTP',parnams),1),...
    parvals(strmatchi('NRELPAR',  parnams),1));
fprintf(fid,'%d %d %d\n',...
    parvals(strmatchi('ICOV',parnams),1),...
    parvals(strmatchi('ICOR',parnams),1),...
    parvals(strmatchi('IEIG',parnams),1));

%% Parameter groups
fprintf(fid,'* parameter groups\n');
for i=1:size(PARGP,1)
   fprintf(fid,'%s %s %g %g %s %g %s\n',...
       PARGP{i,1},PARGP{i,2},PARGP{i,3},PARGP{i,4},PARGP{i,5},...
       PARGP{i,6},PARGP{i,7});
end

%% Parameter data
fprintf(fid,'* parameter data\n');
for i=1:size(PAR,1)
   fprintf(fid,'%s %s %s %g %g %g %s %g %g\n',...
       PAR{i,1},PAR{i,2},PAR{i,3},PAR{i,4},PAR{i,5},...
       PAR{i,6},PAR{i,7},PAR{i,8},PAR{i,9});
end
for i=1:size(PAR,1)
   if strcmpi(PAR{i,2},'TIED')
    fprintf(fid,'%s %s\n',PAR{i,1},PAR{i,10});
   end
end

%% Observation groups
fprintf(fid,'* observation groups\n');
for i=1:size(OBGNME,1)
   fprintf(fid,'%s\n',OBGNME{i});
end

%% Observation data
fprintf(fid,'* observation data\n');
for i=1:size(OBS,1)
   fprintf(fid,'%s %g %g %s\n',...
       OBS{i,1},OBS{i,2},OBS{i,3},OBS{i,4});
end

%% Model command line

% the command that pest executes

fprintf(fid,'* model command line\n');
fprintf(fid,'%s\n',PESTCMD);

%% Model input/output
fprintf(fid,'* model input/output\n');
for i=1:size(TEMPFLE,1)
   fprintf(fid,'%s %s\n',TEMPFLE{i},INFLE{i});
end
for i=1:size(INSFLE,1)
   fprintf(fid,'%s %s\n',INSFLE{i},OUTFLE{i});
end

%% PRIOR rules
fprintf(fid,'* prior information\n');
for i=1:size(PRIOR,1)
   % PILBL PIFAC * PARNAME + PIFAC * log(PARNAME) + ... = PIVAL WEIGHT'\n');
   fprintf(fid,'%s',PRIOR{i,1});  % PRIOR LABEL (no spaces in front
   for j=2:size(PRIOR,2) % Rest of prior info
       if ~isempty(PRIOR(i,j))
           if isnum( PRIOR{i,j}), fprintf(fid,' %g',PRIOR{i,j}); end
           if ischar(PRIOR{i,j}), fprintf(fid,' %s',PRIOR{i,j}); end
       end
   end
   fprintf(fid,'\n');
end

fclose(fid);

fprintf('Ready:\n');
fprintf('pestfiles generated:\n');
fprintf(' %s\n',PESTFILE);
for i=1:length(TEMPFLE), fprintf(' %s',TEMPFLE{i}); end; fprintf('\n');
for i=1:length(INFLE)  , fprintf(' %s',INFLE{i});   end; fprintf('\n');
for i=1:length(INSFLE) , fprintf(' %s',INSFLE{i});  end; fprintf('\n');

