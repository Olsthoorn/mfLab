%ANALYZEREGISGRIDS analyzes regis GRIDS
%
% ToDo: what is this ? (130429)
%
% TO 120401

%% Get afkortingen

dd='..\..\algemeen\';

files={'bestanden.csv';'eenheden.csv';'metadata.csv'};

delim=',';

%% remove ',' between " "

for i=1:length(files)
    fname=[dd files{i}];
    fp=fopen(fname,'r');
    fpout=fopen([fname '.out'],'w');
    while 1
        s=fgets(fp);
        if s==-1, break; end
        I=findstr('"',s);
        if ~isempty(I)
            J=findstr(delim,s);
            if ~isempty(J)
                for j=1:length(J)
                    for ii=1:2:length(I)
                        if J(j)>I(ii) && J(j)<I(ii+1)
                            s(J(j))=';';
                        end
                    end
                end
            end
        end
        fprintf(fpout,'%s',s);
    end
    fclose(fpout);
    fclose(fp);
end

%% scan the files
    
for i=1:length(files)
    fname=[dd files{i} '.out'];
    fp=fopen(fname,'r');
    
    [P,fname    ]=fileparts(fname);
    [P,fname,ext]=fileparts(fname);
    
    s=fgetl(fp); fgets(fp);
    n=length(findstr(delim,s))+1;
    M=textscan(fp,'%s','delimiter',',');
    M=M{1};
    M=reshape(M,n,length(M)/n)';
    fclose(fp);
    eval([fname '=M']);
end

%% getting layer tops
top=dir('*-t-*.mat');

layers.unit=cell(length(top),1); NUNIT=length(layers.unit);

for i=1:NUNIT
    layers.unit{i}=top(i).name(1:findstr('-t-',top(i).name)-1);
end

bot=dir('*-b-*.mat');
kh  =dir('*-kh-*.mat');
c   =dir('*-c-*.mat');
kv  =dir('*-kv-*.mat');

%% allocate space for 3D matrices all units
load(top(1).name);
[NROW,NCOL]=size(A);
layers.top=NaN(NROW,NCOL,NUNIT);

for i=1:NUNIT
    fprintf('%3d loading %s\n',i,top(i).name);
    load(top(i).name);
    layers.top(:,:,i)=A;
end
fprintf('...done\n');

%% sort the layers based on overlap
A=zeros(NROW,NCOL);
tic
for i=1:NUNIT
    for j=(i+1):NUNIT
        I=find(~isnan(layers.top(:,:,i)) & ~isnan(layers.top(:,:,j)));
        if ~isempty(I)
            II=NROW*NCOL*(i-1)+I;
            JJ=NROW*NCOL*(j-1)+I;
            if mean(layers.top(JJ))>mean(layers.top(II))
                fprintf('%d <--> %d\n',i,j)
                A(:,:)=layers.top(:,:,i);
                layers.top(:,:,i)=layers.top(:,:,j);
                layers.top(:,:,j)=A;
                u=layers.unit(i);
                layers.unit(i)=layers.unit(j);
                layers.unit(j)=u;
            end
        end
    end
end
toc
%%
layers.form=cell(NUNIT,1);
layers.long=cell(NUNIT,1);
for i=1:NUNIT
    idx=strmatchi(layers.unit(i),eenheden(:,1));
    layers.form(i)=eenheden(idx,3); layers.form{i}(layers.form{i}(:)=='"')=[];
    layers.long(i)=eenheden(idx,4); layers.long{i}(layers.long{i}(:)=='"')=[];
    I=find(~isnan(layers.top(:,:,i)));
    t=mean(layers.top(NROW*NCOL*(i-1)+I));
    fprintf('%-10s %-50s %-20s top avg =%12.2f\n',...
        layers.unit{i},layers.form{i},layers.long{i},t);
end

        