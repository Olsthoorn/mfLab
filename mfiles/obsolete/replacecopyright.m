function replacecopyright(files)
% TO 091202

cpr={
''
'% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty'
'% under free software foundation GNU license version 3 or later'
};

d=dir(files);

for i=1:length(d)
    
    %% check if file already has copyright notice
    fid=Fopen(d(i).name,'r');
    
    s=fgets(fid);
    already=0;
    while isempty(s) || s(1)~=-1
        if ~isempty(strfind(s,'% Copyright 2009 Theo Olsthoorn'))
            if ~isempty(strfind(s, 'TU-Delft and Waternet'))
                already=1;
                break;
            end
        end
        s=fgets(fid);
    end
    fclose(fid);
    %% if not put the copyright notice in the file
    
    if ~already
        [p basename ext]=fileparts(d(i).name);

        % saveguard the old mfile
        movefile(  [basename ext],...
                   [basename ext '_old'])
        fid =Fopen([basename ext '_old'],'r');  % read from old mfile
        fido=Fopen([basename ext]       ,'w');  % read to new mfile

        s=fgetl(fid);
        found=0;  % to add copyright only once !!
        while isempty(s) || s(1)~=-1
            if~found && ~isempty(strfind(s,'% Copyright 2009 Theo Olsthoorn'))
                found=1;
                for ic=1:length(cpr)
                    fprintf(fido,'%s\n',cpr{ic});
                end
                  fgetl(fid); % skip old second copyright line
                s=fgetl(fid); % get next line
            else
                fprintf(fido,'%s\n',s);
                s=fgetl(fid); % continue copying
            end
        end
        fclose(fid);
        fclose(fido);

        fprintf('File %s done\n',d(i).name);
    else
        fprintf('File %s already has copyright notice ... skipped\n',d(i).name);
    end
end
    
function fid=Fopen(fname,mode)
% fid=Fopen(fname,mode)
% fopen with error control

fid=fopen(fname,mode);
if fid<0
    error('Can''t open file <<%s>> for mode <<%s>>',fname,mode);
end
