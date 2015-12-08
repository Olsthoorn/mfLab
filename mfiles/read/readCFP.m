function cfp=readCFP(fname,pth,cfp)
%READCFP reads MDOFLOW-CFP input file (Conduit flow package)
%
% Example:
%    readCFP([pth fname] ,cfp);
%
% Todo: check and possibly fix it when CFP is fully implemented and tested
%
% TO 090708 090713

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0
fprintf('# MATLAB writeCFP %s\n',datestr(now));
fid=fopen([pth fname],'r');

%1
fprintf(fgets(fid)); % # mode
cfp.mode=fscanf(fid,'%d',1); fprintf(fgets(fid));

if cfp.mode==1 || cfp.mode==3
    
    %2-4
    fprintf(fgets(fid));  % # Data for mode 1
    fprintf(fgets(fid));  % # of nodes, pipes and layers
    cfp.NNODES =fscanf(fid,'%d',1);
    cfp.NPIPES =fscanf(fid,'%d',1);
    cfp.NLAYERS=fscanf(fid,'%d',1);
    fgets(fid);
    %5-6
    fprintf(fgets(fid));  % mean groundwater temperature in pipes
    cfp.Temp=fscanf(fid,'%f',1);
    fgets(fid);

    %7-8 Connections: Node numbers, 6 neighbors and connected 6 pipes
    fprintf(fgets(fid)); % #No  mc mr ml Nb1  Nb2  Nb3  Nb4 Nb5  Nb6 tb1 tb2 tb3 tb4 tb5 tb6
    C=textscan(fid,'%d   %d %d %d  %d %d %d %d %d %d  %d %d %d %d %d %d',cfp.NNODES,'CollectOutput',1);
    cfp.Connections=C{1};
    fgets(fid);

    %% 9-12 GEOMODE
    fprintf(fgets(fid)); % node elevations, two possibilities
    fprintf(fgets(fid)); % 1 node #, elevation above datum (1 line for each node)
    fprintf(fgets(fid)); % 2 number of nodes, distance form vertical centroid (only one line)

    % GEOMODE, see which mode appies:
    p=ftell(fid);
    if fscanf(fid,'%d',1)==cfp.NNODES, cfp.GEOMODE=2; else cfp.GEOMODE=1; end;
    fseek(fid,p,'bof');

    % GEOHEIGHT
    if cfp.GEOMODE==1, % one line per node
        C=textscan(fid,'%f %f',cfp.NNODES,'CollectOutput',1);
        cfp.GEOHEIGHT=C{1};
    else % only one line with distance form cell centroid
        C=textscanf(fid,'%d %f',1,'CollectOutput',1);
        cfp.GEOHEIGHT=C{1};
    end
    fgets(fid);

    %% 13-14   #surface dependent exchange (set 1) or constant exchange (set 0)'
    fprintf(fgets(fid));  % SA_EXCHANGE, CFP-pipe conductance (set 1) or user computes pipe conductance (set 0)
    cfp.SA_EXCHANGE=fscanf(fid,'%d',1); fgets(fid);

    %% 15-22 (Mewton Raphson parameters)
    fprintf(fgets(fid)); % #criterion for convergence,
    cfp.EPSILON=fscanf(fid,'%f',1); fgets(fid);

    fprintf(fgets(fid)); % #maximum number for loop iterations,
    cfp.NITER=fscanf(fid,'%d',1);  fgets(fid);

    fprintf(fgets(fid));  % #parameter of relaxation,
    cfp.RELAX=fscanf(fid,'%f',1); fgets(fid);

    fprintf(fgets(fid));  % #newton raphson print flag,
    cfp.P_NR=fscanf(fid,'%d',1); fgets(fid);

    %% 23-25 Tube data
    fprintf(fgets(fid)); % #data for tube parameters
    fprintf(fgets(fid)); % #no. diameter  tortuosity  roughness   lreynolds treynolds

    C=textscan(fid,'%f %f %f %f %f %f',cfp.NPIPES,'CollectOutput',1);
    cfp.PIPES=C{1};
    fgets(fid);

    %% 26-27 Node data
    fprintf(fgets(fid)); % #node heads (if head unequal -1 the head is fixed
    C=textscan(fid,'%d %d',cfp.NNODES,'CollectOutput',1);
    cfp.HD_OR_FLG=C{1};
    fgets(fid);

    %% 28-29
    fprintf(fgets(fid)); % #exchange terms for flow between continuum and pipe-network
    C=textscan(fid,'%f %f',cfp.NNODES);
    cfp.K_EXCHANGE=C{1};
    fgets(fid);

    if cfp.mode==1, return; end

end

%% 30-31 (only if mode>1)
fprintf(fgets(fid));
fprintf(fgets(fid));
fprintf('Conduit layers must be convertible, LAYCON=3 in BCF or LAYTYPE>0 in LPF or LTHUF>0 in HUF Package\n');

% 32
s=fgets(fid); fprintf(s);
cfp.NCL=sscanf(s,'%d',1);  % total number of conduit layers

%33 34
fprintf(fgets(fid));
s=fgets(fid); fprintf(s);
cfp.CL=sscanf(s,'%d',[1,cfp.NCL]);  % modflow layers that are conduit layers

%35-36
fprintf(fgets(fid));
s=fgets(fid); fprintf(s);
cfp.LTMP=sscanf(s,'%f',1);  % layer temperature

%37-39
fprintf(fgets(fid));
for i=1:cfp.NCL
    fprintf(fgets(fid));
    s=fgets(fid); fprintf(s);
    C=sscanf(s,'%f %f %f',[1,3]);
    cfp.VOID(i)=C(1);
    cfp.LCRITREY_L(i)=C(2);
    cfp.TCRITREY_L(i)=C(3);
end
fclose(fid);
