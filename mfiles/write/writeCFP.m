function writeCFP(basename,cfp)
%WRITECFP writes the Conduit Flow Package input file
%
% Exmample:
%   writeCFP(basenane,cfp)
%
% TO 090708

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename, '.',cfp.ext],'wt');

%0
%fprintf(fid,'# MATLAB writeCFP %s\n',datestr(now));
fprintf(    '# MATLAB writeCFP %s\n',datestr(now));

%1
fprintf(fid,'%s\n%10d\n','# mode',cfp.mode);

if cfp.mode==1 || cfp.mode==3  % something with tubes
    %2-4
    fprintf(fid,'%s\n','# item 2: data for mode 1 conduit pipe system');
    fprintf(fid,'%s\n','# item 3: number of nodes / tubes / layers');
    fprintf(fid,'%18i%12i%10i\n',cfp.NNODES,cfp.NPIPES,cfp.NLAY);

    %5-6
    fprintf(fid,'%s\n','# item 5: temperature in the pipes');
    fprintf(fid,'%25g\n',cfp.temperature);

    %7-8
    fprintf(fid,'# item 7: No  mc mr ml   Nb1 Nb2 Nb3 Nb4 Nb5 Nb6   tb1 tb2 tb3 tb4 tb5 tb6\n');
    fprintf(fid,['%12d%3d%3d%3d',...
                 '%6d%4d%4d%4d%4d%4d',...
                 '%6d%4d%4d%4d%4d%4d\n'],...
                  [cfp.Nvals(:,1:4) cfp.NB cfp.PB]');

    %9-12
    fprintf(fid,'%s\n','# item  9: elevation of conduit nodes. Two possibilites');
    fprintf(fid,'%s\n','# item 10: first: node #  elevation (1 line for each node)');
    fprintf(fid,'%s\n','# item 11: second: GEOHEIGHT (only one line used to assign constant value)');
    fprintf(fid,'%10d%15g\n',cfp.Nvals(:,[1 cfp.iGEOHEIGHT])');

    %13-14
    fprintf(fid,'%s\n','# item 13: surface dependent exchange (set 1) or constant exchange (set 0)');
    fprintf(fid,'%25d\n',cfp.SA_EXCHANGE);

    %15-22 (Mewton Raphson parameters)
    fprintf(fid,'%s\n','# item 15: criterion for convergence');
    fprintf(fid,'%28g\n',cfp.EPSILON);

    fprintf(fid,'%s\n','# item 17: maximum number for loop iterations');
    fprintf(fid,'%26d\n',cfp.NITER);

    fprintf(fid,'%s\n','# item 19: parameter of relaxation');
    fprintf(fid,'%25g\n',cfp.RELAX);

    fprintf(fid,'%s\n','# item 21: newton raphson print flag');
    fprintf(fid,'%25d\n',cfp.P_NR);

    %23-25 Tube data
    fprintf(fid,'%s\n','# item 23: data for tube parameters:');
    fprintf(fid,'%s\n','# item 24: %NO_DIAMETER TORTUOSITY RHEIGHT LCRITERY_P TCRITERY P');
    fprintf(fid,'%14d%6g%10g%10g%8g%10g\n',cfp.Pvals(:,[1 4:8])');

    %26-27 Node data
    fprintf(fid,'%s\n','# item 27: node    heads (if head unequal -1 the head is fixed)');
    fprintf(fid,'%14d%8g\n',cfp.Nvals(:,[1 cfp.iN_HEAD])');

    %28-29
    fprintf(fid,'%s\n','# item 29: node k-exchange terms for flow between continuum and pipe-network ');
    fprintf(fid,'%14d%8g\n',cfp.Nvals(:,[1 cfp.iK_EXCHANGE])');

end

if cfp.mode==2 || cfp.mode==3  % do something with turbulent layers
    %30-32
    fprintf(fid,'%s\n','# item 30. NCL total number of conduit layers');
    fprintf(fid,'%s\n','# item 31. Extra (useless) comment line');
    fprintf(fid,'%22d\n',sum(cfp.CL(cfp.CL~=0)));

    %33-34
    fprintf(fid,'%s\n','# item 33: CL whether or not this layer is a turbulent layertotal number of conduit layers');
    fprintf(fid,'%14d',find(cfp.CL~=0));
    fprintf(fid,'\n');

    %35-36
    fprintf(fid,'%s\n','# item 35: LTEMP mean layer temperature');
    fprintf(fid,'%22g\n',cfp.LTEMP);

    %37-39
    fprintf(fid,'%s\n','#Layer parameters');
    fprintf(fid,'%s\n','      #VOID LCRITREY_L TCRITREY_L');
    fprintf(fid,'%10g %6g %10g\n',cfp.Lvals(cfp.CL~=0,:)');
end

fclose(fid);
