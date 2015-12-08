function [Connections,Pipes,Geoheight,Fixedhd,branch]=poly3DXgrid(xGr,yGr,zGr,branch)

% POLY3DXGRID: Intersect a 3D polyline with a 3D grid to defined tunnel and
%   pipe interactions. Meant for CFP (conduit flow package of MODFLOW 2005)
%
% USAGE:
%    branch=poly3DXgrid(xGr,yGr,zGr,branch)
%
%   xyw=[xNd,yNd,zNd] is a line whos vertices are nodes and whose segments may denote
%   pipes or any internode connections.
%   xGr, yGr, zGr is the grid
%   the output branch additonally contains branch.xgrid,branch.ygrid,branch.zgrid,
%   poly3DXgrid;                                    % Selftest
%   branch=poly3DXgrid(xGr,yGr,zGr,branch);
%
%   TO 070521 070616 090624 090718


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin==0; [Connections,Pipes,Geoheight,Fixedhd,branch]=selftest; return; end

Nx=length(xGr); Ny=length(yGr);
lastpipe=0;

%% in case no polyline given, get it from screen, using z=0
if nargin<=3
        [branch.x,branch.y]=getline(gca); % from the image toolbox
         branch.z=zeros(size(branch.x));
end

%% Loop the branches

for iBr=1:length(branch)
    
    %% Remove points that are double
    I=find(abs(diff(branch(iBr).x))<=eps & abs(diff(branch(iBr).y))<=eps & abs(diff(branch(iBr).z))<=eps);
    if ~isempty(I), branch(iBr).x(I+1)=[]; branch(iBr).y(I+1)=[]; branch(iBr).z(I+1)=[]; end

    hold on
    plot3(branch(iBr).x,branch(iBr).y, branch(iBr).z,'--bs');

%% Intersect grid with branch polyline

    % start with intial branch point
    xNd=branch(iBr).x(1); % new node coordinats after intersection with grid 
    yNd=branch(iBr).y(1); 
    zNd=branch(iBr).z(1);
    fHd=branch(iBr).FH(1); % initialize fixed heads for subnodes
    
    branch(iBr).Pipe=[];    % initialize pipes between cells, use Tubes for realisty indep of mesh

    % ==== intersect the conduits the FD grid
    Np=length(branch(iBr).x);
    if Np>1
        for ip=1:Np-1  % all vectors in interface points
            x1=branch(iBr).x(ip); x2=branch(iBr).x(ip+1);
            y1=branch(iBr).y(ip); y2=branch(iBr).y(ip+1);
            z1=branch(iBr).z(ip); z2=branch(iBr).z(ip+1);

            lam=unique([(xGr(:)-x1)./(x2-x1); (yGr(:)-y1)./(y2-y1); (zGr(:)-z1)./(z2-z1)]);
            I=find(lam>0 & lam<1); % bramch intersection with grid not including the branch point
            
            % intermediate points between (not including the branch given
            % intermediate points) This is guaranteed by I=find(lam>0 & lam<1)
            xt=x1+lam(I)*(x2-x1); xNd=[xNd; 0.5*(xt(1:end-1)+xt(2:end)); x2];
            yt=y1+lam(I)*(y2-y1); yNd=[yNd; 0.5*(yt(1:end-1)+yt(2:end)); y2];
            zt=z1+lam(I)*(z2-z1); zNd=[zNd; 0.5*(zt(1:end-1)+zt(2:end)); z2];
            
            %xx=[xx; xt; x2];  % branch point and cell face intersections
            %yy=[yy; yt; y2];
            %zz=[zz; zt; z2];
            
            N=length(I);
            
            fHd=[fHd; -1*ones(N-1,1); branch(iBr).FH(ip+1);];

            % Tubes properties between these two branch points
            branch(iBr).Pipe=[ branch(iBr).Pipe;...
                [lastpipe+(1:N)', zeros(N,2), ...           % PipeNr FromND ToNd
                  ones(N,1)* branch(iBr).Tube(ip,:)] ];     % Apply tube properties from database
            lastpipe=lastpipe+N;                            % Pipe number counter
        end
        
    end

%%  compute the distance between the new nodes
    branch(iBr).xNd=xNd; branch(iBr).dx=diff(xNd);
    branch(iBr).yNd=yNd; branch(iBr).dy=diff(yNd);
    branch(iBr).zNd=zNd; branch(iBr).dz=diff(zNd);
    branch(iBr).fHd=fHd;
    branch(iBr).dL=sqrt(branch(iBr).dx.^2+branch(iBr).dy.^2+branch(iBr).dz.^2);
    
    %branch(iBr).xx=xx;
    %branch(iBr).yy=yy;
    %branch(iBr).zz=zz;

%% Number of nodes in this branch
    N=length(branch(iBr).xNd);
    
%%  compute grid coordinates of new node, including grid cell number
    branch(iBr).Ix=zeros(N,1);
    branch(iBr).Iy=zeros(N,1);
    branch(iBr).Iz=zeros(N,1);
    branch(iBr).I =zeros(N,1);
    
    for i=1:length(branch(iBr).Ix)
        if xGr(end)>xGr(1), branch(iBr).Ix(i)=find(branch(iBr).xNd(i)>=xGr,1,'last');
        else                branch(iBr).Ix(i)=find(branch(iBr).xNd(i)<=xGr,1,'last'); end
        if yGr(end)>yGr(1), branch(iBr).Iy(i)=find(branch(iBr).yNd(i)>=yGr,1,'last');
        else                branch(iBr).Iy(i)=find(branch(iBr).yNd(i)<=yGr,1,'last'); end
        if zGr(end)>zGr(1), branch(iBr).Iz(i)=find(branch(iBr).zNd(i)>=zGr,1,'last');
        else                branch(iBr).Iz(i)=find(branch(iBr).zNd(i)<=zGr,1,'last'); end
        branch(iBr).I(i) =(branch(iBr).Iz(i)-1)*Ny*Nx+(branch(iBr).Ix(i)-1)*Ny+branch(iBr).Iy(i);
    end

%% compute difference in grid coordinates between branch point neighbors
    branch(iBr).idxyz=[branch(iBr).Ix(2:end)-branch(iBr).Ix(1:end-1),...
                       branch(iBr).Iy(2:end)-branch(iBr).Iy(1:end-1),...
                       branch(iBr).Iz(2:end)-branch(iBr).Iz(1:end-1)];

%% Compute pipe numbers connected to branch node, called PpLR (left and
%  right pipe for each node in a single branch
    branch(iBr).NbLR=zeros(N,2);
    branch(iBr).NbLR(2:end  ,1)=branch(iBr).I(1:end-1);
    branch(iBr).NbLR(1:end-1,2)=branch(iBr).I(2:end  );
    

%% Compute pipe numbers connected to branch node, called PpLR (left and
%  right pipe for each node in a single branch
    branch(iBr).PpLR=zeros(N,2);
    branch(iBr).PpLR(2:end  ,1)=branch(iBr).Pipe(:,1);
    branch(iBr).PpLR(1:end-1,2)=branch(iBr).Pipe(:,1);

%% Compute To and From nodes of pipes (not used by CPF)
    % Number of pipes in a single branch should equal number of nodes minus 1
    branch(iBr).Pipe(:,2)=branch(iBr).I(1:end-1);
    branch(iBr).Pipe(:,3)=branch(iBr).I(2:end  );
                   
end  % branch


%% Join all nodes in a single connections list

%               NdNr         Nb1 .. Nb 6                   Pipe1 .. Pipe6
Connections=[branch(1).I, branch(1).Ix, branch(1).Iy, branch(1).Iz, ...
    branch(1).NbLR, zeros(size(branch(1).I,1),4),...  % connected nodes
    branch(1).PpLR, zeros(size(branch(1).I,1),4),...  % connected pipes
    branch(1).fHd  ,branch(1).zNd] ;   % fixed heads, geoheight
Pipes=branch(1).Pipe;
Geoheight=[branch(1).I branch(1).zNd];

NbCol=4+(1:6); PpCol=NbCol(end)+(1:6); FHCol=PpCol(end)+1; GHCol=FHCol(end)+1;
Connections(:,NbCol)=sort(Connections(:,NbCol),2,'descend');
Connections(:,PpCol)=sort(Connections(:,PpCol),2,'descend');

for iBr=2:length(branch)
    J=intersect(Connections(:,1),branch(iBr).I);  % find already existing nodes in new set
    for j=1:length(J)   % for all repeated nodes
        % Prepare addition of neighbors to existing connections
        ii=find(Connections(:,1)==J(j));           % locate repeated node in existing list
        jj=find(branch(iBr).I   ==J(j));           % locate existing node in new list 
 
        k=find(Connections(ii,NbCol)~=0,1,'last');  % find last occupied position in neighbor list (6 positions)
        m=find(Connections(ii,PpCol)~=0,1,'last');  % find last occupied position in pipe list (6 positions)
        
        % Add the neighbors with repreated node numbers from the new set of
        % connections to previous set
        % CONNECTED NEIGHBORS
        Connections(ii,k+NbCol(1:2))=sort(branch(iBr).NbLR(jj,:),2,'descend'); % put at first non blank position
        % CONNECTED PIPES
        Connections(ii,m+PpCol(1:2))=sort(branch(iBr).PpLR(jj,:),2,'descend'); % put at first non blank position
    end

    [Dummy,M]=setdiff(branch(iBr).I,J);
    Connect=[branch(iBr).I(M), branch(1).Ix(M), branch(1).Iy(M), branch(1).Iz(M),...
        branch(iBr).NbLR(M,:), zeros(length(M),4),...
        branch(iBr).PpLR(M,:), zeros(length(M),4),...
        branch(iBr).fHd(M), branch(iBr).zNd(M)];
    
    Connect(:,NbCol)=sort(Connect(:,NbCol),2,'descend'); % move zeros to end
    Connect(:,PpCol)=sort(Connect(:,PpCol),2,'descend'); % move zeros to end

    Connections=[Connections; Connect];

    Pipes=[Pipes;branch(iBr).Pipe];
end

Fixedhd  =Connections(:,[1,FHCol]);
Geoheight=Connections(:,[1,GHCol]);
        



function [Connections,Pipes,Geoheight,Fixedhd,branch]=selftest
xGr=0:100:1000;
yGr=1000:-100:0;
zGr=0:-10:-100;

%% Branch 1 for self test
branch(1).x=[ 10 200 300 500 700 600]+10*(rand-0.5);
branch(1).y=[800 750 500 100 200 300]+10*(rand-0.5);
branch(1).z= [-10 -90 -30 -60 -10 -50]+10*(rand-0.5);
branch(1).FH=[-1  -1  -1   10  -1  -1];
% Specify the fixed heads right away as they may not vary with stress periods

% must still include wall conductivity
branch(1).P=zeros(length(branch(1).x)-1,5);  % [Diam Tort Roughness LRey TRey]
branch(1).Tube(:,:)=[0.20  1.5  0.001 1000 2000 5; ...
                     0.25  1.2  0.001 1000 2000 5; ...
                     0.30  1.5  0.001 1000 2000 5; ...
                     0.35  1.0  0.001 1000 2000 5; ...
                     0.40  1.5  0.001 1000 2000 5];
    
%% Branch 2 for self test
branch(2).x=[300 350 425 725]+10*(rand-0.5);
branch(2).y=[500 150 200 250]+10*(rand-0.5);
branch(2).z=[-30 -60 -10 -60]+10*(rand-0.5);
branch(2).FH=[-1 -1  -1    5];

branch(2).P=zeros(length(branch(2).x)-1,5);  % [Diam Tort Roughness LRey TRey kWall]
branch(2).Tube(:,:)=[0.20  1.5  0.001 1000 2000  5; ...
                     0.35  1.0  0.001 1000 2000  5; ...
                     0.50  1.5  0.001 1000 2000  5];

%% show branches              
clrs='rbgkmcy';

for iBr=1:length(branch)
    h=plot3(branch(iBr).x,branch(iBr).y,branch(iBr).z,'ro-'); hold on
    set(h,'MarkerFaceColor',clrs(iBr));
    set(gca,'xlim',xGr([1,end]),'ylim',yGr([end,1]),'zlim',zGr([end,1]));
    xlabel('x'); ylabel('y'); zlabel('z');
end

drawnow

%% Compute
[Connections,Pipes,Geoheight,Fixedhd,branch]=poly3DXgrid(xGr,yGr,zGr,branch);

%% Show
for iBr=1:length(branch)
    plot3(branch(iBr).xNd,branch(iBr).yNd,branch(iBr).zNd, [clrs(iBr), 'o-']);
end
