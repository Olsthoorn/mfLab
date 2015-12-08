function HFB=mf_setHFB(varargin)
%MF_SETHFB Sets Horizontal Flow Barrier
%
% Example1: (3 arguments)
%   HFB=mf_setHFB(zone,ILay,c)
%
%   zone is an array with >0 inside HFB and 0 outside
%   ILay the number of layers with HFB
%   c the hydraulic resistance per layer
%   if c is vector, than the values are assumed to belong to the layers
%   in the sequence given by ILay.
%   if length(ILay)>length(c), the last c is used for the rest of the
%   layers
%
% Example2: (4 arguments)
%   HFB=mf_setHFB([xP yP[cP]],xGr,yGr,ILay,c)
%
%   xP yP vertices of HFB,
%   cP resistance per section, will be overwritten by c
%   dims dimension of array like (size(IBOUND))
%   ILay and c as specified above
%
% Example3: (4 arguments)
%   HFB=mf_setHFB([xP yP cP],xGr,yGr,LLay)
%
%   XP and YP as specified above,   
%   cP is the hydraulic resistance at vertices between vertices
%   If last c equals NaN
%   dims and ILay as specified above
%   loop is closed if last cP(end) >0
%   used
%
%   all cP<0 will be set to 0
%   all c==0 will be passed over without output
%
%  if you need to specify different c for different vertices and edges
%  you may call this function repeately with the correct arguments and
%  values.
%
%  TO 110428

switch nargin
    case 3   % USAGE:  HFB=mf_setHFB(zone,ILay,c)
        zone= varargin{1};
        ILay= varargin{2};
        c   = varargin{3};
    case 4   % USAGE   HFB=mf_setHFB([xP yP cP],xGr,yGr,ILay)
        P   = varargin{1}; % cP will be used
        xGr = varargin{2};
        yGr = varargin{3};
        ILay= varargin{4};
    case 5   % USAGE   HFB=mf_setHFB([xP yP cP],xGr,yGr,ILay c)
        P   = varargin{1}(:,1:2); % cP will not be used !!
        xGr = varargin{2};
        yGr = varargin{3};
        ILay= varargin{4};
        c   = varargin{5};
    otherwise
end

if exist('zone','var')
    Ix=find(diff(zone,1,2)~=0);
    RCL=cellIndices(Ix,size(zone),'RCL');
    HFBx = [RCL(:,1),RCL(:,2),RCL(:,1),RCL(:,2)+1];
    
    Iy=find(diff(zone,1,1)~=0);
    RCL=cellIndices(Iy,size(zone),'RCL');
    HFBy = [RCL(:,1),RCL(:,2),RCL(:,1)+1,RCL(:,2)];
    
    n=size(HFBx,1)+size(HFBy,1);
    
    hfb=[HFBx ; HFBy];
    
    HFB=repmat([NaN(n,1) NaN(size(hfb)) NaN(n,1)],[length(ILay),1]);

    u=ones(n,1);
    for iL=1:length(ILay)
        IL=(iL-1)*n+(1:n)';
        HFB(IL,:)=[u*ILay(iL) hfb u*c(min(iL,length(c)))];
    end
else
    hfb=gethfb(P,xGr,yGr);

    [n,m]=size(hfb); u=ones(n,1);

    HFB=repmat(NaN(n,6),[length(ILay),1]); % 6=[L R1 C1 R2 C3 c]

    for iL=1:length(ILay)    % for each layer as specified in call
        IL=(iL-1)*n+(1:n)';  % indices into large HFB array for this layer
        HFB(IL,1)=ILay(iL);    % Layer number in front
        HFB(IL,2:m+1)=hfb;   % hfb next to it
        if exist('c','var'),
            HFB(IL,end)=u*c(min(iL,length(c))); % place or overwrite c
        end
    end
end

%% Remove sections with negative or zero resistance values

HFB=HFB(HFB(:,end)>0,:);   % move where c<=0
HFB(:,end)=1./HFB(:,end);  % use reciprocals for MODFLOW

function hfb=gethfb(P,xGr,yGr)
% GET_hfb: Get the locations of the HFB with c given points P grid xGR,yGr
%
% USAGE:
%   hfb=get_hfb([xP yP [cP]],xGr,yGr)
%
% TO 110529

if size(P,1)<2,
    ME='mf_setHFB:gethfb';
    fprintf('Warning polygon must have at least 2 points: %s\n',ME);
    return;
end % need at least to points

if size(P,2)>2
    if P(end,3)>0 % close HFB boundary
        fprintf('mfLAB:mf_setHFB:gethfb, HFB boundary closed as cP(end)>0\n');
        P(end+1,:)=P(1,:);
    end
end

xP=P(:,1); yP=P(:,2);

[xGr,yGr,xm,ym,Dx,Dy,Nx,Ny]=modelsize(xGr,yGr);
[XM,YM]=meshgrid(xm,ym);

hfb=[];  % [R C R C]  % grows when parsing polyline

x1=xP(1);  % start with first point, xP yP need not be closed !
y1=yP(1);
    
for i=2:length(xP)  % run over each line piece of polygon
    
    x2=xP(i); dx=x2-x1; Ix=find(xm>=min(x1,x2) & xm<max(x1,x2));  % number of xm lines crossed
    y2=yP(i); dy=y2-y1; Iy=find(ym>=min(y1,y2) & ym<max(y1,y2));  % number of ym lines crossed
    
% Looking in x-direction: we will have Iy x values along line
    if ~isempty(Iy), 
        x=(ym(Iy)-y1)*dx/dy+x1; % Iy in number

        Jx=find(diff((XM(Iy,:)<x*ones(1,Nx)),1,2)); % index within submatrix for cells just left of line

        RCL=cellIndices(Jx,[length(Iy),Nx],'RCL');  % RC indices within submatrix
    
        hfbx=[min(Iy)-1+RCL(:,1),RCL(:,2),...
              min(Iy)-1+RCL(:,1),RCL(:,2)+1];
          
       % plot(x,ym(Iy),'.m');  % debug

    else
        hfbx=[];
    end

% Looking in y-direction
    if ~isempty(Ix),
        y=(xm(Ix)-x1)*dy/dx+y1; % Iy in number
    
        Jy=find([diff((YM(:,Ix)>ones(Ny,1)*y),1,1);zeros(size(Ix))]); % index within submatrix for cells just above line
    
        RCL=cellIndices(Jy,[Ny,length(Ix)],'RCL');    % indices within submatrix
  
        hfby=[RCL(:,1)  ,min(Ix)-1+RCL(:,2),...
              RCL(:,1)+1,min(Ix)-1+RCL(:,2)];
          
        % plot(xm(Ix),y,'.');  % debug
        % plot(xm(hfby(:,2)),ym(hfby(:,1)),'m'); % debug
    else
        hfby=[];
    end
      
    if size(P,2)<3
        hfb = [hfb; hfbx ; hfby];
    else
        if ~isempty(hfbx)
            ux=ones(size(hfbx(:,1)));
            hfb=[hfb; [hfbx ux*P(i,3)]];
        end
        if ~isempty(hfby)
            uy=ones(size(hfby(:,1)));
            hfb=[hfb; [hfby uy*P(i,3)]];
        end
    end
    
    x1=xP(i); y1=yP(i);
end