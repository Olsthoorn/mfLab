function POut=line(o,pline)
%GRIDOBJ/LINE get info on all line pieces of polyline pline intersecting grid
%
% Examples:
%   POut=o.line(pline);
%
%   where pline is a polyline. The number of columns in pline must match
%   the dimension of the grid. so pline=[x y z]
%
%   Returns struct POut, in which each element holds the info of the line piece
%   passing through a single cell.
%
%  POut contains
%    x(2) y(2) z(2) ix iy iz I xm ym zm L
%
%   Look at P to see what it contains.
%
% ToDo: replace by lineObj (TO 130428)
%
% SEE ALSO: also cellIndex cellIndices xyzindex inpolygon gridObj gridObj.pline
%
% TO 100830

if size(pline,1)==1, pline=[pline;pline]; end  % need at least two points to make a vector

NLine=size(pline,1)-1; % number of lines in polyline

xGr = o.xGr; Nx=o.Nx;
yGr = o.yGr; Ny=o.Ny;
zGr = o.Z;
nDim = 3;

if size(pline,2)~=nDim % check the nr of columns in pline with the grid dimension
    error('mfLab:linegrid:misMatchPlineGrid',...
        'Dimensions in pline <<%d>> must match dimensions of grid in call to linegrid <<%d>>\n',size(pline,2),nDim);
end

if nDim>=1,
    xGr=sort(xGr(:),1,'ascend')';
    Nx=length(xGr)-1;
end
if nDim>=2,
    yGr=sort(yGr(:),1,'descend');
    Ny=length(yGr)-1;
end
if nDim>=3,

    % Make zGr full 3D array if not already given as such
    if isvector(zGr),  % is zGr a full 3D array or a vector?
        [~,~,zGr]=meshgrid(xGr,yGr,zGr); % make full 3D grid
    end
    %Nz=size(zGr,3);            % Nz+1,
    
    zGr=sort(zGr,3,'descend'); % Highest Z on top
    
    % replace zGr by relative z coordinates and remember original zGr
    zGrOld=zGr; % remember original Z-coordinates  !!
    
    for iz=1:size(zGr,3)
        zGr(:,:,iz) =(zGrOld(:,:,iz)-zGrOld(:,:,end))./(zGrOld(:,:,1)-zGrOld(:,:,end));  % relative z-coordinates
    end

    % replace pline(:,3) by relative z-coordinates and remember original values
    Ix=getindex(pline(:,1),xGr);
    Iy=getindex(pline(:,2),yGr);
    
    Iz=NaN(size(pline(:,1)));
    for i=1:size(pline,1)
        Iz(i)=getindex(pline(i,3),zGr(Iy(i),Ix(i),:));
        
        % pline(i,4) contain the relative coordinates
        pline(i,4)=(pline(i,3)-zGrOld(Iy(i),Ix(i),end))./(zGrOld(Iy(i),Ix(i),1)-zGrOld(Iy(i),Ix(i),end));
    end
end

POut=[];
for iL=1:NLine % number of  lines in polyline
    
    lambda=([]);
    warning('off','mflab:linegrid:divByZeroOk');  % division by zero is ok
    if nDim>=1, lambda=[0, 1, lambda, (xGr(:)'-pline(iL,1))/(pline(iL+1,1)-pline(iL,1))]; end %#ok
    if nDim>=2, lambda=[0, 1, lambda, (yGr(:)'-pline(iL,2))/(pline(iL+1,2)-pline(iL,2))]; end %#ok
    
    % zGr and pline(:,4) are both in relative z-coordinates
    if nDim>=3, lambda=[0, 1, lambda, (zGr(:)'-pline(iL,4))/(pline(iL+1,4)-pline(iL,4))]; end %#ok
    warning('on','mfLab:linegrid:divByZeroOk');   % division by zero
    
    lambda=lambda(lambda>=0 & lambda<=1); lambda(isnan(lambda))=[];
    
    if ~isempty(lambda),

        lambda=unique(lambda);

        if nDim>=1, x=pline(iL,1)+lambda.*(pline(iL+1,1)-pline(iL,1)); end
        if nDim>=2, y=pline(iL,2)+lambda.*(pline(iL+1,2)-pline(iL,2)); end
        if nDim>=3, z=pline(iL,4)+lambda.*(pline(iL+1,4)-pline(iL,4)); end
        
        %plot3(x,y,z,'ro'); % debug

        N=length(lambda)-1; % number of passed cells

        P=repmat(struct('ix',NaN','x',[NaN NaN]),N,1);
  
        % all lambda exist and are between 0 and 1, so no checks are required to
        % see of points fall inside grid, they all are guaranteed to do so
        for i=1:N
            if nDim>=1
                P(i).x     = x([i i+1])';
                P(i).xm    = mean(P(i).x);
                P(i).ix    = getindex(P(i).xm,xGr);
                P(i).L     = sqrt(diff(P(i).x)^2);
                P(i).idx   = P(i).ix;
            end
            if nDim>=2
                P(i).y     = y([i i+1])';
                P(i).ym    = mean(P(i).y);
                P(i).iy    = getindex(P(i).ym,yGr);
                P(i).L     = sqrt(diff(P(i).x)^2+diff(P(i).y)^2);
                P(i).idx   = Ny*(P(i).ix-1)+P(i).iy;
            end
            if nDim>=3
                zBot=   zGrOld(P(i).iy,P(i).ix,   end );          % abs z-coords
                H=-diff(zGrOld(P(i).iy,P(i).ix,[1 end]),1,3);     % abs H
                
                P(i).z     = zBot+z([i i+1])'*H;                    % abs z-coords
                P(i).zm    = mean(P(i).z);                          % bas z-coords
                P(i).L     = sqrt(diff(P(i).x)^2+diff(P(i).y)^2+diff(P(i).z)^2);
%                                zm=mean(z([i i+1]));                              % relative z-coords
                P(i).iz    = getindex(P(i).zm,zGrOld(P(i).iy,P(i).ix,:));   % using releative z-coords
                P(i).idx   = Ny*Nx*(P(i).iz-1)+Ny*(P(i).ix-1)+P(i).iy;
            end
        end

        POut=[POut;P]; %#ok
    end
end
end

function Iu=getindex(u,uGr)
% I=getindex(u,ugr)
% find the index of cell in which u lies given the gridlines UGr
% returns NaN if u is outside the range UGr
% UGr is sorted first
% TO 100830 120407

sgn=sign(uGr(end)-uGr(1));
if sgn<0
  u=-u(:); uGr=-uGr(:);
else
  u=+u(:); uGr=+uGr(:);
end

Iu=floor(interp1(uGr(:),1:length(uGr),u(:)));

Iu(u<=uGr(1))=1;
Iu(u>=uGr(end))=length(uGr)-1;

end