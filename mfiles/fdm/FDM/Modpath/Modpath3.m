function [ptcl]=modpath3(xMesh,yMesh,zMesh,por,qx,qy,qz,plist,times)
% [ptcl]=MODPATH3(xMesh,yMesh,zMesh,por,qx,qy,qz,plist,times)
% particle tracking according to MODPATH, see my PhD 1998
% ptcl = output particles as a struct
% xMesh = vector of cell boundaries
% yMesh = vector of cell boundaries
% qx = qrightface of model
% qy = qfrontface of model
% plist = list of particles [x y code] code is one of {-1,0,1} -1 is backtrack, 1=track, 0=skip this particle
% t = a timeseries or an end time, t=0 is start time assumed. t may be a vector of intermediate times
%
% TO 010811 010821 tested, see directory FDM\FD3D\BC testthesis1 and 2 and runmitiay 011005

Nx=length(xMesh)-1;
Ny=length(yMesh)-1;
Nz=length(zMesh)-1;
meshSize=[Nx,Ny,Nz];

% yMesh numbering and coordinates must be both positive in the same direction
if xMesh(end)<xMesh(1),
   xMesh=fliplr(xMesh(:)');
   qx= -qx(:,end:-1:1,:);
   qy=  qy(:,end:-1:1,:);
   qz=  qz(:,end:-1:1,:);
end
if yMesh(end)<yMesh(1),
   yMesh=flipud(yMesh(:));
   qx=  qx(end:-1:1,:,:);
   qy= -qy(end:-1:1,:,:);
   qz=  qz(end:-1:1,:,:);
end
if zMesh(end)<zMesh(1),
   zMesh=flipud(zMesh(:));
   qx=  qx(:,:,end:-1:1);
   qy=  qy(:,:,end:-1:1);
   qz= -qz(:,:,end:-1:1);
end
QX=zeros(Ny  ,Nx+1,Nz  ); QX(:,2:end-1,:)=qx;
QY=zeros(Ny+1,Nx  ,Nz  ); QY(2:end-1,:,:)=qy;
QZ=zeros(Ny  ,Nx  ,Nz+1); QZ(:,:,2:end-1)=qz;

dxMesh=abs(diff(xMesh));
dyMesh=abs(diff(yMesh));
dzMesh=abs(diff(zMesh));

% STARTING POSITION OF particles / watch for degenerate situations
ptcl=[];
for ip=1:size(plist,1)
   ptcl(ip).co(1,1)=plist(ip,1);
   ptcl(ip).co(2,1)=plist(ip,2);
   ptcl(ip).co(3,1)=plist(ip,3);
   ptcl(ip).dir=plist(ip,end);
   ptcl(ip).ic(1)=max(find(xMesh<=ptcl(ip).co(1)));				% col    of 1st coordinate along xMesh
   ptcl(ip).ic(2)=max(find(yMesh<=ptcl(ip).co(2)));				% row    of 2nd coordinate along yMesh
   ptcl(ip).ic(3)=max(find(zMesh<=ptcl(ip).co(3)));				% layer  of 3rd coordinage along zMesh
   
   ptcl(ip).pth =[0;ptcl(ip).co(1);ptcl(ip).co(2);ptcl(ip).co(3)];	% first point of flowpath
   ptcl(ip).tpth=[0;ptcl(ip).co(1);ptcl(ip).co(2);ptcl(ip).co(3)];	% first point of time point series
   if ptcl(ip).co(1)<min(xMesh) | ptcl(ip).co(1)>max(xMesh) |...
      ptcl(ip).co(2)<min(yMesh) | ptcl(ip).co(2)>max(yMesh) |...
      ptcl(ip).co(3)<min(zMesh) | ptcl(ip).co(3)>max(zMesh),
      ptcl(ip).stop=1;															% stop if outside mesh
   else
      ptcl(ip).stop=0;
	   if ptcl(ip).ic(1)==length(xMesh), ptcl(ip).ic(1)=ptcl(ip).ic(1)-1; end
   	if ptcl(ip).ic(2)==length(yMesh), ptcl(ip).ic(2)=ptcl(ip).ic(2)-1; end
   	if ptcl(ip).ic(3)==length(zMesh), ptcl(ip).ic(3)=ptcl(ip).ic(3)-1; end
   end
end


% FLOW PATH COMPUTATION
for ip=1:length(ptcl)														% compute the flow path
   TIMES=times;
   t=0;
   while ~isempty(TIMES) & t<TIMES(1) & ~ptcl(ip).stop
      ic=ptcl(ip).ic;														% cell with particle	  (col,row,layer)
      co=ptcl(ip).co;														% here is the particle (x,y,z)
      co1=[xMesh(ic(1)  );yMesh(ic(2)  );zMesh(ic(3)  )];		% cell corner with lowest  col row layer
      co2=[xMesh(ic(1)+1);yMesh(ic(2)+1);zMesh(ic(3)+1)];		% cell corner with highest col row layer
      d  =[dxMesh(ic(1)) ;dyMesh(ic(2)) ;dzMesh(ic(3))];			% size of cell [dx,dy,dz]
      q1=[QX(ic(2),ic(1)  ,ic(3)  ); QY(ic(2)  ,ic(1),ic(3)); QZ(ic(2),ic(1),ic(3)  )];		% q at left  and lower cell face
      q2=[QX(ic(2),ic(1)+1,ic(3)  ); QY(ic(2)+1,ic(1),ic(3)); QZ(ic(2),ic(1),ic(3)+1)];		% q at right and upper cell face
      V =por(ic(2),ic(1),ic(3))*prod(d);					% size of cell
      q =q1+(co-co1)./d.*(q2-q1);
      a=(q2-q1)/V;
      
      % compute times to reach other face of cell
		tau=zeros(1,4);
      for i=1:3
         flag(i)=abs(q1(i)-q2(i))<eps;				% v(i) is constant in cell
			next(i)=0;	tau(i)=TIMES(1)-t;			% default v==0;
         if q(i)>eps,
            	if q2(i)>0,			next(i)=1;
               	if flag(i),     			tau(i)=V/q(i)*(co2(i)-co(i))/d(i);
               	else, 						tau(i)=log(q2(i)./q(i))./a(i);
               	end
            	end
         elseif q(i)<-eps,
               if q1(i)<0,			next(i)=-1;
            		if flag(i),		  			tau(i)=V/q(i)*(co1(i)-co(i))/d(i);
            		else, 						tau(i)=log(q1(i)./q(i))./a(i);
	            	end
		         end
         end
      end
      tau(4)=TIMES(1)-t;
      
      % which of the computed times is the right one?
      I=find(min(tau)==tau); taumin=tau(I(1));
      
      % update cell
      if I(end)~=4,							% I(1)==4 means t=TIMES(1) is criterion, point within cell
	      for i=1:length(I)
   	      ptcl(ip).ic(I(i))=ptcl(ip).ic(I(i))+next(I(i));
         end
      end
      
      % move point to new location
      if I(end)==4							% I(1)==4 means t=TIMES(1) is criterion, point lies within cell
	      for i=1:3,								% move point
   	      if flag(i),     ptcl(ip).co(i)=co(i)+d(i)*q(i)/V*taumin;
      	   else            ptcl(ip).co(i)=co(i)+d(i).*q(i)./(q2(i)-q1(i)).*(exp(a(i)*taumin)-1);
         	end
	         if ptcl(ip).ic(i)==meshSize(i)+1 | ptcl(ip).ic(i)==0,
   	         ptcl(ip).stop=1;
      	   end
         end
         t=TIMES(1);				% point is a one of the prespecified times
   	   ptcl(ip).tpth=[ptcl(ip).tpth,[t;ptcl(ip).co]];
         TIMES(1)=[];
      else
	      for i=1:3,								% move point
   	      if flag(i),     ptcl(ip).co(i)=co(i)+d(i)*q(i)/V*taumin;
      	   else            ptcl(ip).co(i)=co(i)+d(i).*q(i)./(q2(i)-q1(i)).*(exp(a(i)*taumin)-1);
         	end
	         if ptcl(ip).ic(i)==meshSize(i)+1 | ptcl(ip).ic(i)==0,
   	         ptcl(ip).stop=1;
      	   end
         end
      	t=t+taumin;			% update time
	      ptcl(ip).pth=[ptcl(ip).pth,[t;ptcl(ip).co]];
      end
   end
end
