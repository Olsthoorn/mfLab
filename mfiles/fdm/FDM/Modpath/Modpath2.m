function [ptcl]=modpath2(xMesh,yMesh,por,qx,qy,plist,tmax)
% [ptcl]=MODPATH2(xMesh,yMesh,por,qx,qy,plist,tmax)
% 2D particle tracking according to MODPATH, see my PhD 1998
% ptcl = output particles as a struct
% xMesh = vector of cell boundaries
% yMesh = vector of cell boundaries
% qx = qrightface of model
% qy = qfrontface of model
% plist = list of particles [x y code] code is one of {-1,0,1} -1 is backtrack, 1=track, 0=skip this particle
% t = stop time, t=0 is start time assumed. t may be a vector of intermediate times
%
% TO 010811

% orientation of mesh??? This may have consequences for the computation below in y-direction

Nx=length(xMesh)-1;
Ny=length(yMesh)-1;
meshSize=[Nx,Ny];

% yMesh numbering and coordinates must be both positive in the same direction
if yMesh(end)<yMesh(1), yMesh=flipud(yMesh(:)); qx=flipud(qx); qy=flipud(qy); end

qx=[zeros(Ny,1),qx,zeros(Ny,1)];
qy=[zeros(1,Nx);qy;zeros(1,Nx)];

dxMesh=abs(diff(xMesh));
dyMesh=abs(diff(yMesh));

while tmax(1)==0 & length(tmax)>0; tmax(1)=[]; end; if isempty(tmax), tmax=365e3; end;

%locate particles / watch for degenerate situations
ptcl=[];
for ip=1:size(plist,1)
   ptcl(ip).co(1,1) =plist(ip,1); ptcl(ip).co(2,1)=plist(ip,2); ptcl(ip).dir=plist(ip,3);
   ptcl(ip).ic(1)=max(find(xMesh<=ptcl(ip).co(1)));
   ptcl(ip).ic(2)=max(find(yMesh<=ptcl(ip).co(2)));
   ptcl(ip).pth=[0;ptcl(ip).co(1);ptcl(ip).co(2)];
   if ptcl(ip).co(1)<min(xMesh) | ptcl(ip).co(1)>max(xMesh) | ptcl(ip).co(2)<min(yMesh) | ptcl(ip).co(2)>max(yMesh)
      ptcl(ip).stop=1;
   else
      ptcl(ip).stop=0;
	   if ptcl(ip).ic(1)==length(xMesh), ptcl(ip).ic(1)=ptcl(ip).ic(1)-1; end
   	if ptcl(ip).ic(2)==length(yMesh), ptcl(ip).ic(2)=ptcl(ip).ic(2)-1; end
   end
end

for ip=1:length(ptcl)
	time=0;	it=1; % initiate time
   while time<tmax(end) & ~ptcl(ip).stop
      ic=ptcl(ip).ic;		% cell with particle
      co=ptcl(ip).co;		% here is the particle
      co1=[xMesh(ic(1)  );yMesh(ic(2)  )];
      co2=[xMesh(ic(1)+1);yMesh(ic(2)+1)];		% upper right corner of cell
      d  =[dxMesh(ic(1));dyMesh(ic(2))];			% size of cell [dx;dy]
      q1=[qx(ic(2),ic(1)  );qy(ic(2)  ,ic(1))];			% q at left  and lower cell face
      q2=[qx(ic(2),ic(1)+1);qy(ic(2)+1,ic(1))];			% q at right and upper cell face
      V =por(ic(2),ic(1))*prod(d);					% size of cell
      q =q1+(co-co1)./d.*(q2-q1);
      a=(q2-q1)/V;
      
      % compute times to reach other face of cell
      for i=1:2
         flag(i)=abs(q1(i)-q2(i))<eps;				% v(i) is constant in cell
			next(i)=0;	tau(i)=tmax(it)-time;		% default v==0;
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
      
      % which of the computed times is the right one?
      I=find(min(tau)==tau); taumin=tau(I(1));
      
      % update cell
      for i=1:length(I)
         ptcl(ip).ic(I(i))=ptcl(ip).ic(I(i))+next(I(i));
      end
      
      % move point to new location
      for i=1:2,								% move point
         if flag(i),     ptcl(ip).co(i)=co(i)+d(i)*q(i)/V*taumin;
         else            ptcl(ip).co(i)=co(i)+d(i).*q(i)./(q2(i)-q1(i)).*(exp(a(i)*taumin)-1);
         end
         if ptcl(ip).ic(i)==meshSize(i)+1 | ptcl(ip).ic(i)==0,
            ptcl(ip).stop=1;
         end
      end
      time=time+taumin;
      ptcl(ip).pth=[ptcl(ip).pth,[time;ptcl(ip).co]];
   end
end
