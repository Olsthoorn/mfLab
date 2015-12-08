%readdxf
% specific for reading polylines from dxf file, only reads polylines
% TO 010825
%
fname='infgeb3.dxf';

fp=fopen(fname,'r');

pline=[]; ipoly=0;

while 1
s=fgetl(fp);

if s<0, break; end

if strcmp(s,'POLYLINE')
   ipoly=ipoly+1;
   pline(ipoly).x=[];
   pline(ipoly).y=[];
   j=0;
elseif strcmp(s,'AcDb2dVertex'),
   j=j+1;
   fgetl(fp); pline(ipoly).x(j)=sscanf(fgetl(fp),'%f');
   fgetl(fp); pline(ipoly).y(j)=sscanf(fgetl(fp),'%f');
end

end

fclose(fp);   

for k=1:2
%join polylines that have been missed
Rcrit=6; i=0;
while 1
   i=i+1; fprintf('\n%2d %2d',length(pline),i);
   if i==length(pline), break; end
   j=i;
   while 1
      j=j+1; fprintf('.')
      if j>length(pline), break; end
      r=sqrt((pline(i).x(end)-pline(j).x(1)).^2+(pline(i).y(end)-pline(j).y(1)).^2);
      if r<Rcrit
         pline(i).x=[pline(i).x,pline(j).x];
         pline(i).y=[pline(i).y,pline(j).y];
         pline(j)=[];
      else
	      r=sqrt((pline(i).x(end)-pline(j).x(end)).^2+(pline(i).y(end)-pline(j).y(end)).^2);
   	   if r<Rcrit
      	   pline(i).x=[pline(i).x,fliplr(pline(j).x)];
         	pline(i).y=[pline(i).y,fliplr(pline(j).y)];
	         pline(j)=[];
   	   else
      		r=sqrt((pline(i).x(1)-pline(j).x(1)).^2+(pline(i).y(1)-pline(j).y(1)).^2);
		      if r<Rcrit
	      	   pline(i).x=[fliplr(pline(i).x),pline(j).x];
   	      	pline(i).y=[fliplr(pline(i).y),pline(j).y];
	   	      pline(j)=[];
   		   else
			      r=sqrt((pline(i).x(1)-pline(j).x(end)).^2+(pline(i).y(1)-pline(j).y(end)).^2);
			      if r<Rcrit
         			pline(i).x=[pline(j).x,pline(i).x];
			         pline(i).y=[pline(j).y,pline(i).y];
         			pline(j)=[];
               end
            end
         end
      end
   end
end
end

boxx=pline(1).x([1,end]);
boxy=pline(1).y([1,end]);
for i=1:length(pline)
   pline(i).boxx=[min(pline(i).x),max(pline(i).x)];
   pline(i).boxy=[min(pline(i).y),max(pline(i).y)];
   boxx(1)=min([boxx(1),pline(i).boxx(1)]);
   boxx(2)=max([boxx(2),pline(i).boxx(2)]);
   boxy(1)=min([boxy(1),pline(i).boxy(1)]);
   boxy(2)=max([boxy(2),pline(i).boxy(2)]);
   % close polygons
   if ( pline(i).x(1) ~= pline(i).x(end) ) | ( pline(i).y(1) ~= pline(i).y(end) )
      pline(i).x(end+1)=pline(i).x(1);
      pline(i).y(end+1)=pline(i).y(1);
   end
end

% make mesh
dx= 5; boxx(1)=floor(boxx(1)/dx)*dx; boxx(2)=ceil(boxx(2)/dx)*dx;
dy=10; boxy(1)=floor(boxy(1)/dy)*dy; boxy(2)=ceil(boxy(2)/dy)*dy;

% cell face coordinates
x=boxx(1):dx:boxx(2);
y=boxy(1):dy:boxy(2);

% cel centere coordinates
xm=0.5*[x(1:end-1)+x(2:end)];
ym=0.5*[y(1:end-1)+y(2:end)];
[XM,YM]=meshgrid(xm,ym);

In=zeros(size(XM));
for i=1:length(pline)
   In=inpolygon(XM,YM,pline(i).x,pline(i).y,In,i);
end

figure
clr='rbgkmcy';
for i=1:length(pline)
   c=max(1,rem(i,length(clr)+1));
   plot(pline(i).x,pline(i).y,['-',clr(c)]); hold on
   I=find(In==i);
   xPoly(i)=XM(I(1)); yPoly(i)=YM(I(1));
   text(xPoly(i),yPoly(i),int2str(i));
end

OpenWatDef

save PlineIGEB2&3 pline x y xPoly yPoly xm ym  In Owater
