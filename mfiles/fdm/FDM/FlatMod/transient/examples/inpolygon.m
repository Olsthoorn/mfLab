function In=inpolygon(X,Y,xp,yp,In,Code)
% In=inpolygon(X,Y,xp,yp,[[In],Code])
%see if points X,Y are in polygon xp,yp, if so return Code in the respective positions of In
% if Code is not given, 1 is used. If In is not given we start with zero matrix of size X
% size of X and Y  and   further that of xp and yp must match
% TO 010825

if nargin<6 | isempty(Code), Code=1; end	% Code given?
if nargin<5 | isempty(In  ), In=zeros(size(X)); end	% In given?
if xp(1)~=xp(end) | yp(1)~=yp(end); xp(end+1)=xp(1); yp(end+1)=yp(1); end  % closed polygon?
dx=diff(xp); dy=diff(yp);
xmin=min(xp); ymin=min(yp); xmax=max(xp); ymax=max(yp);	% box around polygon

I=find( X>=xmin & X<=xmax & Y>=ymin & Y<=ymax);				% consider only points in box
for i=1:length(I)
   ii=I(i);
   J=find((X(ii)>xp(1:end-1) & X(ii)<xp(2:end)) | (X(ii)<xp(1:end-1) & X(ii)>xp(2:end)));
   count=0;
   for j=1:length(J)
      jj=J(j);
      ym=yp(jj)+dy(jj)./dx(jj).*(X(ii)-xp(jj));
      if ym>=Y(ii), count=count+1; end
   end
   if (rem(count,2)), In(ii)=Code; end
end
