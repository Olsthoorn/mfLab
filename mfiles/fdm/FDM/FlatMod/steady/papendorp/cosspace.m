function x=cosspace(x1,x2,n,side)
%function x=cosspace(x1,x2,n)
%function x=cosspace(x1,x2,n,'Left')		left is open
%function x=cosspace(x1,x2,n,'Right')		right is open
% verdeel het stuk x1,x2 in n stukjes met lengteverdeling volgens een cosinus
% TO 001005

if nargin>3
   s=upper(side(1));
else
   s='M';
end

if n>1, n=n-1; else n=1; end
dx=(x2-x1)/n;
xM=(x2-x1)/2;
switch s
case 'M'
   u=pi*([x1:dx:x2]-xM)/(x2-x1);
case 'L'
   u=0.5*pi*([x1:dx:x2]-x1)/(x2-x1);
case 'R'
   u=0.5*pi*([x1:dx:x2]-x2)/(x2-x1);
end
du=cos(0.5*(u(1:end-1)+u(2:end)));
dx=(x2-x1)*du/(sum(du));
x=x1+[0,cumsum(dx)];
