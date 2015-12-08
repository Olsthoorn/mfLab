function y=EdelmanFunc(u,n)
%EDELMANFUNC computes E_n=ierfc(u,n)/ierfc(0,n) used in 1D Edelman solutions
%
% Example:
%    [s,q]=EdelmanFunc(u,n);
%
% See also: EdelmanS EdelmanQ
%
% TO 100225

if  nargin==0, y=selftest;
else
    y=NaN(length(u(:)),length(n(:)));
    for in=1:length(n(:))
        y(:,in)=ierfc(u(:),n(in))/ierfc(0,n(in));
    end
end

function y=ierfc(z,n)
% y=ierfc(z,n) repeated integral of complementary error function
switch n
    case -1,   y=2/sqrt(pi)*exp(-z.^2); return;
    case 0,    y=erfc(z); return;
    otherwise, y=-z/n.*ierfc(z,n-1)+(1/2/n)*ierfc(z,n-2);
end

function y=selftest
% y=selftest -- in case nargin==0
% This generates the tables in of Huisman (1970), p42:47
n=[0 -1 1 2 3];
u=[0:0.002:0.1 0.11:0.01:1.1 1.12:0.01:1.96 2:0.1:3];
y=[u(:),EdelmanFunc(u,n)]; % EdelmanF(voor n=0)