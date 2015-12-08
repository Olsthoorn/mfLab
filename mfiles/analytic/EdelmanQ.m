function [s,q]=EdelmanQ(x,t,a,kD,S,n,L)
%EDELMANQ computes transient 1D flow, constant profile
%
% Example:
%    [s,q]=EdelmanQ(x,t,a,kD,S [,n [,L]]);
%
% INPUT
%    x,t = distance and time vectors
%    a   = constant such that q(0)=a*t^((n-1)/2)
%    kD  = transmissivity, S=storage coefficient
%    n   = power>=0, if left out n=0
%    L   =distance at which head is fixed, if left out L=infinity
% OUTPUT
%    s   = head (Nx,Nt)
%    q   = discharge (Nx,Nt)
%
% See Also: EdelmanFunc EdelmanS
% 
% TO 100225

fprintf('EdelmanQ');

switch nargin
    case {0,1,2,3,4}, [s,q]=selftest; return; 
    case 5, n=0; L=Inf;
    case 6,      L=Inf;
    case 7, b=-sign(L); L=abs(L);
end

TOL=1e-6;

x=x(:)';  % [1,Nx]
t=t(:);   % [Nt,1]

t1=t(t<=0); t=t(t>0);

F=sqrt(S./(4*kD*t)); % [Nt,1]

Cs=ierfc(0,n)/ierfc(0,n-1)*2/sqrt(kD*S)*a*t.^( n   /2)*ones(size(x));
Cq=                                     a*t.^((n-1)/2)*ones(size(x));

s=E(x,F,n);     % [Nt,Nx]
q=E(x,F,n-1);   % [Nt,Nx]

N=100;
if L~=Inf
    for iL=1:N
        ds = E(2*L*iL+x,F,n)  +E(2*L*iL-x,F,n);
        dq = E(2*L*iL+x,F,n-1)-E(2*L*iL-x,F,n-1);
        s=s+b^iL*ds; q=q+b^iL*dq;
        fprintf('.');
        if max(abs(ds(:)))<TOL && max(abs(dq(:)))<TOL, 
            fprintf('%d iterations needed.\n',iL);
            break;
        end
    end
end
s=Cs.*s; if ~isempty(t1); s=[zeros(length(t1),length(x));s]; end
q=Cq.*q; if ~isempty(t1), q=[ones( length(t1),length(x));q]; end


function y=E(x,F,n)
% y=E(x,F,n)  - computes E(u)=ierfc(u,n)/ierfc(0,n), u=x*F
y=ierfc(F*x,n)/ierfc(0,n);

function y=ierfc(z,n)
% y=ierfc(z,n)
% repeated integral of complementary error function
% needed in transient 1D groundwater analysis
% TO 100222
%
switch n
    case -1
        y=2/sqrt(pi)*exp(-z.^2); return;
    case 0
        y=erfc(z); return;
    otherwise
        y=-z/n.*ierfc(z,n-1)+(1/2/n)*ierfc(z,n-2);
end

function [s,q]=selftest
% [s,q]=selftest -- in case nargin<5
kD=200; n=0; a=1; S=0.2; L=70; 
t=logspace(-3,3,51)';         % time vertical    Nt*1
x=logspace(-1,log10(L),20);   % space horizontal 1*Nx
[s,q]=EdelmanQ(x,t,a,kD,S,n,L); % EdelmanF(voor n=0)