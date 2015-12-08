function [s,x,q0] = model(varargin)
%[s,x] = model('zB',zB,'sL',sL,'sR',sR,'alpha',alpha,'kL',kL,'kR',kR,'L',L,'N',N,'x',x);

if nargin==0
    selftest();
    return;
end

[zB,varargin] = getProp(varargin,'zB',-20); % bottom of aquifer
[sL,varargin] = getProp(varargin,'sL', 0);  % head in m above datum (left)
[sR,varargin] = getProp(varargin,'sR', 0);  % head in m above datum (right)
[alpha,varargin] = getProp(varargin,'alpha',0.5); % location jump form kL to kR
[kL,varargin] = getProp(varargin,'kL',10);  % k for 0<x<alpha L
[kR,varargin] = getProp(varargin,'kR',20);  % k for alpha L < x < L
[L ,varargin] = getProp(varargin,'L',200);  % width of model
[x ,varargin] = getProp(varargin,'x',0:100:L);  % measurement locations
[N ,varargin] = getProp(varargin,'N',0.002); % recharge

if ~isempty(varargin)
    error(['Not all items of varargin are used\n',...
           'Remaining types: %s'],sprintfs(cellfun(@class',varargin,'UniformOutput',false)));
end

hL = sL-zB;  % head relative to bottom of aquifer left
hR = sR-zB;  % head relative to bottom of aquifer right

if hL<=0 || hR<=0
    error('%s: thickness of the aquifer must be >0.\nHowever: hL=%.4g, hR=%.4g, check zB=%.4g',...
        mfilename,hL,hR,zB);
end

left  = x<=alpha*L; % logicals, which x belong to left
right = x>=alpha*L; % logicals  which x belong to right

% q0 is outflow at left of model
q0 = ((hL^2-hR^2)-N*L^2*(alpha^2/kL+(1-alpha^2)/kR))/(2*L*(alpha/kL+(1-alpha)/kR));

% heads left and right of alpha L
h(left ) = sqrt(hL^2 - 2*q0/kL*   x(left)   - N/kL *      x(left ).^2 );
h(right) = sqrt(hR^2 + 2*q0/kR*(L-x(right)) + N/kR * (L^2-x(right).^2));

% head relative to datu
s = h + zB;


s = s(:); x=x(:);
end

function [s,x] =selftest()

zB  = -20;
sL  = 0;
sR  = 0;
alpha = 0.5;
kL  = 1;
kR  = 20;
L   = 200;
x   = 0:2:L;
N   = 0.002;
hL  = sL + zB;

[s,x,q0] = model('zB',zB,'sL',sL,'sR',sR,'alpha',alpha,'kL',kL,'kR',kR,'L',L,'N',N,'x',x);

% location of maximum head
xP = -q0/N;  % x at which h=max(h);

% maximum head relative to bottom of aquifer
if xP<alpha*L
    hP = sqrt(hL^2 + q0^2/(kL*N));
else
    hP = sqrt(hL^2 + q0^2/(kL*N)+ 2*q0*L/kR + N*L^2/kR);
end

% maximm head relative to datum
sP = hP+zB;

figure; hold on; xlabel('x [m]'); ylabel('elevation [m]');
plot(x,s,'b',[xP xP],[0 sP],'r--');
fprintf('x_hMax = %.2f m, sMax = %.2f m,  qL = %.4f , qR = %.4f\n',xP,sP, q0,N*L-q0);

end
