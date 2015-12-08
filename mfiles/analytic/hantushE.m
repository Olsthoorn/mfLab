function w=hantushE(u,zScrBot,zScrTop,r,z)
%HANTUSHE computes hantush but for partially penetrating well + self-test
%
%   Ref: Kruseman & De Ridder (1994), p162++
%      There is a solution of early times and one for late times.
%      The essence is that the partial penetration effect is transient.
%
% Example:
%    hantushE();   % selfTest
%    w=hantushE(u,zScrBot,zScrTop,r,z);
%
% See also: hantush hantushn
%
% TO 090101

% selfTest if nargin==0
if nargin==0, selftest(); return; end


zScrBot = abs(zScrBot);
zScrTop = abs(zScrTop);
z       = abs(z);

B = [(zScrBot+z)./r ...
     (zScrTop+z)./r ...
     (zScrBot-z)./r ...
     (zScrTop-z)./r];

% B is always numel(r),4)
% u is always numel(r),numel(t)
w = E(u,B(:,1))-E(u,B(:,2))+E(u,B(:,3))-E(u,B(:,4));

end

function w=E(u,B)
%E computes Hantush's modification of the Theis method for partially penetrating wells (ds_pp)
%
% Example:
%    w=E(u,r/Lambda)
%
%   Ref: Kruseman & De Ridder (1994), p162++
%      There is a solution of early times and one for late times.
%      The essence is that the partial penetration effect is transient.
%
%   The function is embedded in hantushE
%
%   E(u,b/r,d/r,a/r) = M(u,B1)-M(u,B2)+M(u,B3)-M(u,B4)
%
% with:
%   u = r^2S/(4kt)
%   Ss = S/D;
%   B1 = (b+a)/r
%   B2 = (d+a)/r
%   B3 = (b-a)/r
%   B4 = (d-a)/r
%   M(u,B) = int(exp(y)/a*erf(B sqrt(y)), u, Inf)
%
% See also: hantushE hantush hantushn
%
% TO 010409 121001

    if any(u)<=0;
        error('u must be >0');
    end

    B=B(:); % B is always [numel(r,1)]

    NT = size(u,2);  % times
    NP = size(B,1);  % points

    w = NaN(NP,NT);

    for iu=1:NT  % this is: for all times
        y=u(:,iu)*logspace(0,13,5000);
        dy = diff(y,1,2);
        arg = exp(-y)./y.*erf((B*ones(size(y(1,:)))).*sqrt(y));
        w(:,iu)= sum(arg(:,1:end-1).*dy/2 + arg(:,2:end).*dy/2,2);
    end

end

function w=selftest()
%SELFTEST selftest for Hantushs modification for partially penetrating wells
%
% Example:
%   hantushE();
%
% TO 130429

u=logspace(-6,1,71);

B=[0.01 0.05 0.1 0.5 1 2 3 4 5 6];  % see Kruseman & De Ridder (1971)

w=hantush(u,B);

figure;

n=0; leg={};

loglog(1./u,expint(u),'linewidth',3); n=n+1; leg{n}='Theis'; hold on

plot(  1./u,w); xlabel('1/u'); ylabel('W'); hold on

for i=1:length(B)
    n=n+1; leg{n}=sprintf('r/L=%g',B(i)); %#ok
end

grid on; set(gca,'xlim',[1e-1 1e6],'ylim',[1e-5 100],'color','none');

legend(leg,4);

title('Hantush modificationof Theis method for partially penetrating wells');

end