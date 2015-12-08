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

% for all times
for iu = 1:NT
    y      = u(:,iu)*logspace(0,13,5000);
    dy     = diff(y,1,2);
    arg    = exp(-y)./y.*erf((B*ones(NT,1)).*sqrt(y));
    w(:,iu)= sum(arg(:,1:end-1).*dy/2 + arg(:,2:end).*dy/2,2);
end

