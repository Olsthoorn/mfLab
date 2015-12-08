function W=Wh(u,rho)
%W computes Hantush's well function for flow to a well in a semi-confined aquifer
%
% Example:
%    W=Wh(u,rho)
%
% INPUT:
%    u   = vector of r^2S/(4*kD*t)
%    rho = vector of r/lambda (lamda = sqrt(kD*c)
% OUTPUT
%    W   = hantush's well function for semi-confined well extraction
%
%    Can be used also for Theis' W(u) = W(u,0)
%    Steady steate K0(rho) = besselk(0,rho) = 2*W(0,rho)
%
% See also: hantushE hantusn
%
% TO 120114

%% Allow vector input to cover any combination of a list of u and of rho (r/L) values

if nargin==0; selfTest(); return; end

if nargin<2, rho=0; end  % in fact this is Theis

W=NaN(size(u(:)*rho(:)'));

for i=1:length(rho)
    for j=1:length(u)
        W(j,i)=wHantush(u(j),rho(i));
    end
end

%% nested function Wh computes a single value for sclar input
function w=wHantush(u,rho)
%HANTUSH computes well function by numerical integration over log(u)
%
% see also: hantush Wh Wh1 hantushn
%
% TO 120104

% Tackle special cases to remain exact where possible
if u==0                      % Steady state (De Glee)
    w=2*besselk(0,rho);
elseif rho==0                % Theis
    w=expint(u);  % theis
else                         % Hantush
    w=quadgk(@hantush,u,Inf);% quad integration allowing Inf boundary
end

end

function w=hantush(y)
%HANTUSH kernal of hantush's integral
%
%   nested-nested function hantush uses scope of mother function
%
% TO 120104

    rho2=rho/2;
    w=exp(-y-rho2^2./y)./y; % kernel of hantush's intergral
end

end

