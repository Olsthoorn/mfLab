function [X,Y,d]=mf_GMLL2XY(varargin)
%MF_GMLL2XY computes XY and d in m of LL ralative to NC EC
%
%   d=true distance between points and center
%
% Example:
%   [X,Y,d]=mf_GMLL2XY(N,E,NC,EC);
%     XYd =mf_GMLL2XY(N,E,NC,EC);
%   [X,Y,d]=mf_GMLL2XY(NE,NCEC);
%    XYd  =mf_GMLL2XY(NE,NCEC);
%   [X,Y,d]=mf_GMLL2XY(EN)
%   [XYd]  =mf_GMLL2XY(EN);
% 
%   E and N are easting and northing 
%   If ECEC are not defiend (last call) than coordinates are relative to
%   first coordinate and d is cumulative distance, so EC is a path. 
%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110503

R2 = 6371007.2;  % The Earth's authentic radius (Geodetic Union);

switch length(varargin)
    case 1
         N = varargin{1}(:,1);
         E = varargin{1}(:,2);
         NC=N(1); N=N(2:end);
         EC=E(1); E=E(2:end);
    case 2
         N = varargin{1}(:,1);
         E = varargin{1}(:,2);
         NC= varargin{2}(:,1);
         EC= varargin{2}(:,2);
    case 3
         N = varargin{1};
         E = varargin{2};
         NC= varargin{3}(:,1);
         EC= varargin{3}(:,2);
    case 4
         N = varargin{1};
         E = varargin{2};
         NC= varargin{3};
         EC= varargin{4};
    otherwise
        error('incorrect number of input arguments');
end

lam=pi/180*E; lamc=pi/180*EC;
phi=pi/180*N; phic=pi/180*NC;

%X=pi*(E-EC)/180*R2*cos(NM); % NM=center of ifgure.
X=(lam-lamc)*R2 *cos(phic); % * cos(pi*cos(NC)/180);  % NC=center of ifgure.
Y=(phi-phic)*R2;

if nargout~=2
    % cos rule in 3D space vector coordinates
    x=R2*cos(lam).*cos(phi); xc=R2*cos(lamc).*cos(phic);
    y=R2*sin(lam).*cos(phi); yc=R2*sin(lamc).*cos(phic);
    z=R2         .*sin(phi); zc=R2*           sin(phic);
    r=sqrt(x.^2+y.^2+z.^2); rc=sqrt(xc.^2+yc.^2+zc.^2);
    if nargin==1
        x=[xc x]; y=[yc y]; z=[zc z]; r=[rc r];
       theta=acos(dot([x(1:end-1),y(1:end-1),z(1:end-1)],...
                      [x(2:end  ),y(2:end ), z(2:end  )])/(r(1:end-1).*r(2:end)));
        d=R2*theta;
        d=[0; d];    
    else
        dots=NaN(size(x));
        for i=1:size(x,1),
            dots(i)=dot([x(i),y(i),z(i)],...
                        [x(1),y(1),z(1)])./(r(i)*rc);
        end
        theta=acos(dots./(r*rc));
        d=R2*theta;
    end
end

if nargout==1
    X=[X Y d];
end

