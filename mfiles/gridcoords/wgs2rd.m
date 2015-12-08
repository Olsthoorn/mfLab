function [X,Y]=wgs2rd(E,N,~)
%WGS2RD computes Dutch national coords (xRD,yRD) from GoogleEarth WGS coords (E,N)
%
% Example:
%   [xRD,yRD] = wgs2rd(E,N,verbose)
%   [xRD,yRD] = wgs2rd(EN);
%   [xRDyRD]  = wgs2rd(E,N);
%   [xRDyRG]  = wgs2rd(EN);
%
%   use arbitrary third argument to get verbose output
%   produced during internal optimization.
%
% SEE ALSO: rd2wgs rd2gw gw2rd getDinoXSec kmlpath kmlpath2rd
%
% TO 091123  (original fortran code obtained from Peter Vermeulen 2009)

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

verbose=0;

deltax=100;   % step size to compute derivative
deltay=100;   % same in y direction

switch nargin
    case 1
        if size(E,2)~=2
            error('%s: with one argument it must have 2 columns',mfilename);
        end
        N=E(:,2); E=E(:,1);
    case 2
        if ischar(N)
            verbose=1;
            if size(E,2)~=2
                error('%s: with one argument it muse have 2 columns',mfilename);
            end
            N=E(:,2); E=E(:,1);
        end
    case 3
        verbose=1;
    otherwise
        help wgs2rd()
        return;
end

if ~all(size(E)==size(N)), error('%s: size(E) must equal size(N)',mfilename); end
dims=size(E);

X=155000*ones(size(E)); Y=463000*ones(size(E)); 
dXY=NaN(2,numel(E));

while 1
    [E0,N0]=rd2wgs(X,Y);
    
    [dEdx dNdx]=rd2wgs(X+deltax,Y);
    dEdx=(dEdx-E0)/deltax;
    dNdx=(dNdx-N0)/deltax;
    
    [dEdy dNdy]=rd2wgs(X,Y+deltay);
    dEdy=(dEdy-E0)/deltay;
    dNdy=(dNdy-N0)/deltay;

    for i=1:numel(E)
        dXY(:,i)=...
            [dEdx(i) dEdy(i); dNdx(i) dNdy(i)]\...
            [E(i)-E0(i); N(i)-N0(i)];
    end
    
    DX=dXY(1,:)';
    DY=dXY(2,:)';
    [E1,N1]=rd2wgs(X+DX,Y+DY);
    
    % Verbose if nargin>2
    if verbose
     fprintf('DE=%12.4f DN=%12.4f  DX=%12.4g DY=%12.4g\n',...
         [E-E1 N-N1 DX DY]');
    end
    
    if all(abs(DX)<deltax) && all(abs(DY)<deltay),
       X=reshape(X,dims);
       Y=reshape(Y,dims);
       if nargout<2, X=[X(:) Y(:)]; end
       return;
    end
    X=X+DX;
    Y=Y+DY;
end
