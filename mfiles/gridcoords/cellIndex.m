function Idx=cellIndex(varargin)
%CELLINDEX get glocal index of 2D or 3D array using individual coordinate indices
%
% EXAMPLES:
%   Idx=cellIndex(Ix,Iy,Nx,Ny)
%   Idx=cellIndex(Ix,Iy,[Ny Nx])          e.g. Idx=cellIndex(Ix,Iy,size(A))    %2D
%   Idx=cellIndex([Ix Iy],[Ny Nx]);       e.g. Idx=cellIndex(IXYZ ,size(A))    %2D
%   Idx=cellIndex(Ix,Iy,Iz,Nx,Ny,Nz)
%   Idx=cellIndex(Ix,Iy,Iz,[Ny,Nx,Nz])
%       e.g. Idx=cellIndex(Ix,Iy,Iz,size(A)) %3D
%            Idx=cellIndex(Ix,Iy,Iz,gr.size) %3D
%   Idx=cellIndex([Ix Iy Iz],[Ny Nx Nz]); e.g.
%            Idx=cellIndex(IXYZ,size(A))     %3D
%            Idx=cellIndex(IXYZ,gr.size)     %3D
%   Idx=cellIndex([Ix Iy Iz],gridObj);
%   Idx=cellIndex( Ix,Iy,Iz, gridObj);
%
%   Compute the global cell indices Idx of the given indices and array dimension
%   Ix,Iy,Iz are arrays of equal size forming index tuples. Note the order of
%   the Nx and Ny in the forms when individual dimensions are specifed
%   instead of using size.
%   in form 3 and 6 [Ix Iy Iz] must form a vector of with 2 or 3 columns resp.
%
% SEE ALSO IdxMatlab2Modflow cellIndices xyzindex linegrid inpolygon
%
% TO 090317 090509 100521 100830


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

is3D=0;
switch nargin
    case 2 
        if strcmpi(class(varargin{2}),'gridObj')
            % cellIndex([Ix Iy Iz],gridObj)
            is3D=1;
            Ix=varargin{1}(:,1);
            Iy=varargin{1}(:,2);
            Iz=varargin{1}(:,3);
            Nyxz = varargin{2}.size;
            Ny = Nyxz(1);
            Nx = Nyxz(2);
            Nz = Nyxz(3);
        elseif size(varargin{2},2)==3
            % cellIndex([Ix Iy Iz],[Ny Nx Nz]) % notice order of sec arg
            % cellIndex([Ix Iy Iz],size(A))
            % cellIndex([Ix Iy Iz],gr.size)
            is3D=1;
            Ix=varargin{1}(:,1);
            Iy=varargin{1}(:,2);
            Iz=varargin{1}(:,3);
            Ny=varargin{2}(1);
            Nx=varargin{2}(2);
            Nz=varargin{2}(3);            
        else % cellIndex([Ix,Iy],[Ix Iy])
            Ix=varargin{1}(:,1); Ny=varargin{2}(1);
            Iy=varargin{1}(:,2); Nx=varargin{2}(2);
            if size(varargin{2},2)==3
                Iz=varargin{1}(:,3); Nz=varargin{2}(3);
                is3D=1;
            end
        end
    case 3 % cellIndex(Ix,Iy,[Ny Nx])
        Ix=varargin{1}; Iy=varargin{2};
        Ny=varargin{3}(1); Nx=varargin{3}(2);
    case 4
        if strcmpi(class(varargin{4}),'gridObj')
            % cellIndex(Ix,Iy,Iz,gridObj)
            is3D=1;
            Ix=varargin{1};    Iy=varargin{2};     Iz=varargin{3};
            NyNxNz = varargin{4}.size;
            Ny = NyNxNz(1); Nx=NyNxNz(2); Nz=NyNxNz(3);
        else
            if size(varargin{4},2)==1,
                % cellIndex(Ix,Iy,Nx,Ny);
                Ix=varargin{1}; Iy=varargin{2};
                Nx=varargin{3}; Ny=varargin{4};
            else
                % cellIndex(Ix,Iy,Iz,[Ny,Nx,Nz])
                % cellIndex(Ix,Iy,Iz,size(A));
                % cellIndex(Ix,Iy,Iz,gr.size)
                is3D=1;
                Ix=varargin{1};    Iy=varargin{2};     Iz=varargin{3};
                Ny=varargin{4}(1); Nx=varargin{4}(2);
                if length(varargin{4})<3, Nz=1; else Nz=varargin{4}(3); end
            end
        end
    case 6
        is3D=1;
        Ix=varargin{1}; Iy=varargin{2}; Iz=varargin{3};
        Nx=varargin{4}; Ny=varargin{5}; Nz=varargin{6};
    otherwise
        error('Inconsistent number of input arguments, use help cellIndex for help');
end

%% In case of 2D or 3D we always do this

% if any(Ix<1 | Ix>Nx), error('Ix out of bound'); end
% if any(Iy<1 | Iy>Ny), error('Iy out of bound'); end


if isempty(Ix) || isempty(Iy) || isempty(Iz)
    Idx=NaN;
    return;
end

if~is3D
    %5 Allow one or more of the I-vectors to be 1, we'll extend it to the
    %length of the other I-vector
    N = max([length(Ix) length(Iy)]);
    Ix(end+1:N)=Ix(end); % fill Ix vector to length of Iy if necessary
    Iy(end+1:N)=Iy(end); % fill Iy vector to length of Ix if necessary
else
    %5 Allow one or more of the I-vectors to be 1, we'll extend it to the
    %length of the other I-vector
    N = max([length(Ix) length(Iy) length(Iz)]);
    Ix(end+1:N)=Ix(end); % fill Ix vector to length of Iy|Iz if necessary
    Iy(end+1:N)=Iy(end); % same for Iy
    Iz(end+1:N)=Iz(end); % same for Iz
end

Idx=NaN(size(Ix(:)));

% process only the points that are within the grid
if ~is3D
    I = find( Ix>=1 & Ix<=Nx & Iy>=1 & Iy<=Ny);
    Idx(I)=Ny*(Ix(I)-1)+Iy(I);
else 
    I = find( Ix>=1 & Ix<=Nx & Iy>=1 & Iy<=Ny & Iz>=1 & Iz<=Nz);
    Idx(I)=Ny*Nx*(Iz(I)-1)+Ny*(Ix(I)-1)+Iy(I);
end
    
