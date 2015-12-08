function B=mf_Psi(varargin)
% MF_PSI computes stream function based on FLOWRIGHTFACE or FLOWFRONTFACE
%
% USAGE:
%    B=mf_Psi([gr],B);        % add Psi along x-axis for first zx plane
%    B=mf_Psi([gr],B,Iy)      % add Psi along x-axis for zx plane Iy
%    B=mf_Psi([gr],B,ix,'y')  % add Psi along y-axis for zy plane ix
%    B=mf_Psi([gr],B,Iy,'x')  % identical to mf_Psi(B,Iy);
%
% Add the gridObj gr (anywhere in the argument list)
% to correctly compute Psi for cases with confining beds
%
% TO 110319 110404 120902
% TO 130323 (had to dela with variable CBC labels during simulation)

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

% see if grid is given
[gr,varargin] = getType(varargin,'gridObj',[]);

% [avg,varargin] = getWord(varargin,{'avg','ave'});

[B, varargin] = getNext(varargin,'struct',[]);
if isempty(B)
    error('%s: No budget struct. It must be the first argument',mfilename);
end

[i  ,varargin] = getNext(varargin,'double',[]);

[dir,   ~    ] = getNext(varargin,'char','x');

[Ny,Nx,Nz] = size(B(1).term{1});

if isempty(i), i=1; end

%% flow into the model from below, add constant head inflows along bottom row

if dir =='x'
    fprintf('%s: Adding Psi along x-axis through row %d to the budget struct\n',mfilename,i);
else
    fprintf('%s: Adding Psi along y-axis through column %d to the budget struct\n',mfilename,i);
end

for it = length(B):-1:1
    
    switch dir
        case 'x'
            iLbl   = strmatchi({'FLOWRIGHT','FLOW RIGHT'},B(it).label);
            if ~iLbl
                error('%s: label FLOWRIGHT not in budget struct',mfilename);
            end
            
            % First compute Psi as if no confining beds exist
            Q = zeros(Nz+1,Nx-1,numel(i));
            Q(2:end,:,:) = XS(B(it).term{iLbl}(i, 1:end-1, end:-1:1));
            
            Q = flipdim(cumsum(Q,1),1);
            B(it).Psi = Q;
            
            % Add confining beds if they exist
            if ~isempty(gr) && any(gr.LAYCBD)
                Psi = zeros(gr.Nz+1,gr.Nx-1,size(B(it).Psi,3));
                Psi(gr.ITlay,:,:) = B(it).Psi(   1:end-1,:,:);
                Psi(gr.ITcbd,:,:) =       Psi(gr.ITcbd+1,:,:);
                B(it).Psi = Psi;
            end
    
        case 'y'
            iLbl   = strmatchi({'FLOWFRONT','FLOW FRONT'},B(it).label);
            if ~iLbl
                error('%s: label FLOWFRONT not in budget struct',mfilename);
            end
    
            Q = zeros(Nz+1,Ny-1,numel(i));
            Q(2:end,:,:) = YS(B(it).term{iLbl}(1:end-1, i,end:-1:1));
            Q=flipdim(cumsum(Q,1),1);
            B(it).PsiY = Q;
            
            if ~isempty(gr) && any(gr.LAYCBD)
                PsiY = zeros(gr.Nz+1,gr.Ny-1);
                PsiY(gr.ITlay,:,:) = B(it).PsiY(   1:end-1,:,:);
                PsiY(gr.ITcbd,:,:) =       PsiY(gr.ITcbd+1,:,:);
                B(it).PsiY = PsiY;
            end
    end
end
