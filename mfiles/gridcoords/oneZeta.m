function [ztaBot ztaTop] = oneZeta(gr,ZETA,ipln)
%ONEZETA extracts zeta planes (fresh-salt interfaces) from struct ZETA (SWI, salt water intrusion package)
%
% Example:
%     [ztaTop ztaBot] = oneZeta(gr,ZETA [,ipln])
%
%    Extracts zeta planes (interfaces) from struct ZETA (which has the same
%    structure as the that is read in with readBud().
%    ztaBot corresponds with the highest zeta and
%    ztaTop with the lowest.
%    In case the ztas overlaps in more than one layer (inversions)
%    and ztaBot the lowest.  ?
%    Without inversions, both are the same
%
% used in: VerifyWithAmwadu
%
% ToDo: I think this can be done smarter with bsxfun (TO 130428)
%
% SEE ALSO: readbud
%
% TO 120506
%
% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


if nargin<2, error('%s: Insufficient input arguments',mfilename); end

if ~strcmp(class(gr),'gridObj'),
    error('%s: first argument must be gridObj.',mfilename);
end

if nargin<3, ipln=1; end

Nt = length(ZETA); Nxy=gr.Nx*gr.Ny;
delta = 1e-2;

for it=1:Nt
    % valid(:,:,iL) is true of interface is between top and bottom of layer iL
    valid = (ZETA(it).term{ipln}>gr.ZBlay+delta & ZETA(it).term{ipln}<gr.ZTlay-delta);
    
    % because an ipln may have a legal interface in every layer,
    % find lowest position if this ipln
    ztaBot = NaN(gr.Ny,gr.Nx,Nt); % 3D Ny*Nx*Nt (not Nz) to store ztaBot
    for iL=1:gr.Nlay
        I = find(valid(:,:,iL));
        ztaBot(I+Nxy*(it-1)) = ZETA(it).term{ipln}(I+Nxy*(iL-1));
    end
    
    % find highest position of this ipln
    if nargout>1
        ztaTop = NaN(gr.Ny,gr.Nx,Nt); % 3D Ny*Nx*Nt (not Nz) to store ztaTop
        for iL=gr.Nlay:-1:1
            I = find(valid(:,:,iL));
            ztaTop(I+Nxy*(it-1)) = ZETA(it).term{ipln}(I+Nxy*(iL-1));
        end
    end
end

