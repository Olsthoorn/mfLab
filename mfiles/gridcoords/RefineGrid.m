function ANew=RefineGrid(AOld,splitArray,code)
%REFINEGRID refines rows or cols of a FDM model array according to splitArray
%
% Example:
%    ANew=refineRC(AOld,splitArray,code)
%
%    splitArray has the form
%        [i1 i2 i3 i4;
%         N1 N2 N3 N4];
%
%   where i refers to the old column and N the number by which the column must be
%   subdivided.
% code is either 'x', 'y', 'z' 'X' 'Y  or 'Z' to indicated the direction
%   along which the splitArray is to be applied and to indicate whether the
%   array to be split is a grid array or a cell array. A grid array concerns
%   the cell face coordinates and non-grid array the cell center coordinates.
%   The grid array divides the distance betwee adjacent lines according to
%   splitArray.
% A second arbitrary letter in code indicates that the cell properties are
%   to be subdivided and not to be taken as they are.
%   This is done for instance with the cell-size arrays.
% Repeat the procedure for each direction sepaprately. To split an already
%   splitted array, just run the procedure again with the new array and the
%   correct spitarray for the new array. This may get complex though.
%   see explanation at the end of this file.
%
% These function are superseded by methods of the gridObj
%
% SEE ALSO: gridObj modelObj cutoutBCN cutoutXTA RefineBCN RefineGrid JoinBCN JoinLayers removeCBD
%
% TO 100602

code=code(1); % only one-letter code

%% Complete the splitArray to hold all indices, including where no splitting takes place
    
% Permute to aline with x-axis so that all directions can be treated equally
switch lower(code)  % works for 'x' and 'X', 'y' and 'Y', 'z' and 'Z'
    case 'x', fprintf('Refining array in x-direction\n');
              AOld=permute(AOld,[1,2,3]);
    case 'y', fprintf('Refining array in y-direction\n');
              AOld=permute(AOld,[2,1,3]);
    case 'z', fprintf('Refining array in z-direction\n');
              AOld=permute(AOld,[1,3,2]);        
    otherwise
        error('Unknown code, see help for this function');
end


%% Complete the splitArray
if max(splitArray(1,:))>size(AOld,2)
    error('Highest value in top layer of splitArray (%d) exceeds size of Aold (%d)',max(splitArray(1,:)),size(AOld,2));
end

% We only need split(1,size(AOld,2) to tell into how many fields the loca
% pint is divided
split=ones(1,size(AOld,2)); if code=='X' || code=='Y'|| code=='Z', split=split(1:end-1); end
split(splitArray(1,:))=splitArray(2,:);

%% Refinde grid arrays
if code == 'X' || code=='Y' || code=='Z'
    First = AOld(:,1,:);     % remember first
    AOld  = diff(AOld,1,2);  % treat rest as cell values
end

%% Treat as if cell values
ANew = NaN(size(AOld,1),sum(split),size(AOld,3));

k=0;
for j=1:length(split) % = Nr of old cols, row, layers
   if code == 'X' || code=='Y' || code=='Z'
       for jj=1:split(j)
           ANew(:,k+jj,:)=AOld(:,j,:)/split(j); % distance
       end
   else
       for jj=1:split(j)
           ANew(:,k+jj,:)=AOld(:,j,:);  % just copy old value
       end
   end
   k=k+split(j);
end

if code =='X' || code=='Y' || code=='Z', 
    ANew = horzcat(First,ANew);
    ANew(:,2:end,:) = repmat(First,[1,size(ANew,2)-1,1]) + cumsum(ANew(:,2:end,:),2); 
end

%% Permute back
switch lower(code)
    case 'x', ANew=permute(ANew,[1,2,3]);
    case 'y', ANew=permute(ANew,[2,1,3]);
    case 'z', ANew=permute(ANew,[1,3,2]);
end

%% Procedure
% Splitting a layer according to splitArray:
% The size of the array to split is obtained from the size of A.
% The dimension along which to split is obtained from the code. A captical
% code indicates that the array is a grid, that is, it contains cell face
% coordinates instead of cell values as of most fdm arrays.
% A second letter indicates that the array values must be
% divided according to the new splitted arrays. This is implied by the grid
% arrays.
% All arrays are first permuted such that the diretion to be split is the first
% dimension. This simplifies the code, as all directions are thus treated
% equally.
% First the splitArray is completed so that it contains entries for all
% cells.
% The second line gets ones for the cells that are not splitted.
% Then he array of the new cells in made: % NNew=1:sum(IOld)
% INew is a two-line array of size [2,NNew].
% The first line holds the cell indices of the new array.
% The seond line holds the corresponding indices of the old array.
% It is then trivial to generate the new array of the correct size
% and populated it with the values of the cells of the old array.
% This procedure should always work.
%
% To treat all grid arrays equally, we divide them into the coordinates of
% first layer (1,:,:) and the differences, i.e. the cell widths. This array
% is then treated like all others with subdivision implied.
% Once the new arrays are generated, we add the old coordinates to the
% front and apply a cumsum to regain the new coordinates.
% In a last step, all arrays are turned into their original directions.

