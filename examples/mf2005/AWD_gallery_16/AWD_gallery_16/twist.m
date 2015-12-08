function B=twist(B)
%   TWIST  permute 3d array to exchange x and y (repeat to rotage back)
%
%   USAGE B=test(B)
%
%  B must be a 3D array
%
%   TO 110813
    B=permute(B(:,end:-1:1,:),[2,1,3]);
    
%    a=reshape(1:2*3*4,[2,3,4]);