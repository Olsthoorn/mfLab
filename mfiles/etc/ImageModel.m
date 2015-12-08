function [out1,out2]=ImageModel(figFileName,uMdl,vMdl,DX,DY)
%IMAGEMODEL Access to properties of an image relevant to its display.
%
% USAGE:
%   [uMdl,vMdl]=ImageModel(figFileName,DX,DY)
%
% Read the image and set the axes such that the pixel
% or ...
% allow getting the coordinates uMdl vMdl of the model on the immage.
%
% PROCEDURE:
%    rectangle uMdl vMdl has realworld size DX, DY
%    and the LL of the rectangle becomes x=0,y=0.
%
%    The function is used to plot contours on top of a photo of a sandbox
%    model. It may have other applications as well
% INPUTS:
%    figFileName = file name of image
%    uMdl = pixelcoordinates of rectangle to place correctly
%    vMdl = same, for y
%    DX size of rectangle (the model box) in world coordinates
%    DY size of rectangle [the model box)in world coordinates
% OUTPUT:
%    figure number with the axis holding the image with
%    the axis set as specified through the input.
%
% TO 100718

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin==3, DX=uMdl; DY=vMdl; uMdl=[]; vMdl=[]; end

A=imread(figFileName);  % read figure
A=A(end:-1:1,:,:);                    % turn upside down

if nargin<5
    figure; hold on; image(A); set(gca,'ydir','normal'); axis tight;
    fprintf('Click the LL and UR corner of the model in the picture\n');
    [uMdl,vMdl]=ginput(2);
    out1=uMdl;
    out2=vMdl;
end

%% Coordiantes of corners of model picked by ginput
u=[0.5 uMdl(1) uMdl(end) size(A,2)-0.5];
v=[0.5 vMdl(1) vMdl(end) size(A,1)-0.5];


%% Converting to x,y in cm

DU=u(3)-u(2); % [pixels]
DV=v(3)-v(2); % [pixels]

DXDU=DX/DU;
DYDV=DY/DV;

% LL corner of model in world coordinates
x(2)=0;
y(2)=0;

% LL corner of image in real world coordinates
x(1)=x(2)+DXDU*(u(1)-u(2));
y(1)=y(2)+DYDV*(v(1)-v(2));

% UR of model in world coordinates
x(3)=x(2)+DX;
y(3)=y(2)+DY;

% UR of image in world coordinates
x(4)=x(2)+DXDU*(u(4)-u(2));
y(4)=y(2)+DYDV*(v(4)-v(2));

xlim=[x(1) x(4)];
ylim=[y(1) y(4)];

%% Plot the figure in real world coordinates

figure;
image(xlim,ylim,A); set(gca,'ydir','normal');
set(gca,'xlimMode','manual','ylimMode','manual');

if nargin==5, out1=gcf; hold on; end

end
