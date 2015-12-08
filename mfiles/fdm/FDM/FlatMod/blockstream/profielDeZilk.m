% profiel DeZilk, (Stuyfzand, EE')
% script run from runmodel specifying the particular profile

OldPath=path;
p1='z:\projects\Hy171 Zoet zout modellering\Psiphimodel\mfiles';
p2='z:\models\matlabmodels\blockstream';
path(path,p1);
path(path,p2);

ProfileName='EE'

%if ~varexists('xSEE',who); load coordsee; end; 			% brackish and salt interface coordinates

eval(['cd ''',p1,'''']);
if ~varexists(ProfileName,who); rdprofiles; end 				% dune elevation profiles
eval(['cd ''',p2,'''']);

%if ~varexists('xSEE',who); load coordsee; end; 			% brackish and salt interface coordinates

%eval(['Duin=',ProfileName,';']);

% The hydrological system; can immediately be computed analytically using nsecn (see DsnEEStuyf88.m)
%load coodir\coee1981
situation=1998;
switch situation
case 1925,
QZ =[ [     0, -0.4,  0,   0,-1.3,   0,   0,   0,   0];...
      [     0,    0,  0,   0,-1.585,   0,  -1,   0,   0];...
      [     0,    0,  0,   0,   0,   0,   0,   0,   0]];					% zone boundary extractions
%load coordsee1981
case 1957,
QZ =[ [     0,-0.5,  0,   0, -1.57,   0,   0,   0,   0];...
      [     0,   0,  0,   0, -0.0,   0,-1.3,   0,   0];...
      [     0,   0,  0,   0,    0,   0,   0,   0,   0]];					% zone boundary extractions
%load coordsee1981
case 1982,
QZ =[ [     0,-0.5,  0,   0,-1.60,   0,   0,   0,   0];...
      [     0,   0,  0,   0,-0.0,   0,   0,   0,   0];...
      [     0,   0,  0,   0,   0,   0,   0,   0,   0]];					% zone boundary extractions
case 1984,
QZ =[ [     0,-0.5,  0,   0,-1.60,   0,   0,   0,   0];...
      [     0,   0,  0,   0,-0.4,   0,   0,   0,   0];...
      [     0,   0,  0,   0,   0,   0,   0,   0,   0]];					% zone boundary extractions
case 1995,  %VLS gedempt, peil OK -1.25
QZ =[ [     0,   0,  0,   0,-1.64,   0,   0,   0,   0];...
      [     0,   0,  0,   0, -0.4,   0,   0,   0,   0];...
      [     0,   0,  0,   0,    0,   0,   0,   0,   0]];					% zone boundary extractions
case 1996,  %peil OK -0.75
QZ =[ [     0,   0,  0,   0,-1.338,   0,   0,   0,   0];...
      [     0,   0,  0,   0, -0.4,   0,   0,   0,   0];...
      [     0,   0,  0,   0,    0,   0,   0,   0,   0]];					% zone boundary extractions
case 1998,  %Peil OK -0.75, winning OK uit
QZ =[ [     0,   0,  0,   0,-1.42,   0,   0,   0,   0];...
      [     0,   0,  0,   0,    0,   0,   0,   0,   0];...
      [     0,   0,  0,   0,    0,   0,   0,   0,   0]];					% zone boundary extractions
case 2000,  % winning OK uit
QZ =[ [     0,   0,  0,   0,    0,   0,   0,   0,   0];...
      [     0,   0,  0,   0,    0,   0,   0,   0,   0];...
      [     0,   0,  0,   0,    0,   0,   0,   0,   0]];					% zone boundary extractions
end
x=    [     0,  750,1750,2500,3550,4225,5500,6750,8000];					% zone boundaries
h=    [   0,  0, 5.0,  6.5,   3, -0.6, -0.6, -0.6, -1.8, -6];			% prescribed heads
%h=    [   0,  0, 5.0,  6.5,   3,   -0,   -0, -0.0,  -0,  -0];			% prescribed heads
prev=1.1*[1e-3, 1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];				% precipitation-evapotranspiration
c=[   [  10, 0.33,0.25,0.25,0.33,0.50, 100, 100, 250, 500];...
      [ 500,  750,1500,1750,1550,1250,1250,1250,1250,1500];...
      [ 150,  150, 150, 150, 250, 250, 250, 250, 500, 750]];			% resistance of topsystem
kD=[  [ 150,  100, 100, 100, 100, 100, 100, 100, 100, 100];...
      [1500, 1500,1500,1500,1500,1500,1500,1500,1500,1500];...
      [2250, 2250,2250,2250,2250,2250,2250,2250,2250,2250]];			% aquifer system
isfh= [  1,    0,   0,   0,   0,   0,   1,   1,   1,   1];			% h is prescribed on top or not (=prev);
D      =[5;10;5;40;10;60];					% Thickness of analytic aquitard-aquifer sequence
por=0.35;
kLayers=[D(1)./c(1,:); kD(1,:)/D(2); D(3)./c(2,:); kD(2,:)/D(4); D(5)./c(3,:); kD(3,:)./D(6)];
ddeltaB=5/18*0.025
ddeltaS=16/18*0.025-ddeltaB

