%Bruggeman script
% Solves the multilayer solutions in Bruggeman (1999) Analytical Solutions
% of Geohydrological Problems, Elsevier, ISBN 0-444-81829-4
% Cases 710-720, pages 433-449.
% Uses classdef BruggemanMultilayerObj in mflab>mfiles>Analytic
%
% This is solved in an object oriented way.
%
% TO 110803
% 

clear; close all;

%% Showing all cases subsequentely
%  maincase is 710 (one-dimensional) or 720 (axi-symmetric)
% all subcases have been numbered by Bruggeman as follows

subcase710=[1 2    12 13 15 16 17 18 21 22 23 24 25 26 27 28];
subcase720=[1 3 11 12 13 15 16 17    21 22 23 24             31 32 33 34 35];

NLay=5;  % may be any number but less or equal to number of layers defined

%% One dimensional cases

% Precallocate struct with right number of objects using repmat
Brug710=repmat(BruggemanMultilayerObj(NLay),length(subcase710),1);

% Note that every object has the same geohydrological data but will be
% solved differently, according to the case (different boundary conditions)

% then solve and plot each case
for i=1:length(subcase710)
    thiscase=710+subcase710(i)/100;
    Brug710(i)=Brug710(i).solve(thiscase).show;           % solve and plot
end

%% Axially symmetric cases

% Precallocate struct with right number of objects using repmat
Brug720=repmat(BruggemanMultilayerObj(NLay),length(subcase720),1);

% then solve and plot each case
for i=1:length(subcase720)
    thiscase=720+subcase720(i)/100;
    Brug720(i)=Brug720(i).solve(thiscase).show;
end
