%% Verification of the solutions
% the different analytical solutions have been implemented using a class
% called solutionObj. The constructor requires as input a database
% consisting of an array of structs with prescribed fields that define the
% cases, one per line. One can of course also define a single case.
% The constructor only stores the case properties and does not solve the
% case. Solving the case, i.e. computing heads and flows on a grid along
% the cross section is done with the method "solve".
% This method is given the mame of the solution to be used, a matrix with
% three columns holding [t, N, E] and the grid coordinates along the xaxis,
% where L/2<x<L/2. The solution will then yield heads and flow at the end
% of all time steps.The initial head being constant, equal to the mean
% values in the ditches at either side of the cross section. Notice that
% only cases are dealt with where the water level in both ditches are the
% same.

S = solutionObj(P);

Sub = S(1:10);
names = solutions{1:3,:};

tic

%% Solve for all cases in Sub and solutions in names
for j = numel(Sub):-1:1
    for iname = numel(names):-1:1
        name = names{iname};
        S(j).(name) = S(j).solve(name,tne,xGr);
        S(j).(name).plot();
        S(j).(name).setGXG(tne);
    end
end

toc

% plot the GXG for all cases
%% plot GxG from model and analytic
leg=[];
figure; hold on; grid on;

plot([S.GLGDBF],'r--'); leg{end+1}='GLGDBF';
plot([S.GVGDBF],'g--'); leg{end+1}='GVGDBF';
plot([S.GHGDBF],'b--'); leg{end+1}='GHGDBF';

plot([S.GLG],'r');      leg{end+1}='GLG';
plot([S.GVG],'g');      leg{end+1}='GVG';
plot([S.GHG],'b');      leg{end+1}='GHG';

xlabel('Cross section number');
ylabel('head (above datum)');
grid on
title('GxG of the computed cross sections (GGOR tool, MODFLOW)');
legend(leg);
