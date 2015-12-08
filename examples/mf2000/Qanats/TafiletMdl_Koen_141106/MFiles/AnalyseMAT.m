%% Analyse results .MAT files
% Plottable(o,a,b)
%   where o is the table with the information
%   a indicates whether plotting khettaras
%   b is whether you want to save it (0 in this case)
%
% 

display(o);
Length = strmatchi('Length',o.Properties.VariableNames);

if any(Length)
    Plottable(o,1,0)
else
    Plottable(o,2,0)
end