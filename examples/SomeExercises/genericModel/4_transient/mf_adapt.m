%% Generic transient model

%% Explanation
% This is the same model as 1_Steady, with only a minor change to make it
% transient. The change to be made is in setting the value(s) in the column
% transient in the worksheet PER of the accompanying workbook to 1 instead
% of 0. mf_lay2mdl will then automatically retrieve SS and SY from the
% workbook LAY and geneate the 3D model arrays. The model can then be run
% from the initial STRTHD value to see its evaluation over time.
% The visualization is adapted accordingly in mf_analyze.
% Because the steps have been explained in detail in 1_steady, we only give
% the required steps without comments here.
% The obtain a somewhat nicer looking output picture, the well extraction was
% reduced by a factor 5 to -1e5 m3/d (see worksheet PER column with header Q)
% and the cDrn was increased from 1 to 100 m/d (see below, line 36)

%% Steps

basename = 'generic_transient';  % set basename
save('name','basename');         % save it for later retrieval

xGr = 0:5000:75000;  % grid line x-coordinates
yGr = 75000:-5000:0; % gird line y-coordinates

mf_lay2mdl(basename,xGr,yGr); % generate the 3D arrays and gridObj from worksheet LAY

gr = gridObj(xGr,yGr,gr.zGr+200,gr.LAYCBD); % change default top of gr.zGr to +200

well = wellObj(basename,'wells',gr,HK,{'PER','Q'}); % get wells and their flows

IBOUND(:,1,[1 2]) = -1;  % specify fixed-head cells

line = [ 7500 37500   0 ;
        47500 37500 100];  % the position of the drain (polyline)

cDrn = 100.0; % [m/d]

DRN = gr.bcnLine(basename,'DRN',line,cDrn); % generate DRN input

save underneath line
%% Conclusion
% In fact, nothing changes with respect to the 1_steady generic model. The only
% difference is the value in the column transient in worksheet PER.

%% TO 130614
