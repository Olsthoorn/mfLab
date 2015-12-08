function piezom=Pumptest(sheetname)
%PUMPTEST run Hantush or Theis pumping tests
%
% USAGE:
%   Pumptest('datafile','sheetname');
%
% Example:
%    1) grab the a data point or the matchpoint
%    2) shift them so that the points fit the left side of the type curves best.
%    3) Then compute the results
%       s*mp = W = s * 4*pi*T/Q --> mp = 4*pi*T/Q --> T=mp*Q/(4*pi)
%       t*mp = 1/ua = 4Tt/(r^2*Sa) --> Sa = 4T/(r^2*mp)
%
% See also: testanimator animator_loglog
%
% TO 120101

close all;

if ~exist('sheetname','var'), sheetname='Korendijk'; end

TypeCurvesHantush();

%% GetData

[vnams,vals]=getExcelData('../PumpingTests',sheetname,'Hor'); % Sheetname is 'Data'
Qd=vals(:,strmatchi('Q' ,vnams));
td=vals(:,strmatchi('td',vnams));      % column whose header starts with t (time)
rd=vals(:,strmatchi('rm',vnams));
zd=vals(:,strmatchi('zm',vnams));
sd=vals(:,strmatchi('ddnm',vnams)); % column whose header starts with WII_90/16

%% remove empty data
Qd=Qd(~isnan(sd));
td=td(~isnan(sd));
rd=rd(~isnan(sd));
zd=zd(~isnan(sd));
sd=sd(~isnan(sd));

UNR=unique([rd,zd],'rows');

clrs='brgkmcy';

piezom(size(UNR,1)).Q=Qd(1);
for i=1:size(UNR,1)
    j=rem(i,length(clrs));if j==0, j=length(clrs); end
    I=find(UNR(i,1)==rd & UNR(i,2)==zd);
    piezom(i).Q=Qd(1);
    piezom(i).rd=UNR(i,1);
    piezom(i).zd=UNR(i,2);
    piezom(i).td=td(I);
    piezom(i).sd=sd(I);
    piezom(i).tr2=piezom(i).td./piezom(i).rd.^2;
    plot(piezom(i).tr2,piezom(i).sd/piezom(i).Q,[clrs(j) 'o'],...
        'tag','myData',...
        'ButtonDownFcn','animator_loglog start',...
        'EraseMode','xor');
end


%% Steady and moving matchpoints

xM=1; yM=1; % match point

% steady
plot(xM,yM,'ro',...
    'UserData',[xM yM],...
    'MarkerFaceColor','r',...
    'tag','MatchPoint',...
    'ButtonDownFcn',...
    'animator_loglog start',...
    'EraseMode','xor');

% moving
plot(xM,yM,'ko',...
    'DisplayName','MatchPoint',...
    'UserData',[xM yM],...
    'MarkerFaceColor','k',...
    'tag','myData',...
    'ButtonDownFcn',...
    'animator_loglog start',...
    'EraseMode','xor');

get(findobj('tag','MatchPoint'),'UserData');

