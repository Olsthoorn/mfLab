function Vennebulten(varargin)
% Vennebulten % Interactive Workout of delayed yield pummping test Vennebulten
% (Kruseman & De Ridder, 1970).
% USAGE:
% Vennebulten
%
% Boulton type curves will be computed and the data read and plotted,
% together with an inial matchpoint at (1,1).
%
% First grab the a data point or the matchpoint and shift them so that
% the points fit the left side of the type curves best. This is the early
% time drawdown. Then proceed to match the late time. This is done by
% grabing the late Theiscurve and moving it such that the data points best
% fit the entire type curves.
% The compute the results
%  s*mp = W = s * 4*pi*T/Q --> mp = 4*pi*T/Q --> T=mp*Q/(4*pi)
%  t*mp = 1/ua = 4Tt/(r^2*Sa) --> Sa = 4T/(r^2*mp)
% Then get Sy, this equals Sa times the distance between the two theis
% curves which can be read on the graph.
% Finally obtain Boulton's leakage factor. It can be read from the matching
% graph, which has its r/B, with B=sqrt(T/(alpha*Sa)) --> alpha=T/Sa/B^2.
% tw=1/alpha.

T=1; Sa =4*T; Q=4*pi*T; Alpha=0.5;

if nargin==0,
    action='init';
    close all;
    Sy=Sa*1e2;
elseif isnumeric(varargin{1})
    action='update';
    Sy=varargin{1};
elseif ischar(varargin{1})
    action=varargin{1};
    Sy=Sa*1e2;
else
    error('illegal argument <<%s>> in Vennebutten',varargin{1});
end

c=Sy/Alpha; lambda=sqrt(T*c);


switch action
    case 'start'
        set(gcbf,'WindowButtonMotionFcn',[mfilename ' move']);
        set(gcbf,'WindowButtonUpFcn',    [mfilename ' stop']);
        set(gcbf,'WindowButtonDownFcn',  '');
        hy=findobj('tag','Sy');
        set(hy,'UserData',get(gca,'CurrentPoint'));
        return;
    case 'stop'
        set(gcbf,'WindowButtonMotionFcn','');
        set(gcbf,'WindowButtonUpFcn','');
        return;
    case 'move'
        Sy=get(gcf,'UserData');
end

%% Grid

% Make sure we will always have a suitable scale of the axis
r1=1e-6*lambda; r2=1e4 *lambda;

% with sufficient puoints for accuracy
r=logspace(log10(r1),log10(r2),10*(log10(r2)-log10(r1)));

b=1; % use unit layer thickness for convenience, to that k=T
y=[-b 0 1];

% grid cleanup and housekeeping
[r,y,rm,~,~,B,Nr,Ny]=modelsize(r,y); % B=dy=layer thickness

% Make sure we always have a suitable time axis
t1=1e-5;
t2=1e1;
t=logspace(log10(t1/r2^2),log10(t2/r1^2),10*(log10(t2/r1^2)-log10(t1/r2^2)))';

%% The model
kh=[0  ;T/B(2)];
kv=[B(1)/c; Inf];
S =[Sy/B(1); Sa/B(2)];

%% Numerical boundary conditions
IBOUND=ones(Ny,Nr);
IH=zeros(Ny,Nr);
FQ=zeros(Ny,Nr);
FQ(end,1)=Q;

rnge1=floor(log10(min(r)^2/(4*T*max(t))));
rnge2=ceil (log10(max(r)^2/(4*T*min(t))));
rnge =logspace(rnge1,rnge2,10*(rnge2-rnge1));

cutoff = 3; %*gamma;

ua=Sa.*rnge; ua=ua(ua<cutoff);
uy=Sa.*rnge; uy=uy(uy<cutoff);
uB=Sa./(4*T*t)*rm.^2; uB(uB>cutoff)=NaN;

phia=Q/(4*pi*T)*expint(ua);
phiy=Q/(4*pi*T)*expint(uy);

IH(end,:)=Q/(4*pi*T)*expint(Sa./(4*T*t(1))*rm.^2);

phiB=fdm2t(r,y,t,kh,kv,S,IBOUND,IH,FQ,'axial'); phiB=permute(phiB,[3,2,1]);

%% Plot it
I=find(t>=1e-6); % remove first times
q=Q/(4*pi*T);

BB=sqrt(rm/lambda);
iI=round(0.5*length(I));

switch action
    case 'init'
        % setup the plot
        figure; hold on; grid on; fontsize=15;
        xlabel('1/u_a = t/r^2 * (4T/S_y)','fontsize',fontsize);
        ylabel('s/(4\pi T)','fontsize',fontsize); 
        title('Boulton type curves','fontsize',fontsize);
        ylim=[1e-2 1e2]; xlim=[1e-5 1e4];
        set(gca,'xscale','log','yscale','log','ylim',ylim,'xlim',xlim,'fontsize',fontsize);

        set(gcf,'UserData',Sy);
        
        for j=1:size(uB,2)
            plot(1./uB(I,j), phiB(I,j,end)/q,'b', 'tag','SB','erase','xor');  % numeric  Boulton
        end
        plot(    1./ua , phia/q,'r', 'tag','Sa','erase','xor','linewidth',1);  % analytic Sa
        plot(Sy/Sa./uy , phiy/q,'k', 'tag','Sy','erase','xor','linewidth',1,...
            'ButtonDownFcn','Vennebulten start');  % analytic Sy

        %% Put r/B labels
        for j=1:size(uB,2)
            text(1./uB(iI,j)',phiB(iI,j,end)'/q,sprintf('%4g',BB(j)),'clipping','on','tag','text');
        end

        %% GetData

        [vnams,vals]=getExcelData(mfilename,'Data','Hor'); % Sheetname is 'Data'

        td=vals(:,strmatchi('time',vnams));      % column whose header starts with t (time)
        sd=vals(:,strmatchi('WII_90/16',vnams)); % column whose header starts with WII_90/16

        %% remove empty data
        td=td(~isnan(sd))/(24*60);  % minutes to days
        sd=sd(~isnan(sd));          % remove empty data lines

        plot(td,sd,'ro',1,1,'ko','tag','myData','ButtonDownFcn','animator start','EraseMode','xor');
    case 'move'
        hy=findobj('tag','Sy');
        xyP=get(gca,'CurrentPoint');
        xy0=get(hy,'UserData');
        
        % Update Sy and save for next round
        Sy=Sy*xy0(1,2)/xyP(1,2); set(gcf,'UserData',Sy);
        
        set(hy,'xdata',Sy/Sa./uy,'ydata',phiy/q);
        set(hy,'UserData',xyP);
        
        hB=findobj('tag','SB');
        for i=1:length(hB)
            set(hB(i),'xData',1./uB(I,i),'yData',phiB(I,i,end)/q);
        end
        
        hT=findobj('tag','text');
        for j=1:length(hT)
            set(hT(j),'Position',[1/uB(iI,j),phiB(iI,j,end)/q,0],'string',sprintf('%.2f',BB(j)));
        end
end




