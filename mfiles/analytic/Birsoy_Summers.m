function s = Birsoy_Summers(Q,tStart,kD,S,r,time)
%BIRSOY_SUMMERS Variable-discharge well tests and tests in well fields
%   (See Kruseman & De Ridder, 1994, Section 12.1)
%    Birsoy YK & Summers WK (1980) Determination of aquifer parameters from
%    step tests and intermittent pumping data. Ground Water 18, p 137-146.
%
% Example:
%     Birsoy_Summeres(Q,t,kD,S,r,time)
%
% INPUT
%     Q vector of extractioins during the n periods
%     t vector of times, one more than Q marking start and ends of periods
%     kD = transmissivity
%     S  = storage coefficient
%     r  = distance to well
%     time = time for desired output
%
% OUTPUT
%     s = drawdown at times time.
%
%   The drawdown in a step-drawdown test can be converted by this method
%   to that of a continuous test by adapting the time.
%
%   TODO: thoroughly test and document
%
%   TO 010101

% Pumping Flow and Periods

if nargin==0; selfTest(); return; end

Q=Q(:);
tStart=tStart(:);

dQ   = [Q(1); diff(Q)];
Nper = numel(dQ);

s = zeros(size(time));

for iPer=1:Nper
    It = time>tStart(iPer);
    s(It) = s(It) + dQ(iPer)/(4*pi*kD) * log(2.25*kD/(r^2*S) *beta(iPer));
end

    function b= beta(iPer)
        b = ones(size(time));
        for i = 1:iPer
            b(It) = b(It) .* (time(It)-tStart(i)).^(dQ(i)/Q(iPer));
        end
        b = b(It);
    end
end

function s = theis(Q,tStart,kD,S,r,time)
    s = zeros(size(time));
    dQ = [Q(1); diff(Q)];
    
    for iP = 1:numel(Q)
        It = time>tStart(iP);
        s(It) = s(It)+dQ(iP)/(4*pi*kD)*expint(r.^2*S./(4*kD*(time(It)-tStart(iP))));
    end
end

function selfTest()

% Data from Kruseman & De Ridder (1994), p185
Data = [ % n, t(min), Sn, Qn
    1   5 1.38 500
    1  10 1.65 500
    1  15 1.81 500
    1  20 1.93 500
    1  25 2.02 500
    1  30 2.09 500
    2  35 2.68 700
    2  40 2.85 700
    2  45 2.96 700
    2  50 3.05 700
    2  55 3.12 700
    2  60 3.18 700
    2  70 3.29 700
    3  80 3.38 600
    3  90 3.13 600
    3 110 3.15 600
    3 110 3.17 600
    3 130 3.23 600];

%% other data

kD=600;
S =0.1;
r =0.2;

tStart = Data(:,2) / (24*60);
sMeas  = Data(:,3);
Q      = Data(:,4);

time = logspace(-3,0,41);

sBS = Birsoy_Summers(Q,tStart,kD,S,r,time);

% Alternatively do plain superposition with Theis

sTh = theis(Q,tStart,kD,S,r,time);

figure; xlabel('t [d]'); ylabel('s [m]'); title('Birsoy-Summers/Theis');
set(gca,'nextplot','add','xScale','log');

plot(tStart,sMeas,'ko');
plot(time,Q/(4*pi*kD)*sBS,'b');
plot(time,Q/(4*pi*kD)*sTh,'r');
legend('measured','Birsoy','Theis');

end
