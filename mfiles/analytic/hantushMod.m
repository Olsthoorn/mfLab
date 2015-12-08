function [s,tp]=hantushMod(Q,k,Ss,r,z,t,screen,D)
    %HANTUSHMOD computes hantush but for partially penetrating well + self-test
    %
    %   Ref: Kruseman & De Ridder (1994), p168++
    %      There is a solution of early times and one for late times.
    %      This is the solution for late times
    %      The essence is that the partial penetration effect for shot times is transient.
    %      This function will compute the Hantush function for all times.
    %
    %   The short term solution is that of a screen in a half inifnite
    %   space. It can be compared to a point source in half infinite space
    %   See bruggeman (1999) solution 410
    %
    % USAGE:
    %    hantushMod();   % selfTest
    %    [s,tp]=hantushMod(Q,k,Ss,r,t,D);                      % Theis
    %    [s,tp]=hantushMod(Q,k,Ss,r,z,t,[zScrBot,zScrTop]);    % Hantush mod valid t<tp
    %    [s,tp]=hantushMod(Q,k,Ss,r,z,t,[zScrBot,zScrTop],D);  % Hantush mod valid t>tp
    %    [s,tp]=hantushMod(Q,k,Ss,r,z,t, z0);                  % Bruggeman half inf space
    %    [s,tp]=hantushMod(Q,k,Ss,r,z,t, z0,D);                % Bruggeman half aquif D given
    %
    %    the transition time for short and long term validity coincides with the
    %    time the drawdown from the bottom of the aquifer reflects
    %       t=tp = (Ss/k)/20*(2D-(b+z))^2
    %    for z=b=0
    %       tp = D^2*Ss/(5k) that is
    %
    %   To deal with anisotropy use [kh kv] instead of kh
    %
    % See also: hantush hantushn
    %
    % TO 090101

    % selfTest if nargin==0
    if nargin==0, selftest(); return; end

    %% Assert inputs
    if nargin<6
        error('needs at least 6 and at most 8 inputs');
    end
    
    if nargin==6
        D = t;
        t = z;
        screen = [0 D]';
        mode = 'theis';
    elseif numel(screen)==1
        if nargin==7
            mode= 'bruggeman410Inf';
        else
            mode= 'bruggeman410';
        end
    else
        if nargin==7
            mode = 'hantushShort';
        else
            mode = 'hantushLong';
        end
    end
    
    %% All distances are taken relative to top of aquifer
    Q        = abs(Q);
    
    %% anisotropy
    if size(k,2)==2
        kh = k(1); kv=k(2);
        k = (kh*kh*kv).^(1/3);
        r = sqrt(k/kh).*r;
        z = sqrt(k/kv).*z;
        screen = sqrt(k/kv)*screen;
        if exist('D','var')
            D = sqrt(k/kv) * D;
        end
    end
    
    screen  = abs(screen);
    zScrBot = min(screen);
    zScrTop = max(screen);
    z       = abs(z);

    
    switch mode
        case 'theis'
            %%THEIS -- Q/(4*pi*k*D) expint(u), u=r^2*Ss/(4*k*D);
            r = r(:); t=t(:)';
            u = bsxfun(@times,r.^2,Ss./(4*k*t));
            s = Q/(4*pi*k*D) * expint(u);
            tp=0;  % transection point between short and long term validity
        case 'hantushShort'
            r = r(:); t=t(:)';
            u = bsxfun(@times,r.^2,Ss./(4*k*t));
            B = [(zScrTop+z)./r ...
                 (zScrBot+z)./r ...
                 (zScrTop-z)./r ...
                 (zScrBot-z)./r];

            % B is always of size numel(r),4)
            % u is always of size numel(r),numel(t)
            w = M(u,B(:,1))-M(u,B(:,2))+M(u,B(:,3))-M(u,B(:,4));
            s= Q/(8*pi*k*(zScrTop-zScrBot)) * w;
            tp = Inf;

        case 'hantushLong'
            if nargin<7 || ~isnumeric(D)
                error('6th argument must be the aquifer depth');
            end
            r = r(:); t=t(:)';
            u = r.^2*(Ss./(4*k*t));
            w = bsxfun(@plus,expint(u),hantushPP(r,zScrTop,zScrBot,z,D));
            s= Q/(4 * pi * k * D) * w;
            tp = Ss/(20*k)*(2*D-(zScrTop+z))^2;
        case 'bruggeman410Inf'
            s = brug410(Q,k,Ss,r,z,t,zScrTop);
            tp = Inf;
        case 'bruggeman410'
            s = brug410(Q,k,Ss,r,z,t,zScrTop,D);
            tp = Ss/(20*k)*(2*D-(zScrTop+z))^2;
        otherwise
            error('unknown mode see help of this function');
    end

end

function w=M(u,B)
%M computes Hantush's modification of the Theis method for partially penetrating wells
% Valdid for relative short times (ds_pp)
%
% Example:
%    w=M(u,r/Lambda)
%
%   Ref: Kruseman & De Ridder (1994), p162++
%      There is a solution of early times and one for late times.
%      The essence here is that the partial penetration effect is transient.
%
%   The function is embedded in hantushE
%
%   M(u,b/r,d/r,a/r) = M(u,B1)-M(u,B2)+M(u,B3)-M(u,B4)
%
% with:
%   u = r^2S/(4kt)
%   Ss = S/D;
%   B1 = (b+a)/r
%   B2 = (d+a)/r
%   B3 = (b-a)/r
%   B4 = (d-a)/r
%   M(u,B) = int(exp(y)/a*erf(B sqrt(y)), u, Inf)
%
% See also: hantushE hantush hantushn
%
% TO 010409 121001

    if any(u)<=0;
        error('u must be >0');
    end

    B=B(:); % B is always [numel(r,1)]

    NT = size(u,2);  % times
    NP = size(B,1);  % points

    w = NaN(NP,NT);

    for iu=1:NT  % this is: for all times
        y=bsxfun(@times,u(:,iu),logspace(0,13,5000));
        dy = diff(y,1,2);
        arg = exp(-y)./y.*erf(bsxfun(@times,B,sqrt(y)));
        w(:,iu)= sum((arg(:,1:end-1) + arg(:,2:end)).*dy/2,2);
    end

end

function f = hantushPP(r,zScrTop,zScrBot,z,D)
    %%HANTUSHPP -- computes partial penetration for relative long times
    %
    % See Kruseman and De Ridder (1994), p 168
    % USAGE:
    %   w = hantushPP(r,zSscrTop,zSscrBot,z,D)
    %
    % TO 130630
    
    f = 0;
    r = pi*r/D;
    b = pi*zScrBot/D;
    d = pi*zScrTop/D;
    a = pi*z/D;
    N=10000;
    for i=1:N
        df = (1/i)*besselk(0,i*r).*cos(i*a).*(sin(i*b)-sin(i*d));
        f = f + df;
        if all(abs(df)<1e-6), break; end
    end
    if i==N
        error('No convergence after %d terms',N);
    end
    f = f * 4 ./(b-d);
end

function s = brug410(Q,k,Ss,r,z,t,z0,D)
% works for many r, many t but 1 z at a time
    Q = abs(Q); z=abs(z); z0=abs(z0);
    if nargin>7, D=abs(D); end
    
    r=r(:); t=t(:)';
    r2 = r.^2;
    uu = sqrt(Ss./(4*k*t));
    
    %% First is for infinitely deep aquifer
    rho1 = sqrt(r2+(z-z0)^2);
    rho2 = sqrt(r2+(z+z0)^2);
    s  =   bsxfun(@times,Q./(4*pi*k*rho1),erfc(bsxfun(@times,rho1,uu))) + ...
         + bsxfun(@times,Q./(4*pi*k*rho2),erfc(bsxfun(@times,rho2,uu)));
     if nargin<8
         return;
     end
     
    % superposition
    N=1000;
    for i=1:N
        rho1p = sqrt(r2+((2*i*D-z0)-z)^2);
        rho2p = sqrt(r2+((2*i*D+z0)-z)^2);
        rho1m = sqrt(r2+((2*i*D-z0)+z)^2);
        rho2m = sqrt(r2+((2*i*D+z0)+z)^2);
        ds  =   bsxfun(@times,Q./(4*pi*k*rho1p),erfc(bsxfun(@times,rho1p,uu))) + ...
              + bsxfun(@times,Q./(4*pi*k*rho2p),erfc(bsxfun(@times,rho2p,uu))) + ...
              + bsxfun(@times,Q./(4*pi*k*rho1m),erfc(bsxfun(@times,rho1m,uu))) + ...
              + bsxfun(@times,Q./(4*pi*k*rho2m),erfc(bsxfun(@times,rho2m,uu)));
        s = s + ds;
        if max(abs(ds))<1e-6
            break;
        end
    end
end

function selftest()
%SELFTEST selftest for Hantushs modification for partially penetrating wells
%
% Example:
%   hantushE();
%
% TO 130429

Q = 2400;
k = 10;
Ss = 1e-5;
t = logspace(-3,2,51);
screen = [ -50,9,-60.1];
z = 0;
z0= -55;
r = logspace(1,3,6);
D = 300;

[s0,tp0]=hantushMod(Q,k,Ss,r,t,D);           % theis
[s1,tp1]=hantushMod(Q,k,Ss,r,z,t,screen);    % Hantush modified, valid t<tp
[s2,tp2]=hantushMod(Q,k,Ss,r,z,t,screen,D);  % Hantush modified, valid t>tp
[s3,tp3]=hantushMod(Q,k,Ss,r,z,t,z0);        % Bruggeman half inf space
[s4,tp4]=hantushMod(Q,k,Ss,r,z,t,z0,D);      % Bruggeman half aquif D given

figure;
defaults = {'xGrid','on','yGrid','on','nextPlot','add','xlim',t([1 end]),'xScale','log'};

leg={'theis','hantushShort','hantushLong','3DpoinInf','3DpointD'};

ax = axes(defaults{:});

xlabel(ax,'t [d]'); ylabel(ax,'s [m]'); hold on
title(ax,sprintf('%s test hantush variants, partially penetrating wells',mfilename));

plot(t,s0);
plot(t,s1);
plot(t,s2,'o');
%plot(t,s3,'s');
%plot(t,s4,'x');

fprintf('tp different situations is: ');
fprintf('theis=%.4g hShort=%.4g, hLong=%.4g, brug3DInf=%.4g brug3D=%.4g\n',...
    [tp0 tp1 tp2 tp3 tp4]);
fprintf('\n');

legend(leg,4);
plot([tp2 tp2],get(ax,'ylim'),'k');
plot([tp4 tp4],get(ax,'ylim'),'r');

end