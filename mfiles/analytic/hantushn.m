function drawdown=hantushnWork(Q,r,t,Sacc,S,c,T,N)
%HANTUSHN Analytical N-layer solution of flow in semiconf multiple aquifer system
% Applies Laplace transform with back-transformation from Laplace space
% according to Stehfest. Implies the solution given by Hemker and Maas
% (1987)
%
% Example:
%    hantushn();  % selfTest
%    s = hantushn(Q,r,t,Sacc,S,c,T[,N]);
%
% The solution has been used for the interpretation of the pumping test
% De Zilk, Pastoorslaan by Jochem Fritz and Rob Rijsen. Sept 3-14 2007.
%
% INPUTS:
%    Q = extraction vector (positive extractions will create positive dd')
%    r = distance vector
%    t = time vector
%
%    drawdown will be in [length(Aquif,length(r),length(t)] format
%    If instead r=r_iz_t list, the drawdown will be in list format
%    providing the drawdown for every coordinate pair
%    in that case t is dummy !
%
%    Sacc = Storage coefficient of aquitards
%    S = Storage coefficient of aquifers
%    c = aquitard resistance vector one more than aquifers
%    T = Transmissivity vector
%    N  = Stehfest's parameter, default 10
%
% The first resistance layer is on top of the first aquifer by definition.
% apply Inf to close off the given zero drawdown on top of this layer.
%    if length(c)and Sacc must be the same
%    if length(c)=length(T)+1, this means that there is a resistance
%    layer also at the bottom of the system and a fixed zero drawdown
%    below it.
%    If, however, length(T)==length(c), then the lowest layer is an
%    aquifer closed at its bottom.
%
% OUTPUT: dd(nLay,nR,nT)
%
% See also: stehfest
%
% TO 090329 090330

    fprintf('Running Hantushn\n');

    %% Selftest ??
    if nargin<7, drawdown=selfTest; return; end

    % Housekeeping
    Q=Q(:); t=t(:); Sacc=Sacc(:); c=c(:); S=S(:); T=T(:);

    if length(Sacc)~=length(c), error('length(Sacc)~=length(c)'); end
    if length(S   )~=length(T), error('length(S)~=length(T)');    end
    if length(Q   )~=length(T), error('length(Q)~=length(T)');    end
    
    % compute Stehfest's coeffient
    if ~exist('N','var'), N=10; end % Stehfest's default N
    v = Stehfest(N);

    %% Compute drawdown
    % T, r and t together form an array of values (nL,nR,nt)    
    drawdown=zeros(length(T),length(r),length(t));
    for ir=1:length(r)
        for it=1:length(t)
            drawdown(:,ir,it)=ddOnePoint(Q,r(ir),t(it),Sacc,S,c,T,v);
        end
    end
end

function dd=ddOnePoint(Q,r,t,Sacc,S,c,T,v)
    %DDONEPOINT -- compute dd for all layers for this r and t
    %
    
    s=zeros(size(T));
    for iStehfest=1:length(v)
        p=iStehfest*log(2)/t;
        d=p*S./T;              % number of aquifers
        b=sqrt(p*Sacc.*c);     % number of aquitard
        if length(T)==1
            eii=(b(1).*coth(b(1)))./(c(1).*T);
            eij=(b(2).*coth(b(2)))./(c(2).*T);
            A=eii+eij+d;
        else
            bcothb=b.*coth(b);
            bsinhb=b./sinh(b);
            if length(c)>length(T) % aquiard and zero drawdown at bottom
                eii=  bcothb(1:end-1)./(c(1:end-1).*T); % links to overlying  aquitard
                eij=  bcothb(2:end)  ./(c(2:end)  .*T); % links to underlying aquitard
                fii=  bsinhb(2:end-1)./(c(2:end-1).*T(2:end));
                fij=  bsinhb(2:end-1)./(c(2:end-1).*T(1:end-1));
            else % no bottom aquitard
                eii=  bcothb       ./(c.*T);
                eij= [bcothb(2:end)./(c(2:end).*T(1:end-1));0];
                fii=  bsinhb(2:end)./(c(2:end).*T(2:end));
                fij=  bsinhb(2:end)./(c(2:end).*T(1:end-1));
            end
            A=diag(eii+eij+d)-diag(fii,-1)-diag(fij,+1);
        end

        % this is suggested by Hemker in his papers, but seems to give
        % not perfectly the same results
        %         if false
        %            D = diag(T)^(1/2) * A * diag(T)^(-1/2);   
        %            [R,W]=eig(A);  %D
        %            V= diag(T)^(-1/2)*R;
        %         else
           [V,W]=eig(A);  %D
        %         end
        
        s=s + v(iStehfest)/(2*pi*p) * V * diag(besselk(0,r*sqrt(diag(W)))) * V^(-1) * (Q./T);
    end
    dd=s*log(2)/t;

end

function dd=selfTest()
    %SELFTEST -- fires when nargin==0, sovles problem stated in figure 2
    % of Hemker and Maas (1987)

    c    = [1000 1500 1000 4000 20000];   % resistance of aquitards
    Sacc = [   3     5   3    2 1]*1e-3; % Stor coef of aquitards
    T    = [2000  1500 500 2000];        % T of aquifers
    S    = [  10     4   1    3]*1e-4;   % Stor coef of aquifers
    Q    = [0, 10000, 0, 0];             % Extracion from aquifers

    r =logspace(1,log10(6000),40);
    t =[1e-5 1e-4 1e-3 1e-2 1e-1 1 10];

    dd=hantushn(Q,r,t,Sacc,S,c,T);

    close all

    %% =============fig 3 Hemker and Maas (1987) ===========
    defaults = {'nextplot','add','ydir','reverse','xGrid','on','yGrid','on'};
    
    figure('name','dd versus distance Hemker and Maas (1987) fig 2','pos',screenPos(0.8,0.5));
    
    ax(1) = subplot(1,2,1, defaults{:},'xscale','lin','yscale','lin','xlim',[0   6000],'ylim',[0      1]);
    ax(2) = subplot(1,2,2, defaults{:},'xscale','log','yscale','log','xlim',[10 10000],'ylim',[0.001 10]);
    
    xlabel(ax(1),'r [m]'); ylabel(ax(1),'dd [m]');
    xlabel(ax(2),'r [m]'); ylabel(ax(2),'dd [m]');
    title(ax(1), {'H/M 87, fig 3a, dd vs r';...
                ['t=[',sprintf('%s',sprintf(' %g',t)) ']']});
    title(ax(2), {'H/M 87, fig 3b, dd vs r';...
                ['t=[',sprintf('%s',sprintf(' %g',t)) ']']});
    
    for it=length(t):-1:1
        for iL=1:2 %size(dd,1)
            plot(ax(1),r,dd(iL,:,it),[mf_color(iL) mf_linetype(it)],'lineWidth',2);
            plot(ax(2),r,dd(iL,:,it),[mf_color(iL) mf_linetype(it)],'lineWidth',2);
        end
    end
end

function v=Stehfest(N)
%% Stehfest v-coefficients, need to be calculated only once
    v=NaN(N,1);
    for i = 1:N
        dum=0;
        for k=floor((i+1)/2):min(i,N/2)
            dum=dum+k^ (N/2)*factorial(2*k)/(factorial(N/2-k)*factorial(k)*factorial(k-1)*factorial(i-k)*factorial(2*k-i));
        end
        v(i)=(-1)^(i+N/2)*dum;
    end
end
