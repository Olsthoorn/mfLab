function F = CJp335_6(varargin)
%CJP335_6 -- computes formula Carslaw & Jaegr, p335, 6, temp around constant temp rod
% The solution for the temperature around a rod with circular cross section
% of radius Rp in infinite space initially at zero temp. The temperature in
% the rod is swichted to unity at t=0 and stays so afterwards.
% The solution is given as an integral between 0 and infinitity.
%
% Computing the integral is not trival because the integrand
% oscillates and has a singularity at u=0. The procedure to compute the
% integral partly borrowed from Peng, Yeh and Yang (2002) Advances in Water
% Resources, Vol 25, 663-675.
%
% The procedure is as follows:
% 1) estimate the first root. Doesn't have to be exact.
% 2) estimate a value of u << first root to cover the intergral over the
% singularity analytically.
% 3) estimate the upper bound of the integration using the exponal damping,
%    making sure the exp(-tau u^2)<UTOL
% 1) Integrate first part 0->u1 uing analytical fomrmula form Peng et al.
% 2) Integrate second part u1->u2 by quad with given tolerance
% 3) Integrate third  part u2->u3 by quad with given tolerance
% add the three contributions:
% Compute F=1+2/pi sum(q);
%
% The results are as accurate as the table by Peng et al. (2002).
%
% USAGE: F = CJp335_6(rho,tau,show)
%    rho = r/Rp
%    tau = kappa/R^2 *t
%      and kappa=lambda/rhoc, the diffusivity of the medium around the rod.
%    show = string 'show' or omitted, triggers showing itself
%
%  The necessary subfunctions are included in this m-file
%
% TO 140131

    if nargin==0
        r      = 2.0;  % m
        Rp     = 0.6;  % m
        rho    = r/Rp; % [-]
        t      = 3600; % s
        kappa  = 7.73e-7; % m2/s
        tau    = t * kappa/Rp.*2;
        F      = CJp335_6(rho,tau,'show');
        return;
    else
        [rho, varargin] = getNext(varargin,'double',[]);
        [tau, varargin] = getNext(varargin,'double',[]);
        [show, ~]       = getWord(varargin,'show');
        if isempty(tau)
            error('need at least two input arguments, rho and tau');
        end
    end
    
    %% Get the roots
    UTOL  = 1e-30;

    u1   = pi/(rho-1);  % estimate of first root
    uMax = sqrt(-log(UTOL)/tau);  % end of integration when damping < UTOL
    uRt  = [u1/1000 u1  uMax];
    if uRt(end)<uRt(1)
        uRt = uRt(end) ./ [3 2 1];
    elseif uRt(end)<uRt(2)
        uRt(2)=uRt(end)/2;
    end
    q    = zeros(size(uRt));    
    
    % Integrate from 0 to uRt(1) across singularity
    q(1) = log(rho)*atan(besselj(0,uRt(1))/bessely(0,uRt(1))); % singularity

    % For the two remaining intervals use quad with given accuracy
    TOL   = 1e-8;               % integration tolerance
    TRACE = true;               % show quad convergence
    func = @(u) Fu(u,rho,tau);  % define anonymous function   
    for i=numel(uRt):-1:2
        q(i) = quad(func,uRt(i-1),uRt(i),TOL,~TRACE);
    end
    
    % Assemble fhe function
    F = max(0,1+2/pi*sum(q));

%    fprintf('u = %10g %10g %10g\n',uRt);
%    fprintf('q = %10g %10g %10g\n',q);
%    fprintf('F = %10g\n',F);
    

    % If desired show the function
    if show
        fsz = 14;
        figure('name','integrand+roots','pos',screenPos(0.6));
        axes('nextPlot','add','fontsize',fsz);
        xlabel('u','fontsize',fsz);
        ylabel('integrand','fontsize',fsz);
        title(sprintf('Intgrand of Carslaw and Jaeger(1959) p335 eq6 for \\rho= %g, \\tau=%g',rho,tau));
        u = linspace(uRt(1),uRt(end),5000);  % entire span of uRt
        plot(u,-func(u));
        plot(uRt,zeros(size(uRt)),'ro');
    end

%    fprintf('rho = %10g tau = %10g nRoots = %3d F = %10.5f\n',rho,tau,numel(uRt),F);
%    fprintf('\n');
end

function uNew = newtonRoot_u(u,rho,tau)
%NEWTONROOT -- find next root of Fu
% USAGE uNew = newtonRoot_u(u,rho,tau)
%      rho = r/Rp;
%      tau = t * kappa / Rp^2      (kappa = lambda/(rhoc) = diffusivity)
%
% SEE ALSO: fu, facc
%
% TO 140131

    TOL = 1e-10;

    MAXITER = 50;
    for i=1:MAXITER
        [fa,fu_fa] = Facc(u,rho,tau);
        if isnan(fa),
            uNew = u;
            return;
        end
        uNew = u - fu_fa;
        if abs(uNew-u)<TOL || isnan(uNew)
            return;
        else
            u=uNew;
        end
    end
    uNew=NaN; % No convergence
end

function [fa,fu_fa] = Facc(u,rho,tau)
%FACC -- computes derivative of Fu (Carslaw and Jaeger, p335, eq 6 below
%integral
%  USAGE:   [fa, fu_fa] = Facc(u,rho,tau)
%     fa = derivative
%     fu_fa = fu/fa, needed in Newton Raphson
%
% SEE ALSO fu, newtonRoot_u
%
% TO 140131

    [fu, J0u, Y0u, J0ru, Y0ru, Du, etu2_u] = Fu(u,rho,tau);

    J1u  = besselj(1,u);
    Y1u  = bessely(1,u);
    J1ru = besselj(1,rho.*u);
    Y1ru = bessely(1,rho.*u);

    fa = -(2*tau.*u+1./u).*fu + etu2_u .* (...
        +1./Du   .*( - rho*J1ru.*Y0u - J0ru.*Y1u + rho*Y1ru*J0u + Y0ru*J1u)...
        +2./Du.^2.*( (J0ru.*Y0u -Y0u.*J0u) .* (J0u.*J1u + Y0u.*Y1u) )...
        );

    fu_fa = fu./fa;
    
end

function [fu,J0u,Y0u,J0ru,Y0ru,Du,etu2_u] = Fu(u,rho,tau)
% Computes Fu, Carlaw and Jaeger, p335, eq 6, below the integral
% USAGE: [fu,J0u,Y0u,J0ru,Y0ru,Du,etu2_u] = Fu(u,rho,tau)
%
%  fu     = function value
%  J0u    = J0(u);      Y0u  = Y0(u)
%  J0ru   = J0(rho*u);  Y0ru = Y0(rho*u)
%  Du     = J0(u).^2 + Y0(u).^2
%  eth2_u = exp(-tau*u.^2)./u
%
% SEE ALSO facc, newtonRoot_u
%
% TO 140131

    J0u    = besselj(0,u);
    Y0u    = bessely(0,u);
    J0ru = besselj(0,rho*u);
    Y0ru = bessely(0,rho*u);
    Du   = J0u.^2+Y0u.^2;
    etu2_u =  exp (-tau *u.^2)./u;
    
    %notice that Peng (2002) has interchanged the central factor, causing the
    %sign to be opposite.
    
    fu = etu2_u .* (J0ru.*Y0u - Y0ru.*J0u)./Du;
    
end
