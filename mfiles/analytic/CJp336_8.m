function H = CJp336_8(varargin)
%CJP335_6 -- computes formula Carslaw & Jaegr, p336, eq8,the heat loss
% through the circumference of a rod at constant unit temp since t=0 in a
% medium that was at zero temperature initially.
% The temperature in the medium can be computed using CJp335_6.
%
% Computing the integral is not trival because the integrand
% has a singularity at u=0. The procedure to compute the
% integral borrowed from Peng, Yeh and Yang (2002) Advances in Water
% Resources, Vol 25, 663-675.
%
% The procedure is as follows:
% 1) take a small value u1 over which the integral is computed analytically
% 2) take a large value u2 as the upper bound of the integration.
%     u1 = UTOL  and u2 is from exp(-tau u^2)<UTOL, i.e. that the 
%     damping exponent is less than UTOL.
% 3) Integrate the integral between u1 and u2 by quad with given tolerance
%    and add the three contributions:
% 4) Compute H=8\lambda/pi sum(q);
%
% USAGE: F = CJp336_8(lambda,tau)
%    lambda = heat conductivity of the medium [W/m/K]
%    tau = kappa/R^2 *t
%      and kappa=lambda/rhoc, the diffusivity of the medium around the rod.
%
% TO 140201

    [lambda,varargin] = getNext(varargin,'double',[]);
    [tau, ~         ] = getNext(varargin,'double',[]);
    if isempty(tau) || isempty(lambda)
        error('need at least two input arguments, lambda and tau');
    end
    
    %% Get the integration interval
    UTOL  = 1e-10;

    u(1) = UTOL;                  % start of integration interval
    u(2) = sqrt(-log(UTOL)/tau);  % end of integration when damping < UTOL
    q    = zeros(1,2);    
    
    % Integrate from 0 to uRt(1) across singularity
    q(1) = -pi/2*atan(besselj(0,u(1))/bessely(0,u(1)));

    % Use quad with given accuracy
    TOL   = 1e-8;               % integration tolerance
    TRACE = true;               % show quad convergence
    
    func = @(u) exp (-tau *u.^2)./u./(besselj(0,u).^2+bessely(0,u).^2);

    q(2) = quad(func,u(1),u(2),TOL,~TRACE);
    
    H = 8*lambda/pi*sum(q);

end
