function F = logQuad(func,u1,u2,TOL,NSTEPS)
   %%LOGQUAD -- trapezium integration but integrates with exponent
   % func must be positive and declining
   %
   % USAGE: A = logQuad(func,u1,u2[,TOL[,NSTEPS]])
   %   func = function handle as in quad
   %   u1 and u2 integration interval
   %   TOL abs tolerance default 1e-6
   %   NSTEPS max nr of iterations default 100
   %
   % TO 140201
   
   if nargin<5, TOL=1e-6; end
   if nargin<6, NSTEPS=100; end

    A=zeros(1,NSTEPS);
    A(1) = int(func,u1,u2,1,TOL);

    for i=2:NSTEPS
        A(i) = int(func,u1,u2,i,TOL); 
        if abs(A(i)-A(i-1))<TOL
            F=A(i);
            return;
        end
    end
    
    fprintf('No convergence in %d steps\n',NSTEPS);
end
function A = int(func,u1,u2,n,TOL)
    du = (u2-u1)/n;
    y = max(TOL,-func(u1:du:u2));
    lambda = du./log(y(1:end-1)./y(2:end));
    dA = lambda.*y(1:end-1).*(1-exp(-du./lambda));
    A  = sum(dA);
end
    
