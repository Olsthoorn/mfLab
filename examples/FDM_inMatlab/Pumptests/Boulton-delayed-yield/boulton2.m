function W=boulton2(t,r,T,alpha,Sa,Sy)
% argument of Boulton type curves
% TO 120113 % not correct !!

lamb=sqrt(T/(alpha*(Sa+Sy)));

rho=r/lamb;

eta=(Sa+Sy)/Sa;

t=t(:); r=r(:)';

W=NaN(length(t),length(r));

for i=1:length(r)
    for j=1:length(t)
        W(j,i)=2*quad(@boultonArg,1,2);
    end
end

    function w=boultonArg(x)
        mu1=alpha*t*eta*(1-x.^2)/2;
        mu2=alpha*t*sqrt(eta^2*(1+x.^2).^2-4*eta*x.^2);
        
        w=2*besselj(0,rho*x)./x.*(1-exp(-mu1).*...
            (cosh(mu2)+alpha*t*eta*(1-x.^2)./(2*mu2).*sinh(mu2)));
    end
end