function W=boulton(t,r,he,z,kv,alpha,Sa,Sy)
% argument of Boulton type curves
% TO 120113 % not correct !!

rho=r/he;
eta=(Sa+Sy)/Sa;
zeta=z/he;
gamma=kv/(alpha*Sa*he);

t=t(:); r=r(:)';

W=NaN(length(t),length(r));

for i=1:length(r)
    for j=1:length(t)
        W(j,i)=2*quadgk(@boultonArg,0,Inf);
    end
end

    function b=boultonArg(u)
        gthu=gamma*u.*tanh(u);
      
        f1=alpha*t*(eta+gthu)/2;
        f2=alpha*t*(eta-gthu)/2;
        f3=alpha*t*sqrt((gthu+eta).^2-4*gthu)/2;

        Ft=exp(-f1).*(cosh(f3)+f2./f3.*sinh(f3));
        
        b=besselj(0,u*rho)./u.*(1-cosh(u*(1-zeta))./cosh(u).*Ft);
    end
end