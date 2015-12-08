% Test the analytical solution
% TO 151122

% We will compare here the steady state solution of the numerical and the
% analytical solution. The latter steady state solution as function of x along
% the cross section, we don't know for given seepage, we only know it for
% given head in the underlying aquifer. What we will thus do is compute the
% numeical model for given seepage, extract from it the head in the
% regional aquifer and apply that to the steady-state analytical solution.
% The heads in both solution should than almost perfectly match.
%
% TO 151126

% This mfile is run after mf_analyze so that the heads of the numerical
% model and all other information of the models is at hand.

%% Numerial steady state solution
% We assume here that the last heads in the numerical model are steady
% state, which we achieve by setting the net recharge for the numerical
% model equal to a given value for at least the last simulation year, while
% also the seepage is constant. We can then extract the steady state heads
% from the numecal model as
Area = sum(gr.AREA .* active,2);

hxNum = H(end).values;
hNum  = sum(H(  end).values,2)./sum(gr.AREA,2);
phiN   = sum(Phi(end).values,2)./sum(gr.AREA,2); % numerical head in regional aquifer
qNum  = sum(B(end).term{strmatchi('WELLS',B(end).label)}(:,:,end),2) ./ Area;
rchN  = sum(B(end).term{strmatchi('RECH' ,B(end).label)}(:,:,  1),2) ./ Area;
etrN  = sum(B(end).term{strmatchi('ET'   ,B(end).label)}(:,:,  1),2) ./ Area;
stoN  = sum(B(end).term{strmatchi('STOR' ,B(end).label)}(:,:,  1),2) ./ Area;
QNum  = sum(B(end).term{strmatchi('HEADDEPB',B(end).label)}(:,:,1),2) + ...
        sum(B(end).term{strmatchi('RIV'  ,B(end).label)}(:,:,  1),2);
Nnum  = rchN - etrN;

cNum  = sum(gr.DZ(:,:,2)./VKCB ,2)     ./ Area; 
kDNum = sum(gr.DZ(:,:,1).*HK(:,:,1),2) ./ Area;
%SyNum = sum(SY(:,:,1), 2)              ./ Area;
SyNum = sum(SS(:,:,1).*gr.DZ(:,:,1),2) ./ Area;

% Compare this with the given seepages
qdiff =  sqrt(mean((q - qNum).^2));
hLR   =  hDitch(:,end);
Nend  = -diff(tne(end,2:3)); % steady state recharge

% try various phi until q equals the desired value

%% Verify the analytic solution

% Compute hx and hbar analytically using the solution with the fixed phi
% and the phi that was computed with Modflow that was fed with a fixed q.

xAn  = bsxfun(@minus, b, gr.xm); % coordinates start at xAn=0 for xNum=b

% Analytical solution for given phi (take phiN from numeric model)
hbar = (phiN + Nend.*c) - ( (phiN + Nend.*c) - hLR) .* ...
       lambda./b .* sinh(b./lambda) ./ ...
       (lambda ./c .* w./D1 .* sinh(b./lambda) + cosh(b./lambda));
   q = (phiN-hbar)./c;
   

% Analytic solution for given phi (taking phiN from numeric model)
hx = bsxfun(@times,  phiN + Nend.*c, active) - ...
     bsxfun(@times, (phiN + Nend.*c-hDitch(:,end)) ./ ...
         ((lambda./c .* w./D1).*sinh(b./lambda) + cosh(b./lambda)), ...
         cosh( bsxfun(@times, xAn, 1./lambda) ) );

% Analyic solution for the outflow (taking phiN from numeric model) 
Q = - (phiN + Nend.*c - hDitch(:,end)) .* sinh(b./lambda) ./ ...
    ((w./D1).*sinh(b./lambda) + (lambda./kD1) .* cosh(b./lambda));

%% Show steady-states solution

figure;
axes(defaults{:});
xlabel('x [m]'); ylabel('head [m]');
title('steady-state head analytic (blue) vs numerical (red)');
I = [1 2];
plot(gr.xm,hx(I,:),'b');
plot(gr.xm,H(end).values(I,:,:),'r');

%% Compare results
display([hbar hNumeric(:,end)]);
display([Q    QNum(    :,end)]);
