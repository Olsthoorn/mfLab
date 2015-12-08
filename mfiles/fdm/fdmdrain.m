function [Phi,Q,QL]=fdmdrain(dL,D,k,temp,FH,FQ)
%FDMDRAIN a block-centred steady-state finite difference model for flow in drains
%
% Example:
%    [Phi,Q,Qx]=fdmdrain(dL,D,k,temp,QInit,FH,FQ)
%
% INPUT:
%    dL  = section lengths coordinates
%    D  = pipe diameter per section (nodes)
%    k  = wall roughness (m)
%    temp = celcius
%    FH = fixed heads (NaN for ordinary points), at least one non zero
%     the imaginary part (if present) is the cicumferal resisance in days
%    FQ = fixed nodal flows (m3/d), may be all zero, inflow positive, per node
% OUTPUT:
%    Phi= computed heads in nodes
%    Q  = computed inflows of nodes (m3/d)
%    QL = flow between nodes (Nx-1) (m3/d)
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 090428 090508


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

g=9.81; % gravity

FQ=FQ(:)./24/3600; % change from m3/d to m3/s
FH=FH(:);
D =D (:);
dL=dL(:); NL=length(dL);
A =pi*D.^2/4;
D5=D.^5;

Ifh = find(~isnan(FH) & imag(FH)==0);  % fixed heads
I   = find( isnan(FH) | imag(FH)> 0);  % heads to be computed (active cells)
Ighb= find(~isnan(FH) & imag(FH)> 0);  % nodes with general head boundary

Phi=FH;
R=8.*0.025*0.01./(g*pi^2).*(dL(1:end-1)./D5(1:end-1)+dL(2:end)./D5(2:end))/2; % Initial section resistance

% General head boundary, having circumferal resistance c stored in imag of FH
Cghb = zeros(size(FH));
if ~isempty(Ighb)
    Cghb(Ighb)=pi*D.*dL./imag(FH(Ighb)); % abs(imag(FH( )) is c (resisance of circumference in d)
end

m=0;
while 1
    M=sparse([1:NL-1,  2:NL]',...
             [2:NL  ,1:NL-1]',...
            -[1./R  ;  1./R] ,...
             NL,NL,3*NL);   % System matrix
    Mdiag= -sum(M,2);       % Main diagonal
    
    if ~isempty(I) && ~isempty(Ifh)
        Phi(I)=spdiags(Mdiag(I)+Cghb(I),0,M(I,I))\...
            (FQ(I)-M(I,Ifh)*FH(Ifh)+Cghb(I).*FH(I)); % solve
    end
    Q =spdiags(Mdiag,0,M)*Phi;		% Nodal inflows all nodes
    QL=-diff(Phi)./R;              % Internonal piple flows

    R1=8/(g*pi^2)*abs(QL).*...
    (Colebrook(abs(QL)./A(1:end-1),D(1:end-1),k,temp).*dL(1:end-1)./D5(1:end-1)+...
     Colebrook(abs(QL)./A(2:end  ),D(2:end  ),k,temp).*dL(2:end  )./D5(2:end  ))/2;
 m=m+1;
 fprintf('Iteration %d ',m');
 fprintf(' %g',Colebrook(abs(QL)./A(1:end-1),D(1:end-1),k,temp))
 fprintf('\n');

    if max(abs(log(R./R1)))<0.01; % allow 3% dev of friction coeff,
        break;
    else
       R=R1; % Section resistance
    end
end

Q =Q *24*3600; % back to m3/d
QL=QL*24*3600; % back to m3/d

function lamb=Colebrook(vm,D,k,Temp)
% Colebrook's pipe friction coefficient

    nu=0.0005103/((Temp+43.103).^(1.5017));
    
    Re=max(2000,abs(vm).*D/nu);
    
    lamb=0.02;    
    while 1
        lamb1=0.25./(log10(1./(0.4.*Re.*sqrt(lamb))+k/3.7./D).^2);  
        if max(abs(lamb1./lamb-1))<0.0001
            break;
        end
        lamb=lamb1;
    end
    lamb=lamb1;
