function [Section,phi,q,sTop,X,x,kD,c,H,n,Q]=nsecn(x,kD,c,H,n,Q,X)
%NSECN solves multi layer flat analytical model (steady state)
%
% Example:
%    [Section,phi,q,sTop,X,x,kD,c,H,n,Q]=nsecn(x,kD,c,H,n,Q,X);
%
% USAGE: [phi,q,s]=nsecn(x,kD,c,h,n,Q,X)
% USAGE: [phi,q,s[,X[,x[,kD[,c[,h[,n[,Q]]]]]]=nsecn(x,kD,c[,h[,n[,Q[,X]]])
%
% Analytic one-dimensional Mazure solution for nLayers and n-connected sections (using Matrix functions).
% Left and right most sections run to -inf and +inf resp.
%
% input:
%  x, (nNod by  1  ) vector of coordinates of intersection points
% kD, (nLay by nSec) matrix of transmissivity values
%  c, (nLay by nSec) matrix of vertical resistance values of overlaying aquitards
%  h, ( 1   by nSec) vector of given head on top of each section
%  n, ( 1   by nSec) vector of given recharge on top of each section
%  Q, (nNod by nSec) matrix of nodal injections [L2/T]
%  X,  a vector of point where head phi and flow will be computed
%
%  output:
%  phi computed heads [H]    a (nLay by length(X)) matrix
%  q   computed flows [L2/T] a (nLay by length(X)) matrix
%  s   [L/T] is the downward positive seepage rate through top of each layer, a (nLay by length(X) matrix
%  X   is sorted X or 10 points within each inner section when omitted or empty in the input
%  The other outputs equal the inputs, x is sorted and kD, c,h, Q are augmented if necessary
%     these outputs may be utilised to fill out simplified inputs or generate X
%  if Q is empty,  0 will be used for all nodes and layers, a single value per layer is accepted
%  if h is empty,  0 will be used for all sections
%  if n is empty,  0 will be used for all sections
%  if c or kD are empty an error message is generated
%  if colmns of c and kD will be added to the right, but a value for each layer is necessary
%  if number of columns of c, kD, h or Q are > than follows from length(x) an error message is generated 
%
% See also: hantushn shownsecn
%
% T.N.Olsthoorn 990505 001113 001203 (neerslag)

if nargin<7; X=[]; end
if nargin<6; Q= 0; end
if nargin<5; n= 0; end
if nargin<4; H= 0; end

if nargin<3 | isempty( c); error('ERROR nsecn:  c must be specified!'); end
if nargin<2 | isempty(kD); error('ERROR nscen: kD must be specified!'); end
if nargin<1; help('nsecn'); return; end

[x,kD,c,H,n,Q,nLay,nSec,nNod,X]=checkInput(x,kD,c,H,n,Q,X);		% nNod is number of internal nodes = nSec-1
% all points including of the outer sections to infinity!
x=[   x(1),        x,     x(end)      ]; Nx=length(x); %  later we replace x(1) and x(end) bij -Inf and Inf resp.

J=eye(nLay);

% MidSection points are used to compute relative coordinates within sections
% system info for all sections
for iSec=1:nSec
   a=1./(kD(:,iSec).* c(:,iSec));
   b=1./(kD(:,iSec).*[c(2:nLay,iSec);inf]);
   Section(iSec).A=diag(a+b)-diag(a(2:nLay),-1)-diag(b(1:nLay-1),1);	% system matrix
	Section(iSec).RA=sqrtm(Section(iSec).A);   								% matrix roof of system matrix
   Section(iSec).T=diag(kD(:,iSec));											% T matrix
   Section(iSec).h=ones(nLay,1)*H(iSec);										% given heads
   Section(iSec).n=[1;zeros(nLay-1,1)]*n(iSec);								% given precipitation surplus vector
   if iSec<nSec, Section(iSec).Q=Q(:,iSec); end								% injection at right node of section
   Section(iSec).xL=x(iSec);														% left x of section
   Section(iSec).xR=x(iSec+1);													% right x of section
   Section(iSec).xM=0.5*(x(iSec+1)+x(iSec));									% center x of section
   Section(iSec).b =0.5*(x(iSec+1)-x(iSec));									% half width of section (0 for left and right sections)
   Section(iSec).eA=expm(-Section(iSec).b*Section(iSec).RA);			% because xR=b and xL=-b only eA and eB needed
   Section(iSec).eB=expm(+Section(iSec).b*Section(iSec).RA);
  	Section(iSec).fA=Section(iSec).T*Section(iSec).eA*Section(iSec).RA;	% because xR=b and xL=-b only fA and fB are needed
   Section(iSec).fB=Section(iSec).T*Section(iSec).eB*Section(iSec).RA;
   Section(iSec).LTn=eye(nLay)/Section(iSec).A/Section(iSec).T*Section(iSec).n;	% specific precip. surplus
   if iSec==1 | iSec==nSec
      Section(iSec).qnAbs=Section(iSec).T*Section(iSec).RA*Section(iSec).LTn; % store postive value only
   else
      Section(iSec).qnAbs = Section(iSec).T*Section(iSec).RA*tanhm(Section(iSec).b*Section(iSec).RA)*Section(iSec).LTn;
   end
end
x(1)=-Inf; x(end)=Inf; Section(1).xL=-Inf; Section(end).xR=Inf;

% Generating and filling the total coefficient matrix
C=zeros(2*nLay*(nSec-1),2*nLay*nSec);			% Coefficient matrix
R=zeros(2*nLay*(nSec-1),1);							% Right hand side vector
for iSec=1:nSec-1
	   m1=1+(iSec-1)*2*nLay; m2=m1+2*nLay-1; n1=m1; n2=n1+2*nLay-1;
   C(m1:m2,n1:n2)=[+Section(iSec).eA,+Section(iSec).eB;+Section(iSec).fA,-Section(iSec).fB];
	R(m1     :m1+  nLay-1)=R(m1     :m1+  nLay-1)-Section(iSec).h;
    R(m1+nLay:m1+2*nLay-1)=R(m1+nLay:m1+2*nLay-1)-Section(iSec).Q;
    R(m1+nLay:m1+2*nLay-1)=R(m1+nLay:m1+2*nLay-1)-Section(iSec).qnAbs;
end
for iSec=2:nSec
   m1=1+(iSec-2)*2*nLay; m2=m1+2*nLay-1; n1=1+(iSec-1)*2*nLay; n2=n1+2*nLay-1;
	C(m1:m2,n1:n2)=[-Section(iSec).eB,-Section(iSec).eA;-Section(iSec).fB,+Section(iSec).fA];
   R(m1     :m1+  nLay-1)=R(m1     :m1+  nLay-1)+Section(iSec).h;
   R(m1+nLay:m1+2*nLay-1)=R(m1+nLay:m1+2*nLay-1)-Section(iSec).qnAbs;
end

% Solve the system, using all layers and leaving out the outer nLay columns
% These have no freedom, because the sectons run to infinity. COEF(1:nLay) and COEF(end-nLay:end) must be zero.
COEF=C(:,nLay+1:end-nLay)\R;	COEF=[zeros(nLay,1);COEF;zeros(nLay,1)];

phi=zeros(nLay,length(X));			% heads
q  =zeros(nLay,length(X));			% flow [L2/T]
s  =zeros(nLay,length(X));			% seepage upward [L/T]
N  =zeros(nLay,length(X));		% dummy testen neerslag
for iSec=1:nSec
   k=2*nLay*(iSec-1)+1;		l=k+nLay-1;		% k:l indices in COEF matrix for this section
   
   I=find(X>=Section(iSec).xL & X<=Section(iSec).xR);
   u=X(I)-Section(iSec).xM;	% relative coordinates in section
      
   A =Section(iSec).A;
   RA=Section(iSec).RA;
   T =Section(iSec).T;
   b =Section(iSec).b;
   h =Section(iSec).h;
   LTn=Section(iSec).LTn;
   for i=1:length(u)
	   % compute head and flow for this point
   	C1=expm(-u(i)*RA)*COEF(k:l);
	   C2=expm(+u(i)*RA)*COEF(k+nLay:l+nLay);
   	C3=RA*C1;
      C4=RA*C2;
      phi(:,I(i))	=	C1+C2+h;				% head withour effect precipitation
      q  (:,I(i))	=	T   *( C3-C4);		% flux without effect precipitation
	   s	(:,I(i))	=	T*RA*(-C3-C4);		% dq/dx without effect precipitation
      N  (:,I(i)) =  Section(iSec).n;	% infiltration needed for seepage
      if any(Section(iSec).n~=0)			% in case precipitation
	      if iSec==1
   	      phi(:,I(i))=phi(:,I(i))+  (J-expm(+u(i)*RA))*LTn;
      	   q(  :,I(i))=q(  :,I(i))+T*RA*expm(+u(i)*RA) *LTn;
      		s(  :,I(i))=s(  :,I(i))+T* A*expm(+u(i)*RA) *LTn;
	      elseif iSec==nSec
   	      phi(:,I(i))=phi(:,I(i))+  (J-expm(-u(i)*RA))*LTn;  
      	   q(  :,I(i))=q(  :,I(i))-T*RA*expm(-u(i)*RA) *LTn;
         	s(  :,I(i))=s(  :,I(i))+T* A*expm(-u(i)*RA) *LTn;
	      else
   	      phi(:,I(i))= phi(:,I(i))+ (J-coshm(+u(i)*RA)/coshm(b*RA))*LTn;
      	   q(  :,I(i))=q(  :,I(i))+T*RA*sinhm(+u(i)*RA)/coshm(b*RA) *LTn;
         	s(  :,I(i))=s(  :,I(i))+T* A*coshm(+u(i)*RA)/coshm(b*RA) *LTn;
	      end
      end
      % compute seepage through top of each aquifer, starting at lowest aquifer upward
      %   fi=[h(iSec);phi(:,i)];
%   for iLay=1:nLay, ss(iLay,i)=(fi(iLay)-fi(iLay+1))/c(iLay,iSec); end  % checking correctness of s
	end
end
%s=flipud(cumsum(flipud(s)));
sTop=s-N; for i=nLay-1:-1:1, sTop(i,:)=sTop(i,:)+sTop(i+1,:); end
sTop=-sTop;	% line is negative if water flows downward and vice versa
i=0; % comment

function [x,kD,c,h,n,Q,nLay,nSec,nNod,X]=checkInput(x,kD,c,h,n,Q,X);
nNod=length(x);		% nodes without -Inf and Inf not yet added at this point
if x(  1)==-Inf; x=x(2:end  ); end;		% strip -Inf if necessary;
if x(end)==+Inf; x=x(1:end-1); end;		% strip +Inf if necessary;
nSec=length(x)+1;		% number of sections
nLay=size  (kD,1);	% number of layers

fprintf('Running NsecN %s\n',datestr(now));
fprintf('Number of sections is %d as derived from x input.\n',nSec);
fprintf('Number of aquifers is %d as derived from transmissivity input.\n',nLay);
fprintf('Values will be computed for %f<=x<=%f\n',min(X),max(X));
fprintf('Checking input...');

x=sort(x(:)');
if isempty(X);
   J=0:9;  dx=diff(x)/10;
   for i=1:nNod-1;
      X=[X,[x(i)+J*dx(i)]];
   end
else
   X=X(:)';
end

if isempty(Q); Q=0; end
if isempty(h); h=0; else  h=h(1,:); end

if size( c,1)~=nLay; error('kD and c have different number of layers!'); end
if size( c,2)>nSec; error('Too many  c-sections, compared to number of x-values!'); end
if size(kD,2)>nSec; error('Too many kD-sections, compared to number of x-values!'); end
if size( h,2)>nSec; error('Too many  h-sections, compared to number of x-values!'); end
if size( n,2)>nSec; error('Too many  n-sections, compared to number of x-values!'); end
if size( Q,2)>nSec-1; error('Too many Q-points , must be <= nSections-1! (=Q between sections'); end

ns=size(kD,2); for is=ns+1:nSec; kD( :,is)=kD(:,ns); end;		% augment kD if necessary
ns=size( c,2); for is=ns+1:nSec;  c( :,is)= c(:,ns); end;		% augment  c if necessary
ns=size( h,2); for is=ns+1:nSec;  h( 1,is)= h(1,ns); end;      % augment  h if necessary
ns=size( n,2); for is=ns+1:nSec;  n( 1,is)= n(1,ns); end;      % augment  h if necessary
nl=size( Q,1); for il=nl+1:nLay;	 Q(il, :)= Q(nl,:); end;		% augment  Q is necessary
ns=size( Q,2); for is=ns+1:nSec-1;Q( :,is)= Q(ns,:); end;		% augment  Q is necessary

fprintf('ok.\n');

function Z=sinhm(z)
	Z=(expm(z)-expm(-z))/2;
function Z=coshm(z)
   Z=(expm(z)+expm(-z))/2;
function Z=tanhm(z)
   Z=(expm(z)-expm(-z))/(expm(z)+expm(-z));
   
function [x,kD,c,h,n,Q]=testinput
% TEST 1
   x =[  -2000 -1500 -1000 -500  0    500 1000 1500 2000  ];
h = [  -2     2    -2    2   -2    2   -2    2   -2    2];
n = [1e-3  1e-3  1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3];
c =[[ 250   250   250  250  250  250  250  250  250  250];...
    [ 250   250   250  250  250  250  250  250  250  250];...
    [ 250   250   250  250  250  250  250  250  250  250]];
kD=[[ 500   500   500  500  500  500  500  500  500  500];...
    [ 500   500   500  500  500  500  500  500  500  500];...
    [ 500   500   500  500  500  500  500  500  500  500]];
Q =[[    0     0     0    0    0     0     0    0  0    ];...
    [    0     0     0    0    0     0     0    0  0    ];...
    [    0     0     0    0    0     0     0    0  0    ]];
