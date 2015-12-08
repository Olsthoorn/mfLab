function Y=cosspace(x1,x2,N,mode)
% Y=COSSPACE(x1,x2,N,mode)
% produces Y as a series in N steps from x1 to x2 such that  d follows cosine
% mode may be omitted or be set to 'l' or 'r' to get the left or right half of the cosine
% TO 000528

d1=-pi/2; d2=pi/2;

if nargin>3
	switch upper(mode(1))
	case {'L'}
		d2=0;   
	case {'R'}
   	d1=0;
	end
end


r=[d1:(d2-d1)/(N+1):d2]; r=r(2:end-1);
alfa=(x2-x1)/sum(cos(r));
Y=[x1,x1+cumsum(alfa*cos(r))];