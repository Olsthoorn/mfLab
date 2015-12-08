function Y=dspace(x1,x2,d,N)
% Y=DSPACE(x1,x2,d,N)
% produces Y as a series in N steps from x1 to x2 such that  d=(y(i+1)-(y(i))/(y(i)-y(i-1))=constant
% to produce a series between 0 and 10 in 5 steps such that d=1.5
% Y=dspace(0,10,1.5,5);
r=d.^[0:N-1];
alfa=(x2-x1)/sum(r);
Y=[x1,x1+cumsum(alfa*r)];