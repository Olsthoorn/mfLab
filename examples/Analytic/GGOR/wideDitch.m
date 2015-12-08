%function [h1,h2,h1_avg,h2_avg,x1,x2]=wideDitch(tne)
%     [h1,h2,h1_avg,h2_avg,x1,x2]=wideDitch(tne)
%
% Analytic multi-layer solution on top of semi-confined aquifer separated
% by leaking semi-confined bed with given supply in second aquifer.
% Situation
%
%                               |||||||||||| q1 |||||||||||||||||||
%           wide ditch       |  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv |
%                            |                                      |
%    |~~~~~~~~  c11~~~~~~~~~~~                                      |
%    |          T11                         T21                     |
%    |//////////c12/////////////////////////c22//////////////////// |
%    |          T12                         T22                     |
%    |                          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ |
%    |                          |||||||||||  q2 ||||||||||||||||||| |
%    |==== x1 ===>                                   <== x2 ========|
%    |<-------- b1 --------->|<-------------- b2 ------------------>|
%
%
%
% Analytical solutions both have constant layer thickness and prescribed
% upward flux through the confining bed.
%
% tne is the recharge and ET during the day. The head at the start of the
% first day (t=0) is set equal to hLR the ditch level on both sides.
%
% TO 101227

%% tne has [time Precip and ET]. Note that these values are during the days

%% Initiaize solution specific
N=0.01;
qu=-0.003;

NLay=2;

h1_hat=-1;
h2_hat= 1;

kD1 = [10; 50];
kD2 = [10; 50];
c1  = [  5; 20; Inf];
c2 =  [1e6; 20; Inf];
b1  =  8;
b2  = 50;
q1   = [0; qu];
q2   = [N; qu];
 
T1   =diag(kD1); T1_m1=T1^(-1);
T2   =diag(kD2); T2_m1=T2^(-1);

%I   = eye(NLay);

A1=-diag( 1./(kD1(2:NLay  ).*c1(2:NLay)),-1)+...
          +diag( 1./(kD1(1:NLay  ).*c1(1:NLay))+1./(kD1(1:NLay).*c1(2:NLay+1)), 0)+...
          -diag( 1./(kD1(1:NLay-1).*c1(2:NLay)),1);
A2=-diag( 1./(kD2(2:NLay  ).*c2(2:NLay)),-1)+...
          +diag( 1./(kD2(1:NLay  ).*c2(1:NLay))+1./(kD2(1:NLay).*c2(2:NLay+1)), 0)+...
          -diag( 1./(kD2(1:NLay-1).*c2(2:NLay)),1);

A1_m1=A1^(-1);
A2_m1=A2^(-1);

sqrtA1 = sqrtm(A1); sqrtA1_m1=sqrtA1^(-1); % sqrtA_m1= sqrtA^(-1);
sqrtA2 = sqrtm(A2); sqrtA2_m1=sqrtA2^(-1); % sqrtA_m1= sqrtA^(-1);

sinhm1=funm(b1*sqrtA1,@sinh); % sinhm_m1=sinhm^(-1);
coshm1=funm(b1*sqrtA1,@cosh); % coshm_m1=coshm^(-1);

sinhm1_m1=sinhm1^(-1); % sinhm_m1=sinhm^(-1);
coshm1_m1=coshm1^(-1); % coshm_m1=coshm^(-1);

sinhm2=funm(b2*sqrtA2,@sinh); % sinhm_m1=sinhm^(-1);
coshm2=funm(b2*sqrtA2,@cosh); % coshm_m1=coshm^(-1);

sinhm2_m1=sinhm2^(-1); % sinhm_m1=sinhm^(-1);
coshm2_m1=coshm2^(-1); % coshm_m1=coshm^(-1);

B = (coshm1*sqrtA1_m1*sinhm1_m1*T1_m1+coshm2*sqrtA2_m1*sinhm2_m1*T2_m1)^(-1);
Q=B*(h2_hat-h1_hat+A2_m1*T2_m1*q2-A1_m1*T1_m1*q1);

h1_avg=h1_hat+A1_m1*T1_m1*(q1+Q/b1);
h2_avg=h2_hat+A2_m1*T2_m1*(q2-Q/b2);

x1=0:1:b1;
x2=0:1:b2;

h1=zeros(NLay,length(x1));
for i=1:length(x1)
    h1(:,i)=h1_hat+funm(x1(i)*sqrtA1,@cosh)*sqrtA1_m1*sinhm1_m1*T1_m1*Q+A1_m1*T1_m1*q1;
end

h2=zeros(NLay,length(x2));
for i=1:length(x2)
    h2(:,i)=h2_hat-funm(x2(i)*sqrtA2,@cosh)*sqrtA2_m1*sinhm2_m1*T2_m1*Q+A2_m1*T2_m1*q2;
end

figure; hold on; grid on; xlabel('x [m]'); ylabel('head [m]');
title('wide ditch');

leg='';
plot(      x1,h1(1,:),'r'  ); leg{end+1}='h_d(1)';
plot(      x1,h1(2,:),'r.-'); leg{end+1}='h_d(2)';

plot(b1+b2-x2,h2(1,:),'b');   leg{end+1}='h_2(1)';
plot(b1+b2-x2,h2(2,:),'b.-'); leg{end+1}='h_2(2)';

plot(x1,      h1_hat   *ones(size(x1)),'r','linewidth',2); leg{end+1}='h_{ditch}';
plot(x1,      h1_avg(1)*ones(size(x1)),'rx');              leg{end+1}='h_{dAvg}(1)';
plot(x1,      h1_avg(2)*ones(size(x1)),'r+');              leg{end+1}='h_{dAvg}(2)';

%plot(b1+b2-x2,h2_hat   *ones(size(x2)),'b','linewidth',2); leg{end+1}='h_{right}';
plot(b1+b2-x2,h2_avg(1)*ones(size(x2)),'bx');              leg{end+1}='h_{2Avg}(1)';
plot(b1+b2-x2,h2_avg(2)*ones(size(x2)),'b+');              leg{end+1}='h_{2Avg}(2)';

legend(leg(2:end));

