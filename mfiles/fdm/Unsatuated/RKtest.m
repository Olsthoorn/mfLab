   

x= (-1:0.01:1)';
a = [1.1275      0.35018     -0.29907      0.02289       -0.262      -1.7502     -0.28565     -0.83137];
y = polyval(a,x);
ya = polyval(polyder(a),x);

xlim = [min(x) max(x)];

figure; axes('nextplot','add','xGrid','on','yGrid','on','xlim',xlim,'fontSize',12);
xlabel('x'); ylabel('y'); title('integration of a polinomial, using Runge Kutta 3rd order');
plot(x,y);
plot(x,ya);
legend('y','yAcc');


t  = -1; dt= 2.0;
h0 = polyval(a,t);
plot(t,h0,'bo');

hEnd = RungeKutta(t,h0,dt);

plot(t+dt,hEnd,'ro');


% %% Newton Raphson to accurately estimate the next ht
% %  Pre estimate ht(1) using Newton Raphson method
% if iOuter == 1
%      for iNR=1:3  % Newton Raphson iterations
%         hNROld     = ht(1);
%         if iOuter==1
%             dTheta    = Ss(strtH(1)) * dh;
%             ht(1)      = ht(1) - yNR(hNROld)/yacNR(hNROld);
%             htEnd(1)   = strtH(1) + (ht(1) - strtH(1))/implicity;
%         else
%             dTheta     = (h2theta(htEnd(1)) - h2Theta(hStrt(1)));
%             dh         = htEnd(1) - ht(1);
%         end
%         %fprintf(' %g',ht(1));
%         if abs(ht(1)-hNROld)<hClose, break; end
%      end
%      htEnd(:)  = strtH + (ht-strtH)/implicity;
%      %fprintf('\n');
% end
% 
