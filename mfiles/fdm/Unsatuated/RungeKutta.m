function h = RungeKutta(t,h0,dt)
    % recursieve Runge-Kutta integratie, 3rd order, O(h^4)
    % Call 3rd order RK function to estimate the next point
    % Verifies that the direct esitmation equals the one done in two
    % stpes of half the time step length each. If not it calls itself
    % with halving the step. This recursive function automatically
    % subdivides the distance to the extent necesassary to attain the
    % required accuracy
    % TO 150928
    dh = 0.00001;
    
    h1 = RK3(t,h0,dt);           plot([t,t+dt],[h0,h1],'ro-');
    h2 = RK3(t,h0,dt/2);         plot([t,t+dt/2],[h0,h2],'bo-');
    h3 = RK3(t+dt/2,h2,dt/2);    plot([t+dt/2,t+dt],[h2,h3],'go-');
    Delta = (h3-h1);
    if abs(Delta)>dh
        h1 = RungeKutta(t      ,h0,dt/2);
        h2  = RungeKutta(t+dt/2,h1,dt/2);
        h   = h2;
    else
        h=h3;
    end
 end 
    function h = RK3(t,h0, dt)
        %% 3rd order Runge-Kutta integration (Abramowitz and Stegun (1964) 25.5.8)
        k1 = dt * hAcc(t, h0);
        k2 = dt * hAcc(t+dt/2, h0 + k1/2);
        k3 = dt * hAcc(t+dt   ,h0 - k1 + 2*k2);
        h  = h0 + k1/6 + 2*k2/3 + k3/6;
    end

    function hac = hAcc(x,~)
         a = [1.1275      0.35018     -0.29907      0.02289       -0.262      -1.7502     -0.28565     -0.83137];
         %y = polyval(a,x);
        hac = polyval(polyder(a),x);
    end
    
% % Voor Newton Raphson
% function y = yNR(h)
%     y = Cz(1) * (h-ht(2)) + Ccauchy(1) * (h - hCauchy(2)) - P(it)  ...
%         +dz(1)/dtau(itau) * (h2theta(htEnd(1)) - h2theta(strtH(1))) ...
%         +1/dtau(itau) * (Wpond(h) - Wpond(strtH(1)));
%     end
%     function yacc = yacNR(h)            
%             yacc = Cz(1) + Ccauchy(1) + dz(1)/dtau(itau) * Ss1(h) + Cpond(h)/dtau(itau);
%     end
