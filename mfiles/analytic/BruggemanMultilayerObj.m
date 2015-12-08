classdef BruggemanMultilayerObj
%BRUGGEMANMULTILAYEROBJ class deff for objects of analytical solutions by Bruggeman (1999)
%
% Class defefinition to solve all multilayer solutions from Bruggeman (1999)
% of Geohydrological Problems, Elsevier. ISBN 0-444-81829-4
% Cases 710-720, pages 433-449.
% See BruggemanMultilayer in mflab>examples>Analytic for usage
%
% See also: hantushn stehfest
%
% TO 110803
    properties
        maincase
        subcase   
        s  % [title]
        b  % [ m ] width of the system (distance to right hand boundary
        z0 % [ m ] top of aquifer aquitard system
        NLay
        t  % [ d ];
        x  % [ m ]; distance
        r  % [ m ]; distance in axi-symmetric flow cases
        d  % [ m ] aquitard thickness
        D  % [ m ] aquifer thickness
        q  % [m2/d] discharge
        Q  % [m2/d] discharge flat
        h  % [m3/d] extraction radial
        k  % [m/d] conductivity
        kD % [m2/d] transmissivity
        c  % [ d ] resistance confining beds
        w  % [m/d] entry resistance
        A;     A_m1
        sqrtA; sqrtA_m1
        cS;    cS_m1
        T;     T_m1;  % trasmissivity matrix
        S;     S_m1       % storage ceofficient matrix
        H;     H_m1       % entry resistance matrix
        ATm1q        % simple vector in many solutions
        sqrtATm1Q    % simple vector in many solutions
        sqrtATm1q    % same used in axial symmetric solutions (720 series)
        snhm
        snhm_m1
        cshm
        cshm_m1
        DZ
        Z
        Iz
        Phi
        Qx
    end
    methods
        function obj=BruggemanMultilayerObj(nlay)
            %BRUGGEMANMULTILAYEROBJ constructor for objects of analytical solutions by Bruggeman (1999)
            %
            %  Example:
            %     BruggemanMultilayerObj(nlay);
            %     parameters are preset below
            %     results published in Stromingen (april 2013).
            %
            % For convenience and as short-cut we put the data here. We may
            % better put data somewhere else, like an excel workbook and
            % read it in form there. This can be done using the function
            % GetExcelData in mflab or using readxls function in Matlab.
            % By setting NLay we use just NLay layers, even though more are
            % defined in this overview.
            %
            % See also hantushn stehfest
            %
            obj.b  = 50;  % width of the system (distance to right hand boundary
            obj.z0 = 0;  % top of aquifer aquitard system
            obj.d  = [   5   10    5  10   2   10    5    2     5   10]; % aquitard thickness
            obj.D  = [  15   20   15  20  20   15   10   10    15   20]; % aquifer thickness
            obj.q  = [   5   -2    1  -1   3   10   -5    7     1   -6]*0.001; % discharge
            obj.Q  = [   1    2    0  -3   1   -2    1    2   -10    3]; % discharge flat
            obj.h  = [   1   -1    0  -3   1   -2    1    2   -10    3]; % extraction radial
            obj.k  = [   2    2    5  10  15    1    5   10   0.5   10]; % conductivity
            obj.c  = [100    10  200  40 100    5   10    3    10  100]; % resistance confining beds
            obj.w  = [   1    5    1   2   2    3    1    1     2    1]; % entry resistance
            obj.S  = [   1    1    1   1   1    1    1    1     1    1]*1e-3; % elastic storage
            
            if nargin==0, nlay=5; end  % choose n layers
        
            obj.NLay=nlay;        
            obj.d=obj.d(1:nlay);
            obj.D=obj.D(1:nlay);
            obj.q=obj.q(1:nlay);
            obj.Q=obj.Q(1:nlay);
            obj.h=obj.h(1:nlay);
            obj.k=obj.k(1:nlay);
            obj.c=obj.c(1:nlay);
            obj.w=obj.w(1:nlay);
            obj.S=obj.S(1:nlay);

            coshm= @(x) (expm(x)+expm(-x))/2;
            sinhm= @(x) (expm(x)-expm(-x))/2;

            obj.kD = obj.k.*obj.D;         % transmissivity

            obj.DZ=reshape([obj.d;obj.D],[2*length(obj.d),1]); obj.Z=obj.z0-[0; cumsum(obj.DZ)];
            obj.Z=obj.Z(1:2*nlay+1);
            obj.Iz=sort([1:nlay 1:nlay]);
            obj.t=logspace(-2,1,41);

            %% Derive matrices

            obj.T  = diag(obj.kD(1:nlay));   obj.T_m1=obj.T^(-1);                    % transmissivity matrix
            obj.S_m1= diag(1./obj.S);                                                % storage ceofficient matrix
            obj.H   = diag(1./(obj.k(1:nlay).*obj.w(1:nlay))); obj.H_m1= obj.H^(-1); % entry resistance matrix

            % system matrix
            obj.A=-diag( 1./(obj.kD(2:nlay  ).*obj.c(2:nlay)),-1)+...
              +diag( 1./(obj.kD(1:nlay  ).*obj.c(1:nlay))+1./(obj.kD(1:nlay).*[obj.c(2:nlay) Inf]), 0)+...
              -diag( 1./(obj.kD(1:nlay-1).*obj.c(2:nlay)),1);
          
            obj.sqrtA  = sqrtm(obj.A);
            obj.ATm1q=obj.A\(obj.T_m1*obj.q(1:nlay)');              % simple vector in many solutions
            obj.sqrtATm1Q=(obj.sqrtA\(obj.T_m1*obj.Q(1:nlay)'));    % simple vector in many solutions
            obj.sqrtATm1q=(obj.sqrtA\(obj.T_m1*obj.q(1:nlay)'));    % same used in axial symmetric solutions (720 series)

            obj.snhm=sinhm(obj.b*obj.sqrtA); obj.snhm_m1=obj.snhm^(-1); % need these in solutions
            obj.cshm=coshm(obj.b*obj.sqrtA); obj.cshm_m1=obj.cshm^(-1); % need these in solutions

            obj.cS=diag(obj.c(1:nlay).*obj.S(1:nlay)); obj.cS_m1=obj.cS^(-1); % transient (710.01 and 710.03)
                    
        end
        
        function obj=solve(obj,caseid)
            % Works for one-diemsional and axially symmetric cases and
            % deals with zero dimensional cases (only vertical flow 710.01,
            % 710.02 and 710.03) separately
                        
            % extract maincase and subcase from caseid
            if ischar(caseid)
                caseid=round(str2double(caseid));
            end
            obj.maincase=floor(caseid);
            obj.subcase =round(100*(caseid-obj.maincase));
            
            fprintf('Solving Bruggeman Multilayer case %d.%d\n',obj.maincase,obj.subcase);
            
            % define inline functions
            coshm= @(x) (expm(x)+expm(-x))/2;
            sinhm= @(x) (expm(x)-expm(-x))/2;
            
            % use nlay for convenience of writing
            nlay=obj.NLay;

            
            switch obj.maincase    %% Solution 710.01
                case 710
                    
                    % Dealin withg 710.01, 710.02 and 710.03 separately
                    if obj.subcase<10
                        Nt=length(obj.t);

                        obj.A=-diag( 1./(obj.S(2:nlay ).*obj.c(2:nlay)),-1)+...
                              +diag( 1./(obj.S(1:nlay  ).*obj.c(1:nlay))+1./(obj.S(1:nlay).*[obj.c(2:nlay) Inf]), 0)+...
                              -diag( 1./(obj.S(1:nlay-1).*obj.c(2:nlay)),1);

                        obj.A_m1=obj.A^(-1);
                        htop=obj.h; htop(2:end)=0;
                    
                        switch obj.subcase
                            case 1
                                % voldoet aan de differentiaalvergelijking, maar de dimensie van A is hier
                                % 1/d want de vector h heeft elementen h/(c S).

                                obj.s={ 'Bruggeman (1999) Solution 710.01. Transient.';
                                    'n-layer system under open water with sudden rise h of its level';
                                    'phi=phi(t). Only vertical flow in semi-pervious layers.';
                                    };
                                obj.Phi=NaN(nlay,Nt);
                                for it=1:length(obj.t)
                                    obj.Phi(:,it)=(eye(nlay)-expm(-obj.A*obj.t(it)))*obj.A_m1*obj.cS_m1*htop';
                                    obj.Qx = [];
                                end

                            case 2 %% Solution 710.02

                                obj.s={ 'Bruggeman (1999) Solution 710.01. Steady state.';
                                    'n-layer system under open water with sudden rise h of its level';
                                    'phi=phi(t). Only vertical flow in semi-pervious layers.';
                                    };

                                obj.Phi=obj.A_m1*obj.cS_m1*htop';
                                obj.Qx  =[];
                        end
                        return; % finished, skip rest
                    end
                    
                    % Coordaninates for one-dimensional cases
                    L=obj.b;
                    switch obj.subcase
                        case {1,2,11,12,13,14,15,16,17}
                            obj.x=sinespace(0,L,100,0,pi/2);
                        case {18}
                            obj.x=unique([sinespace(0,L,100,pi/2,0) sinespace(L,2*L,100,0,pi/2)]);
                            obj.x=obj.x(:)';
                        case {21,22,23,24,26,27,28}
                            obj.x=sinespace(0,L,100,pi/2,0);
                        case {25}
                            obj.x=sinespace(0,L,100,0,pi/2);
                        otherwise
                            error('Unknown case %d.%d\n',obj.maincase,obj.subcase);
                    end
                    Nx=length(obj.x);
                    obj.Phi = NaN(nlay,Nx);  % preallocate
                    obj.Qx  = NaN(nlay,Nx);  % preallocate
                    
                    % rest of the one-dimensional subcases
                    switch obj.subcase
                        case 12 %% Solution 710.12

                            obj.s={ 'Bruggeman (1999) Solution 710.12. Steady state.';
                                'All aquifers with open boundary. Sudden drawdown of the surface water level,';
                                'which is kept constant thereafter. phi=phi(obj.x)=drawdown';
                                };

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=      expm(-obj.x(ix)*obj.sqrtA)*obj.h';
                                obj.Qx( :,ix)=obj.T*expm(-obj.x(ix)*obj.sqrtA)*obj.sqrtA*obj.h';
                            end

                        case 13  %% Solution 710.13

                            obj.s={ 'Bruggeman (1999) Solution 710.13. Steady state.';
                                'All aquifers with open boundary. The plane at x=0 ketp at zero head.';
                                'Constant vertical infiltration q into the aquifers.';
                                }; % phi=phi(obj.x)=drawdown';

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=     (eye(nlay)-expm(-obj.x(ix)  *obj.sqrtA))*obj.ATm1q;
                                obj.Qx( :,ix)=-obj.T*expm(-obj.x(ix)*obj.sqrtA)*obj.sqrtA  *obj.ATm1q;
                            end

                        case 15 %% Solution 710.15

                            obj.s= { 'Bruggeman (1999) Solution 710.14. Steady state.';
                                 'Fully penetrating well screens in all aquifers at x=0.';
                                 'Sudden discharges which are kept constant theraafter.';
                                }; % phi=phi(obj.x)=drawdown';

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=       expm(-obj.x(ix)*obj.sqrtA)*          obj.sqrtATm1Q/2;
                                obj.Qx(:,ix) = obj.T*expm(-obj.x(ix)*obj.sqrtA)*obj.sqrtA*obj.sqrtATm1Q/2;
                            end

                        case  16  %% Solution 710.16

                            obj.s={ 'Bruggeman (1999) Solution 710.16. Steady state.';
                                'All aquifers with open boundary with entrance resistance.';
                                'Sudden drawdown of the surface water level which i kept constant thereafter.';
                                }; % phi=phi(obj.x)=drawdown';

                            sqrtApHm1Hh=(obj.sqrtA+obj.H)\(obj.H*obj.h');

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=       expm(-obj.x(ix)*obj.sqrtA)*          sqrtApHm1Hh;
                                obj.Qx(:,ix)  =obj.T*expm(-obj.x(ix)*obj.sqrtA)*obj.sqrtA*sqrtApHm1Hh;
                            end

                        case  17  %% Solution 710.17

                            obj.s={ 'Bruggeman (1999) Solution 710.17. Steady state.';
                                'All aquifers with open boundary with entrance resistance.';
                                'Constant infiltration. Zero head at x=0.'
                                }; % phi=phi(obj.x)=drawdown';

                            B=(obj.sqrtA+obj.H)\(obj.H*obj.ATm1q);

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)= obj.ATm1q-expm(-obj.x(ix)*obj.sqrtA)*B;
                                obj.Qx( :,ix)=    -obj.T*expm(-obj.x(ix)*obj.sqrtA)*obj.sqrtA*B;
                            end

                        case 18   %% Solution 710.18

                            obj.s={ 'Bruggeman (1999) Solution 710.18. Steady state.';
                                'Infiltration into aquifers with q(obj.x)=q for |x|<b en 0 for |x|>b.';
                                'Flux=0 for x=0.';
                                }; % phi=phi(obj.x)=drawdown';

                            bexpm  = expm(-obj.b*obj.sqrtA);

                            for ix=1:length(obj.x)
                                if abs(obj.x(ix))<=obj.b
                                    obj.Phi(:,ix)=(eye(nlay)-coshm(obj.x(ix)*obj.sqrtA)*          bexpm)*obj.ATm1q;
                                    obj.Qx( :,ix)= obj.T*    sinhm(obj.x(ix)*obj.sqrtA)*obj.sqrtA*bexpm *obj.ATm1q;
                                else
                                    obj.Phi(:,ix)=      expm(-abs(obj.x(ix))*obj.sqrtA)*          obj.snhm*obj.ATm1q;
                                    obj.Qx( :,ix)=obj.T*expm(-abs(obj.x(ix))*obj.sqrtA)*obj.sqrtA*obj.snhm*obj.ATm1q;
                                end
                            end

                        case 21   %% Solution 710.21

                            obj.s={ 'Bruggeman (1999) Solution 710.25. Steady state.';
                                'Given drawdown h for x=b. Flux=0 for x=0.';
                                 }; % phi=phi(obj.x)=drawdown';

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=       coshm(obj.x(ix)*obj.sqrtA)*          obj.cshm_m1*obj.h';
                                obj.Qx( :,ix)=-obj.T*sinhm(obj.x(ix)*obj.sqrtA)*obj.sqrtA*obj.cshm_m1*obj.h';
                            end

                        case 22  %% Solution 710.22

                            obj.s={ 'Bruggeman (1999) Solution 710.22. Steady state.';
                                'Given drawdown h for x=b. Zero drawdown at x=0.';
                                 }; % phi=phi(obj.x)=drawdown';

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=       sinhm(obj.x(ix)*obj.sqrtA)*          obj.snhm_m1*obj.h';
                                obj.Qx( :,ix)=-obj.T*coshm(obj.x(ix)*obj.sqrtA)*obj.sqrtA*obj.snhm_m1*obj.h';
                            end

                        case 23  %% Solution 710.23

                            obj.s={ 'Bruggeman (1999) Solution 710.23. Steady state.';
                                'Zero head at x=b. Zero flux at x=0. Constant infitlration.';
                                 }; % phi=phi(obj.x)=drawdown';

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=(eye(nlay)-coshm(obj.x(ix)*obj.sqrtA)*     obj.cshm_m1)*obj.ATm1q;
                                obj.Qx( :,ix)=obj.T*sinhm(obj.x(ix)*obj.sqrtA)*obj.sqrtA*obj.cshm_m1 *obj.ATm1q;
                            end

                        case 24 %% Solution 710.24

                            obj.s={ 'Bruggeman (1999) Solution 710.24. Steady state.';
                                'Fully penetrating well screens in all aquifers at x=0, -b, 3b, -3b etc.';
                                'Constant but different discharges q. Flux=0 for x=0.';
                                 }; % phi=phi(obj.x)=drawdown';

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=       coshm(obj.x(ix)*obj.sqrtA)*          obj.snhm_m1*obj.sqrtATm1Q/2;
                                obj.Qx( :,ix)=-obj.T*sinhm(obj.x(ix)*obj.sqrtA)*obj.sqrtA*obj.snhm_m1*obj.sqrtATm1Q/2;
                            end

                        case 25 %% Solution 710.25

                            obj.s={ 'Bruggeman (1999) Solution 710.25. Steady state.';
                                'Fully penetrating well screen in all aquifers in the middle of a strip with width 2b.';
                                'Constant but different dischargfes q. Drawdown=0 for x=0 and x=2b.'
                                }; % phi=phi(obj.x)=drawdown';

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=       sinhm(obj.x(ix)*obj.sqrtA)*          obj.cshm_m1*obj.sqrtATm1Q/2;
                                obj.Qx( :,ix)=-obj.T*coshm(obj.x(ix)*obj.sqrtA)*obj.sqrtA*obj.cshm_m1*obj.sqrtATm1Q/2;
                            end

                        case 26

                            obj.s={ 'Bruggeman (1999) Solution 710.25. Steady state. ';
                                'All aquifers with entrance resistances w_i to open water with a constant drawdown.';
                                'Zero flux at x=0. phi=phi(obj.x)=drawdown.'
                                };

                            F_m1=(obj.H_m1*obj.sqrtA*obj.snhm+obj.cshm)^(-1);

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=       coshm(obj.x(ix)*obj.sqrtA)*          F_m1*obj.h';
                                obj.Qx( :,ix)=-obj.T*sinhm(obj.x(ix)*obj.sqrtA)*obj.sqrtA*F_m1*obj.h';
                            end

                        case  27  %% Solution 710.27

                            obj.s={ 'Bruggeman (1999) Solution 710.27. Steady state.';
                                'All aquifers with entrance resistance w_i to oen water at x=b.';
                                'Constant drawdown of the open water level. Zero drawdown at x=0.'
                                }; % phi=phi(obj.x)=drawdown';

                            F_m1=(obj.H_m1*obj.sqrtA*obj.cshm+obj.snhm)^(-1);

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=       sinhm(obj.x(ix)*obj.sqrtA)*          F_m1*obj.h';
                                obj.Qx( :,ix)=-obj.T*coshm(obj.x(ix)*obj.sqrtA)*obj.sqrtA*F_m1*obj.h';
                            end

                        case 28  %% Solution 710.28

                            obj.s={
                                'Bruggeman (1999) solution 710.28. Steady state.';
                                'Alle aquifers with entrance resistance w to open water at x=b.';
                                'Vertical infiiltration q. Flux=0 for x=0.'
                                };  % phi=phi(obj.x)=drawdown';

                            F_m1=(obj.H_m1*obj.sqrtA*obj.snhm+obj.cshm)^(-1);

                            for ix=1:length(obj.x)
                                obj.Phi(:,ix)=(eye(nlay)-coshm(obj.x(ix)*obj.sqrtA)*      F_m1)*obj.ATm1q;
                                obj.Qx( :,ix)= obj.T*sinhm(obj.x(ix)*obj.sqrtA)*obj.sqrtA*F_m1 *obj.ATm1q;
                            end
                        otherwise
                            error('Unknown solution %d.%d',obj.maincase,obj.subcase);
                    end
                    return;
                    
                case 720  % Axial symmetric cases
                    
                    % axially symmetric coordinates
                    R=obj.b;
                    switch obj.subcase
                        case {1,2,3}
                            obj.r=sinespace(0,R,100,0,pi/2);
                        case {11,12,13,14,15,16,17}
                            obj.r=sinespace(R,2*R,100,0,pi/2);
                        case {21,22,23,24}
                            obj.r=sinespace(0,R,100,pi/2,pi);
                        case {31,32,33}
                            obj.r=sinespace(0,R,100,0,pi/2);
                        case 34
                            obj.r=sinespace(0,R,100,0,pi/2);
                        case 35
                            obj.r=unique([sinespace(0,R,100,pi); sinespace(R,2*R,100,pi/2)]);
                            obj.r=obj.r(:)';
                        otherwise
                            error('Unknown case %d.%d\n',obj.maincase,obj.subcase);
                    end
                    Nr=length(obj.r);
                    obj.Phi = NaN(nlay,Nr);  % preallocate
                    obj.Qx  = NaN(nlay,Nr);  % preallocate
                    
                    % eigen values/vectors and bessel functions needed
                    [V,D_]=eig(R*obj.sqrtA); V_m1=V^(-1);
                    K0R=V*diag(besselk(0,diag(D_)))*V_m1; K0R_m1=K0R^(-1);
                    K1R=V*diag(besselk(1,diag(D_)))*V_m1; K1R_m1=K1R^(-1);
                    I0R=V*diag(besseli(0,diag(D_)))*V_m1; I0R_m1=I0R^(-1);
                    I1R=V*diag(besseli(1,diag(D_)))*V_m1; I1R_m1=I1R^(-1);

                    % axially symmetric subcases
                    switch obj.subcase
                        case 1 %% Solution 720.01

                            obj.s={
                                'Bruggeman (1999) solution 720.01. Steady state.';
                                'Vertical infiiltration q(obj.r). q(obj.r)=q for 0<r<=R and 0 for r>R.'
                                };  % phi=phi(obj.x)=drawdown';

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                if obj.r(ir)<=R
                                    obj.Phi(:,ir)=(eye(nlay)-          R*obj.sqrtA*K1R*v*diag(besseli(0,d_))*v_m1)*         obj.ATm1q;
                                    obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*R*obj.sqrtA*K1R*v*diag(besseli(1,d_))*v_m1*obj.sqrtA*obj.ATm1q;
                                else
                                    obj.Phi(:,ir)=                     R*I1R*v*diag(besselk(0,d_))*v_m1*          obj.sqrtATm1q;
                                    obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*R*I1R*v*diag(besselk(1,d_))*v_m1*obj.sqrtA*obj.sqrtATm1q;
                                end
                            end

                        case 3    %% Solution 720.03

                            obj.s={
                                'Bruggeman (1999) solution 720.03. Steady state.';
                                'Fully penetrating wells in al aquifers at r=0.';
                                'Constant but different discharges Q'
                                };  % phi=phi(obj.x)=drawdown';

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=                     v*diag(besselk(0,d_))*v_m1*          obj.T_m1*obj.q';
                                obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*v*diag(besselk(1,d_))*v_m1*obj.sqrtA*obj.T_m1*obj.q';
                            end

                        case 11   %% Solution 720.11

                            obj.s={
                                'Bruggeman (1999) solution 720.11. Steady state.';
                                'All aquifers with open boundary. Constant drawdown of the open water level.'
                                };  % phi=phi(obj.x)=drawdown';

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=                     v*diag(besselk(0,d_))*v_m1*          K0R_m1*obj.h';
                                obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*v*diag(besselk(1,d_))*v_m1*obj.sqrtA*K0R_m1*obj.h';
                            end

                        case 12  %% Solution 720.12

                            obj.s={
                                'Bruggeman (1999) solution 720.12. Steady state.';
                                'Constant lowering of the polder level around a circular basin.';
                                'Zero drawdown at r=R.'
                            };  % phi=phi(obj.x)=drawdown';

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=           (eye(nlay)-K0R_m1*v*diag(besselk(0,d_))*v_m1)*          obj.h';
                                obj.Qx( :,ir)=-2*pi*obj.r(ir)*obj.T*K0R_m1*v*diag(besselk(1,d_))*v_m1 *obj.sqrtA*obj.h';
                            end


                        case 13  %% Solution 720.13

                            obj.s={
                                'Bruggeman (1999) solution 720.12. Steady state.';
                                'All aquifers with open boundary and zero head at r=R. Constant infiltration q.'
                                };  % phi=phi(obj.x)=drawdown';

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=           (eye(nlay)-K0R_m1*v*diag(besselk(0,d_))*v_m1)*          obj.ATm1q;
                                obj.Qx( :,ir)=-2*pi*obj.r(ir)*obj.T*K0R_m1*v*diag(besselk(1,d_))*v_m1 *obj.sqrtA*obj.ATm1q;
                            end

                        case 14    %% Solution 720.14

                            obj.s={
                                'Bruggeman (1999) solution 720.14. Steady state.';
                                'Fully penetrating circular well screens with unilateral discharges q in all';
                                'aquifers at r=R.';
                                };  % phi=phi(obj.x)=drawdown';

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=                     v*diag(besselk(0,d_))*v_m1*          K1R_m1*obj.sqrtATm1Q;
                                obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*v*diag(besselk(1,d_))*v_m1*obj.sqrtA*K1R_m1*obj.sqrtATm1Q;
                            end

                        case 15   %% Solution 720.15

                            obj.s={
                                'Bruggeman (1999) solution 720.15. Steady state.';
                                'All aquifers with entrance resistance w to open water with a constant drawdown h at r=R.';
                                'aquifers at r=R.';
                                };  % phi=phi(obj.x)=drawdown';

                            F_m1=(obj.H_m1*obj.sqrtA*K1R+K0R)^(-1);
                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=                     v*diag(besselk(0,d_))*v_m1*          F_m1*obj.h';
                                obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*v*diag(besselk(1,d_))*v_m1*obj.sqrtA*F_m1*obj.h';
                            end

                        case 16   %% Solution 720.16

                            obj.s={
                                'Bruggeman (1999) solution 720.16. Steady state.';
                                'All aquifers with entrance resistance w to open water with zero drawdown at r=R';
                                'Constant lowering of the polder level around a circular basin.'
                                };  % phi=phi(obj.x)=drawdown';

                            F_m1=(obj.H_m1*obj.sqrtA*K1R+K0R)^(-1);
                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                   obj.Phi(:,ir)=           (eye(nlay)-v*diag(besselk(0,d_))*v_m1*          F_m1)*obj.h';
                                   obj.Qx( :,ir)=-2*pi*obj.r(ir)*obj.T*v*diag(besselk(1,d_))*v_m1*obj.sqrtA*F_m1 *obj.h';
                            end

                        case 17  %% Solution 720.17

                            obj.s={
                                'Bruggeman (1999) solution 720.17. Steady state.';
                                'All aquifers with entrance resistance w to open water with zero drawdown at r=R';
                                'Constant infiltration q.'
                                };  % phi=phi(obj.x)=drawdown';

                            F_m1=(obj.H_m1*obj.sqrtA*K1R+K0R)^(-1);
                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=           (eye(nlay)-v*diag(besselk(0,d_))*v_m1*         F_m1)*obj.ATm1q;
                                obj.Qx( :,ir)=-2*pi*obj.r(ir)*obj.T*v*diag(besselk(1,d_))*v_m1*obj.sqrtA*F_m1*obj.ATm1q;
                            end

                        case 21 %% Solution 720.21

                            obj.s={
                                'Bruggeman (1999) solution 720.21. Steady state.';
                                'Given drawdown h for r=R. All aquifers with open boundary to the surface water.';
                                'Flux=0 at r=0.'
                                };  % phi=phi(obj.x)=drawdown';

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=                      v*diag(besseli(0,d_))*v_m1*          I0R_m1*obj.h';
                                obj.Qx( :,ir)=-2*pi*obj.r(ir)*obj.T*v*diag(besseli(1,d_))*v_m1*obj.sqrtA*I0R_m1*obj.h';
                            end


                        case 22   %% Solution 720.22

                            obj.s={
                                'Bruggeman (1999) solution 720.22. Steady state.';
                                'All aquifers with open boundary to the surface water with zero head at r=R.';
                                'Flux=0 at r=0. Constant infiltration `.'
                                };  % phi=phi(obj.x)=drawdown';

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=          (eye(nlay)-v*diag(besseli(0,d_))*v_m1*          I0R_m1)*obj.ATm1q;
                                obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*v*diag(besseli(1,d_))*v_m1*obj.sqrtA*I0R_m1 *obj.ATm1q;
                            end

                        case 23  %% Solution 720.23

                            obj.s={
                                'Bruggeman (1999) solution 720.23. Steady state.';
                                'All aquifers with entrance resistance w to open water with constant drawdown.';
                                'Zero flux at r=0.'
                                };  % phi=phi(obj.x)=drawdown';

                            F_m1=(obj.H_m1*obj.sqrtA*I1R+I0R)^(-1);
                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=                      v*diag(besseli(0,d_))*v_m1*          F_m1*obj.h';
                                obj.Qx( :,ir)=-2*pi*obj.r(ir)*obj.T*v*diag(besseli(1,d_))*v_m1*obj.sqrtA*F_m1*obj.h';
                            end

                        case 24  %% Solution 720.24

                            obj.s={
                                'Bruggeman (1999) solution 720.24. Steady state.';
                                'All aquifers with entrance resistance w to open water with constant drawdown.';
                                'Zero head at r=R. Zero flux at r=0. Constant infiltration q.'
                                };  % phi=phi(obj.x)=drawdown';

                            F_m1=(obj.H_m1*obj.sqrtA*I1R+I0R)^(-1);
                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=          (eye(nlay)-v*diag(besseli(0,d_))*v_m1*          F_m1)*obj.ATm1q;
                                obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*v*diag(besseli(1,d_))*v_m1*obj.sqrtA*F_m1 *obj.ATm1q;
                            end

                        case 31  %% Solution 720.31

                            obj.s={
                                'Bruggeman (1999) solution 720.31. Steady state.';
                                'Fully penetrating wells in all aquifers at r=0.';
                                'Constant but different discharges Q. Open water boundary at r=R with zero drawdown.';
                                };  % phi=phi(obj.x)=drawdown';

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=                     (v*diag(besselk(0,d_))*v_m1          -v*diag(besseli(0,d_))*v_m1*          I0R_m1*K0R)*obj.T_m1*obj.q';
                                obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*(v*diag(besselk(1,d_))*v_m1*obj.sqrtA+v*diag(besseli(1,d_))*v_m1*obj.sqrtA*I0R_m1*K0R)*obj.T_m1*obj.q';
                            end

                        case 32   %% Solution 720.32

                            obj.s={
                                'Bruggeman (1999) solution 720.32. Steady state.';
                                'Fully penetrating wells in all aquifers at r=0.';
                                'Constant but different discharges Q. Closed boundary (flux=0) at r=R.';
                                };  % phi=phi(obj.x)=drawdown';

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=                     (v*diag(besselk(0,d_))*v_m1          +v*diag(besseli(0,d_))*v_m1*          I1R_m1*K1R)*obj.T_m1*obj.q';
                                obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*(v*diag(besselk(1,d_))*v_m1*obj.sqrtA-v*diag(besseli(1,d_))*v_m1*obj.sqrtA*I1R_m1*K1R)*obj.T_m1*obj.q';
                            end

                        case 33  %% Solution 720.33

                            obj.s={
                                'Bruggeman (1999) solution 720.33. Steady state.';
                                'Fully penetrating wells in all aquifers at r=0 with constant but different discharges Q.';
                                'Bunndaries with entrance resistance w to open water at r=R.'
                                };  % phi=phi(obj.x)=drawdown';

                            F_m1=(obj.sqrtA*I1R+obj.H*I0R)^(-1);
                            G   =(obj.sqrtA*K1R-obj.H*K0R);
                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                obj.Phi(:,ir)=                     (v*diag(besselk(0,d_))*v_m1          +v*diag(besseli(0,d_))*v_m1*          F_m1*G)*obj.T_m1*obj.q';
                                obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*(v*diag(besselk(1,d_))*v_m1*obj.sqrtA+v*diag(besseli(1,d_))*v_m1*obj.sqrtA*F_m1*G)*obj.T_m1*obj.q';
                            end

                        case 34  %% Solution 720.34

                            obj.s={
                                'Bruggeman (1999) solution 720.34. Steady state.';
                                'All aquifers with fully penetrating cylindrical well screens wih radius R';
                                'and constant but different discharges q. Zero flux at r=0.'
                                };  % phi=phi(obj.x)=drawdown';

                            obj.Phi=zeros(nlay,Nr);
                            obj.Qx  =zeros(nlay,Nr);

                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                if obj.r(ir)<=R
                                    obj.Phi(:,ir)=                      R*K0R*v*diag(besseli(0,d_))*v_m1          *obj.T_m1*obj.q';
                                    obj.Qx( :,ir)=-2*pi*obj.r(ir)*obj.T*R*K0R*v*diag(besseli(1,d_))*v_m1*obj.sqrtA*obj.T_m1*obj.q';
                                else
                                    obj.Phi(:,ir)=                     R*I0R*v*diag(besselk(0,d_))*v_m1*          obj.T_m1*obj.q';
                                    obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*R*I0R*v*diag(besselk(1,d_))*v_m1*obj.sqrtA*obj.T_m1*obj.q';
                                end
                            end

                        case 35 %% Solution 720.35

                            % This one needs careful check ! May be Bruggeman constains error.

                            obj.s={
                                'Bruggeman (1999) solution 720.35. Steady state.';
                                'Circular polder with radius R, surrounded by a polder wiht different leel, or isolated';
                                'reservoir in open water.'
                                };  % phi=phi(obj.x)=drawdown';

                            h1=obj.h; h1(2:end)=0;
                            
                            for ir=1:length(obj.r)
                                [v,d_]=eig(obj.r(ir)*obj.sqrtA); v_m1=v^(-1); d_=diag(d_);
                                if obj.r(ir)<=R
                                    obj.Phi(:,ir)=          (eye(nlay)-R*obj.sqrtA*K1R*v*diag(besseli(0,d_))*v_m1)*          h1';
                                    obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*R*obj.sqrtA*K1R*v*diag(besseli(1,d_))*v_m1 *obj.sqrtA*h1';
                                else
                                    obj.Phi(:,ir)=                     R*obj.sqrtA*I1R*v*diag(besselk(0,d_))*v_m1*          h1';
                                    obj.Qx( :,ir)=2*pi*obj.r(ir)*obj.T*R*obj.sqrtA*I1R*v*diag(besselk(1,d_))*v_m1*obj.sqrtA*h1';
                                end
                            end
                        otherwise
                            error('Unknown case %d.%d',obj.maincase,obj.subcase);
                    end
                    obj.x=obj.r;  % need obj.x for plotting x-axis (r-axis)
                otherwise
                    error('Unknown maincase %d',obj.maincase);
            end
        end

        function obj=show(obj)          
            % SHOW multilayer solution
            % obj.show      plots the multilayer solution obtrainde by
            % Produces three subplots witth
            % heads, cross section with heads and stream lines and
            % discharge
            % works for onedimensional and axially symmetric cases
            % Deals with zero dimensional cases (only vertical flow 710.01, 710.02 and 710.03 separately
            

            % Dealng with 710.01, 710.02 and 710.03 separately
            if obj.maincase == 710
                if obj.subcase<10
                    figure;
                    switch obj.subcase
                        case 1
                            plot(obj.t,obj.Phi); grid on; title(obj.s);
                            xlabel('time [d]'); ylabel('head [m]');
                        case 2
                            close(gcf);
                            disp(obj.s);
                            disp(obj.Phi);
                    end
                    return;
                end
            end
            
            % Figure with readable screensize
            scrsz=get(0,'screensize'); scrsz([3 4])=round(0.6*scrsz([3 4]));
            figure('position',scrsz);
            
            % Colors for cross section plot
            yellow = [1   1   0.5];  % soft yellow (aquifers)
            grey   = [0.8 0.8 0.8];  % soft grey   (aquitards)

            Nx=size(obj.Phi,2); Nz=size(obj.Phi,1); Iz_=sort([1:Nz 1:Nz]);

            %% First plot just head lines
            subplot(3,1,1); hold on; grid on;
            plot(obj.x,obj.Phi); grid on;
            title(obj.s); xlabel('x [m]'); ylabel('head [m]');

            %% Second plot the iso head lines and streamline in cross section
            subplot(3,1,2,'color',yellow); hold on;
            title('Cross section with heads and streamlines');

            Psi=flipud(cumsum(flipud(obj.Qx)));

            contour(obj.x,obj.Z,[zeros(1,Nx);obj.Phi(obj.Iz,:)],'b');  % heads
            contour(obj.x,obj.Z,[Psi(Iz_,:);zeros(1,Nx)],'r');   % streamlines
            xlabel('x [d]'); ylabel('elevation [m]');
            legend('headlines','streamlines');

            %% plot aquitards
            for i=1:Nz
                fill(obj.x([1 end end 1]),obj.Z([2*i-1 2*i-1 2*i 2*i]),grey,'facealpha',0.7);
            end

            axis('tight');

            %% Third plot, discharge: Qx directly and form head gradient,
            %  computed right here
            
            subplot(3,1,3); hold on
            title('Discharge directly and from head gradient');
            
            Qx2=-diff(obj.Phi,1,2).*(obj.kD'*(1./diff(obj.x,1,2))); % Qx from head gradients
            xm=0.5*(obj.x(1:end-1)+obj.x(2:end));  % coordinates for head gradients

            if obj.maincase==720, % if axially symmetric
                Qx2=(ones(size(obj.kD'))*2*pi*xm) .* Qx2;
            end

            plot(obj.x,obj.Qx ,'r+');  % discharge directly from solution
            plot(   xm,    Qx2,'bx');  % discharge from head gradient
            if obj.maincase==710
                xlabel('x [m]'); ylabel('Q [m^2/d]'); % dimension m2/d
            else
                xlabel('r [m]'); ylabel('Q [m^3/d]'); % dimension m3/d
            end
        end
    end
end
