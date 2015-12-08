function MS = sorbedMass(C,ISOTHM,Vlay,RHOB,SP1,SP2)
% MS = grid.sorbeMass(C,ISOTHM,Vlay,RHOB,SP1,SP2) computes total sorbed mass in cells IDX of grid
% C is dissolve concentration
% C, RHOB,SP1 and SP2 are wither strings or vectors contaning a subset the arrays,
% by calling this function
% sorbedMass(C(IDX),ISOTHEM,gr.Vlay(IDX),RHOB(IDX),SP1(IDX),SP2(IDX).
% The indexed inputs may also be scalars.
% IDX global indices of grid cells to be included in mass computation
% ISOTHM type of sorption
% SP1 = primary   sorption coefficient depending on ISOTHM
% SP2 = secondary sorption coefficient depending on ISOTHM
%  ISOTHM = 1: Linear sorption, SP1 = Kd, SP2 = not used
%  ISOTHM = 2: Freundlich     , SP2 = Kf, SP2 = power a
%  ISOTHM = 3: Langmuir       , SP1 = Kl, SP2 = S_saturated sites [kg/kg]
% TO 121123

    switch ISOTHM
            case 0
                MS = sum(RHOB.*Vlay.*CS);
            case 1, % Linear
                MS = sum(RHOB.*Vlay.*SP1.*C);
            case 2, % Freundlich
                MS = sum(RHOB.*Vlay.*SP1.*SP2.*C.^(SP2-1));
            case 3, % Langmuir
                MS = sum(RHOB.*Vlay.*SP1.*SP2./(1+SP1.*C.^2));                                        
            otherwise
                error('%s: Unknown ISOTHM <<%d>> must be 1..0',mfilename,ISOTHM);
    end
end
