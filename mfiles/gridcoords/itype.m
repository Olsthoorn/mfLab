classdef itype  < int32
    % ITYPE for SSM package of MT3DMS and SEAWAT
    % used for named referencing instead of using the numbers in writeSSM
    % see Zheng & Wang (1999), p122.
    % writeSSM won't work without it.
    %
    % TO 120421
    enumeration
        CHD  (1)
        WEL  (2)
        FLUX (2)
        DRN  (3)
        RIV  (4)
        GHB  (5)
        MLS  (15)
        CCC  (-1)
        MNW  (27)
    end
end
