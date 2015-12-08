function profile = saltProfile(z,zMid,sigma,cMin,cMax)
%saltProfile computes a conc profile at elevations z using erf function for gradual trans..
%  from cMin to cMax using zMid as elevation of average concentration and
%  sigma [L] to define the trans... between the cMin and cMax
%
% USAGE: profile = saltProfile(z,zMid,sigma,cMin,cMax)
%        profile = saltProfile(z,-50,25,50,10000)
%
%  z may be array or vector or scalar
% TO 150712

if isa(z,'gridObj')
    z = z.ZM;
else
    z = z(:);
end

profile    = cMin + (cMax - cMin) * 0.5 *(1+erf( -(z - zMid)/sigma));

end

