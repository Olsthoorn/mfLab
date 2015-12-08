function [Lat,Lon, address] = geoAddress(address)
%GEOADDRESS retrieves the Lat Lon coordinates of a legal google maps address
%
% Example:
%    [Lat,Lon] = geoAddress(address); 
%
% Hint: if the address is to be composed of pieces, put those in a struct
%
% TO 121008

URL = 'http://maps.googleapis.com/maps/api/geocode/xml?address=';

if iscell(address)
    for i=2:numel(address)
        address{1} = [address{1} ', ' address{i}];
    end
    address = address{1}; address(address==' ')='+';
elseif ischar(address)
    address(address==' ')='+';
end

URL = [URL address '&sensor=false'];

A = urlread(URL);

I = regexp(A,'<geometry>'); A=A(I(1):end);
I = regexp(A,'<lat>');      A=A(I(1):end); Lat =sscanf(A,'<lat>%f',1);
I = regexp(A,'<lng>');      A=A(I(1):end); Lon =sscanf(A,'<lng>%f',1);
