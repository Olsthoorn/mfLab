%ANALYZE analyzes layers of NHI
%
% NHI is Dutch National Hydrologic Instrumetn (www.NHI.nu)
% REGIS is regional GIS of subsurface of the Netherlands (www.dinoloket.nl)
% AGV is Amstel Gooi and Vecht waterboard model
%
% TO 120401

clayers(xm,ym,AHN,'AHN');
clayers(xm,ym,ANIFCT,'ANIFCT');
clayers(xm,ym,HK,'HK');
clayers(xm,ym,RCH,'RCH');
clayers(xm,ym,THK,'THK');
clayers(xm,ym,KD,'kD');
clayers(xm,ym,VCONT,'VCONT');
clayers(xm,ym,STRTHD,'STRTHD');
clayers(xm,ym,Z(:,:,2:end),'Z');

