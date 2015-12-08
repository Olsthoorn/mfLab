function RIV = setRIVdepth(RIV,Z,LAYCBD)
%SETRIVDEPTH sets layer number equal to the layer in which the RIV bottom resides.
%
% Example:
%    RIV = setRIVdepth(RIV,Z) 
%
% TO 120430

iLAY=2; iROW=3; iCOL=4; %iSP=1; 

if ~all(LAYCBD(1:end-1)==1),
    error(['This routine cannot handle situations if not every aquifer has a confining bed below it\n.',...
        'Consider adapting it using gr.LAYCBD accordingly.']);
end

for i=1:size(RIV,1)
    z = Z(RIV(i,iROW),RIV(i,iCOL),:);
    iZ   = find(z>RIV(i,end),1,'last');
    if ~isempty(iZ)
        RIV(i,iLAY) =floor((iZ+1)/2);  % assuming all aquifers have a confining bed below 
    end
end


