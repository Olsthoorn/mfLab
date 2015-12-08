function well=setWellCout(well,C,iComp)
% well=setWellCout(well,C,iComp)  % adds Cout extracted by the sceen to
% wells in well.
% well is a list of wells (may be a single well) of class wellObj
% C is the struct obtained from readMT3D (see readMT3D)
% iComp is the concentration compoment in C (may be omitted)
% if wells do not yet have their injection concentration set,
% it will be made, size well.Q, and fille with NaNs
%
% The well Q and injection conc Q.C necessary for MT3DMS and SEAWAT are set
% when the wells are generated using corresponding conlumns in the PER
% worksheet.
%
% SEE ALSO: wellObj
%
% TO 120423
%
    if nargin<3,iComp=1; end
    if nargin<2
        error('mfLab:wellCout:arg1NotAWellObj',...
            'wellCout: not enough input arguments, you need (well,C[,iComp]), i.e. Conc struct, Budget struct not needed\n');
    end
    for iw=1:numel(well)
        [nC,nt]=size(well(iw).C);
        if isempty(well(iw).C),
            well(iw).C   =NaN(nC,nt);
            well(iw).Cout=NaN(nC,nt);
        else
            well(iw).Cout=zeros(nC,nt);
        end    
        for it=1:length(C)
           well(iw).Cout(iComp,it)= sum(well(iw).fQ.*C(it).values(well(iw).idx))/sum(well(iw).fQ); %.*B(it).term{iLbl}(o.idx);
        end
    end
end
