function o=set2meanQ(o)
    %% well=well.set2meanQ --- set Q to mean of all non NaN and non zero values
    %  TO 120512
    
    for iw=1:length(o)
        Q = o(iw).Q; Q(isnan(Q))=0;
        I= find(Q ~= 0);        
        if ~isempty(I)
            o(iw).Q(I) = sum(Q(I).*o(iw).Dt(I))/sum(o(iw).Dt(I));
        end
    end
end
