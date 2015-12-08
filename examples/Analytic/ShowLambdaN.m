%SHOWLAMBAN -- shows  1./diag(D) = lambda decoupled for all solutiions of
%seris Bruggeman 710.
%
% TO 131120
basename = 'ShowLambdaN';

scenario = 1;

Brug = 'Brug710';
Brug = 'Brug720';
            

switch scenario
    case 1
        [THdr,T] = getExcelData(basename,'T','Hor','fill',false);
        [cHdr,C] = getExcelData(basename,'C','Hor','fill',false);

        Lambda = NaN(size(T));

        for i=1:size(T,2)
            A = sysmat(T(~isnan(T(:,i)),i),C(~isnan(C(:,i)),i));
            [~,D] = eig(A);
            lambda = sqrt(1./diag(D));
            Lambda(1:numel(lambda),i) = lambda;
        end
        
        sumLambda= Lambda;
        sumLambda(isnan(sumLambda))=0; 
        sumLambda=sum(sumLambda,1);
    case 2     
        if ~exist(Brug,'var')
            BruggemanMultilayer;
        end
            eval(['B = ' Brug]);

            Lambda= NaN(20,numel(B));
            for i=1:numel(B)
                [~,D] = eig(B(i).A);
                lambda = 1./diag(D);
                Lambda(1:numel(lambda),i) = lambda;
            end
            Lambda = Lambda(~all(isnan(Lambda),2),:);
end

for i=1:size(Lambda,1)
    fprintf(' %10.0f',Lambda(i,:)); fprintf('\n');
end
