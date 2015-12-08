function ZETA=interpNaNs(ZETA)
% ZETA=interpNaNs(ZETA)
% interpolates NaNs based on average of non NaN neigbors
% continues until all NaN's have been removed

fprintf('running: interpNaN --> Interpolating ZETA to get rid of NaNs ...\n');

[NROW,NCOL,NZ]=size(ZETA);

for i=1:NZ
    fprintf('   zeta layer %d\n',i);
    zeta=ZETA(:,:,i);
    if all(isnan(zeta))
        error('Can''t interpolate matrix layer %d with only NaNs',i);
    end

    neighbor=NaN(NROW,NCOL,8);  % neighbors put in third dimension

    k=0;
    while 1
        k=k+1;
        neighbor(:,:,:)=NaN;
        neighbor(2:end  ,2:end  ,1)=zeta(1:end-1,1:end-1);
        neighbor(2:end  , :     ,2)=zeta(1:end-1, :     );
        neighbor(2:end  ,1:end-1,3)=zeta(1:end-1,2:end  );
        neighbor( :     ,1:end-1,4)=zeta( :     ,2:end  );
        neighbor(1:end-1,1:end-1,5)=zeta(2:end  ,2:end  );
        neighbor(1:end-1, :     ,6)=zeta(2:end  , :     );
        neighbor(1:end-1,2:end  ,7)=zeta(2:end  ,1:end-1);
        neighbor( :     ,2:end  ,8)=zeta( :     ,1:end-1);

        A=neighbor; B=zeros(size(A));
        A( isnan(neighbor))=0;
        B(~isnan(neighbor))=1;
        M=sum(A,3)./sum(B,3);

        zeta(any(isnan(neighbor),3))=M(any(isnan(neighbor),3));
        fprintf('Iteration %d: %d cells are still NaN\n',k,sum(isnan(zeta(:))));
        if all(~isnan(zeta(:))),
            ZETA(:,:,i)=zeta;
            break;
        end
    end
end


