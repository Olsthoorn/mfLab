function timprs=MT3D_timprs(PERnams,PERvals)
% timprs=mt3D_timprs(PERnams,PERvals)
%
% times of head output, to be used to force MT3D output
% to be synchronically with that of the heads or drawdown
%
% TO 120524

NPER=size(PERvals,1);

%% Variable from the PER sheet to compute ouput times for MT3DMS
%Column                      Variable
iPERLEN = strmatchi('PERLEN',PERnams);
iNSTP   = strmatchi('NSTP'  ,PERnams);
iTSMULT = strmatchi('TSMULT',PERnams);
iIHDDFL = strmatchi('IHDDFL',PERnams);
iHdpr   = strmatchi('Hdpr'  ,PERnams);
iDdpr   = strmatchi('Ddpr'  ,PERnams);
iHdsv   = strmatchi('Hdsv'  ,PERnams);
iDdsv   = strmatchi('Ddsv'  ,PERnams);

%% tstp frequency
tstp_frequency = PERvals(:,iIHDDFL).*any(PERvals(:,[iHdpr iDdpr iHdsv iDdsv]),2);
tstp_frequency(tstp_frequency==0)=Inf;

%% Collect:
% First compute the times within the stress period and always include the
% end of the stess period. Notice that the user may vary the size, the
% number of steps and the multiplication factor on a per stress period
% basis.
N=0;
for it=NPER:-1:1
    T{it}=getT(PERvals(it,iPERLEN),PERvals(it,iNSTP),PERvals(it,iTSMULT));
    N=N+length(T{it});
end

%% Compute timprs
% The t=0 is excluded. Compute the times by adding the local times. For
% each stress period the frequency is used that pertains to it
timprs=zeros(1,NPER); k=0; t=0;
for it=1:NPER
    I = find(rem(1:PERvals(it,iNSTP),tstp_frequency(it))==0);
        
    if isempty(I)
        I= PERvals(it,iNSTP);
    end
    
    timprs(k+(1:length(I)))=t+T{it}(I);
    k=k+length(I);

    %% Add stress period time if necessary
    if I(end)~=PERvals(it,iNSTP)
        timprs(k+1)=t+PERvals(it,iPERLEN);
        k=k+1;
    end
    t=t+PERvals(it,iPERLEN);
end

end

function T=getT(perlen,nstp,tsmult)
% T=getT(perlen,nstp,tsmult)
%
% Compute the time within a stress period given its length, the number of
% time steps and the multiplication factor
%
% TO 120524

series=cumprod(tsmult*ones(1,nstp))/tsmult;
T=cumsum(perlen*series/sum(series));

end