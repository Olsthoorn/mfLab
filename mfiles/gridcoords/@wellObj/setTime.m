function o=setTime(o,time,Cin)
    % well=well.setTime(time)     --- add values for t and Dt to wells
    % well=well.setTime(time,Cin) --- add values for t and Dt and Cin to wells
    %  well is a wellObj or an array of wellObj
    %  %O 120512
    Dt_=diff(time); Nt=length(Dt_);
    for iw=1:length(o)
        if length(o(iw).Q)~=Nt
            error('wellObj:setTime:QandDTincompatible',...
                'length(well(iw).Q)=%d ~= length(Dt)=%d',length(o(iw).Q),length(Dt_));
        end
        o(iw).t=time;
        o(iw).Dt=Dt_;
    end
    if nargin==3
        if length(Cin)==Nt
            Cin = Cin(:)';
        elseif length(Cin)==1
            Cin = ones(1,Nt)*Cin;
        else
            error('%s: Cin must be sclar or of length Nt = length(time)-1 = %d, but is of length %d',...
                mfilename,Nt,length(Cin));
        end
        
        for iw=1:length(o)
            o(iw).C = Cin;
        end
    end
end
