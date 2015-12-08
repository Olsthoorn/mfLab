function o = monotonic(o,varargin)
%%LINEOBJ/MONOTONIC makes line object monotonically increasing or
%%decreasing
%
% USAGE o = lineObj.monotonic('up'|'down')
%
% TO 131217

[dir, ~] =getProp(varargin,'dir','down');

upward = strmatchi(dir,'up');

for io = numel(o):-1:1

    zm = [o(io).P.zm];
    sm = o(io).sCell;
    
    plot(sm,zm,mf_color(io),'lineWidth',1);
   
    zmOk = zm;
    smOk = sm;
    N = numel(zmOk);
    while true
        if upward
            keep = [true diff(zmOk)>0];
        else
            keep = [true diff(zmOk)<0];
        end
        Nnew = sum(keep);
        
        zmOk = zmOk(keep);
        smOk = smOk(keep);
        
        if Nnew==N
            break;
        else
            if Nnew==0
                error('line cannot be made monotonically up or down');
            else
                N=Nnew;
            end
        end
    end
    
   
   zmNew = interp1(smOk,zmOk,sm);
   dz    = zmNew(:) - zm(:);
   
   for i=1:numel(sm)
       o(io).P(i).z = o(io).P(i).z + dz(i);
   end
   
   plot(sm,[o(io).P.zm],mf_color(io),'marker','^','markerfaceColor','r');
   
end
   
   
    
    
    
   