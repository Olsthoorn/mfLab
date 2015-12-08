function catchment = getRech(catchment,PrPER,PERLEN,S,varargin)
%GETRECH -- computes the recharge and runoff for all catchments including
% the influx of upstream catchments
%
% USAGE: getRech(catchments,PrPER,PERLEN,S[,'waterType',type[,'lines',Nlines]])
% catments = array of cathment structs containing fields
%      name
%      upStr    % names of the upstream catchments (cell array)
%      initial loss [ m ], immediate loss after a shower
%      rf       % [ - ] runoff factor a fraction of (P+Rin)-initialLoss,
%      gf       % [1/d] groundwater discharge factor
%      S0       % [ m ]initial storage in catchment
% PrPER  = [ m ] recharge per stress period
% PERLEN = [ d ] length of stress period
% S = scenarios, whether it will be shown or not
% balanceType|bType = one of 'grw','sw','tot' triggers showing the water budget
%             of groundwater, surface water of both together
% Nlines    = is number of lines to plot (default all); Inf is allowed
%
% Rin  = sum(Rout from upstream catchments if any)
% Gin  = sum(Gin  from upstream catchments if any)
%
% Pr   = PrPer/PERLEN
% Loss = min(Rin*PERLEN + PrPER - initialLoss,initialLoss)/PERLEN
% Rout = rf  * (Rin + Pr - Loss);
% Rch  = Rin + Pr - Loss - Rout;
% Gout(t) = gf * St(t-1); Gout(1)=gr*S0;
% dSdt = Gin+Rch-Gout;  dS = dSdt * dt
% swBal= Pr+Rin-Loss-Rout-Rch;
% gwBal= Gin+Rch-Gout-dSdt;
% alBal= Pr+Rin+Gin-Loss-Rout-Gout=dSdt
%
% The function recursively traveres the cathments tree. It runs through its
% upstream catchments before comptuing its own water budget terms, runoff
% and and recharge, which it will store in its UserData struct.
% After running this function, every catchment knows water budget terms
% stored in its UserDat under appropriate field names.
%
% This routine is short but advanced as recursive functions tend to be.
% Make sure the catchments have the required fields in their UserData.
%
% TO 12108

balanceType = getProp(varargin,{'water','balance'},'');
balanceType = lower(balanceType(1));
NLines    = getProp(varargin',{'N','lines'},numel(PERLEN));
NLines    = min(numel(PERLEN),NLines);

for ic=1:numel(catchment)
    % run over its upstream catchments first
    upstr = catchment(ic).UserData.upstr;
    Rin= zeros(size(PERLEN));
    Gin= zeros(size(PERLEN));
    if ~isempty(upstr)        
        Iu = strmatchi(upstr,{catchment.name},'exact');
        if Iu(1)
            for i = Iu
                if ~catchment(i).UserData.visited
                    catchment(i) = getRech(catchment(i),PrPER,PERLEN);
                end
                Aratio = catchment(i).area/catchment(ic).area;
                Rin = Rin + Aratio * catchment(i).UserData.Rout;                    
                Gin = Gin + Aratio * catchment(i).UserData.Gout;
            end
        end
    end
    
    %% at this point we are at a leaf of the catchment tree and must make up and
    %  store the budget terms of this catchment.
    % This will be done in terms of m/d, as rates.
    % So In = Out + Storage

    S0    = catchment(ic).UserData.S0;
    ILoss = catchment(ic).UserData.initialLoss; % per stress period
    rf    = catchment(ic).UserData.rf;
    gf    = catchment(ic).UserData.gf;
    
    Loss  = min(PrPER + Rin.*PERLEN , ILoss)./PERLEN; % PrPER and ILoss are per stress period
    Pr    = PrPER./PERLEN;
    Rout  = (Pr + Rin - Loss) * rf;       % runoff of surface water
    Rch   =  Pr + Rin - Loss - Rout; % not Gin is added to Recharge (in fact what adds to storage)!
   
    % Compute the storage in the catchment from the inputs and from that
    % compute the groundwater outflow
    St     = NaN(size(Pr));
    Gout   = NaN(size(Pr));
    StPrev  = S0;
    for i=1:numel(Pr)
        Gout(i)= gf*(StPrev+0.5*Rch(i));
        St(i)  = StPrev + Rch(i) + Gin(i) - Gout(i);
        StPrev = St(i);
    end
    dSdt = diff([S0;St]);
       
    catchment(ic).UserData.Pr         = Pr;
    catchment(ic).UserData.Gin        = Gin;
    catchment(ic).UserData.Rout       = Rout;
    catchment(ic).UserData.Rin        = Rin;
    catchment(ic).UserData.Loss       = Loss;
    catchment(ic).UserData.Rch        = Rch;
    catchment(ic).UserData.Gout       = Gout;
    catchment(ic).UserData.St         = St;
    catchment(ic).UserData.dSdt       = dSdt;
    catchment(ic).UserData.alBalance  = Pr+Rin-Loss-Rout+Gin-Gout-dSdt;
    catchment(ic).UserData.swBalance  = Pr+Rin-Loss-Rout-Rch;
    catchment(ic).UserData.gwBalance  = Rch+Gin-Gout-dSdt;
    catchment(ic).UserData.visited    = true;    
end

if isempty(balanceType)
    return;
end

%% Print water balance for verification'
if S.Control
    switch balanceType
        case 'g'
            gwNames = {'Rch','Gin','Gout','dSdt','gwBalance'};
            showBalance(catchment,gwNames,NLines);
        case 's'
            swNames = {'Pr','Rin','Loss','Rout','Rch','swBalance'};
            showBalance(catchment,swNames,NLines);
        case 'a'
            alNames = {'Pr','Rin','Loss','Rout','Gin','Gout','dSdt','alBalance'};
            showBalance(catchment,alNames,NLines);
        otherwise
            error('Don''t understand balance type <<%s>> use ''g'', ''s'' or ''all''',balanceType);
    end
end

function showBalance(catchment,fldNams,Nlines)
    %SHOWBALANCE -- shows water budget terms for groundwater, surface water
    % or both
    %
    % USAGE: showBalacne(fldNames,Nlines)
    %
    % TO140207
    
    for ic=1:numel(catchment)
        fprintf('%s of cathment(%d), name= ''%s''\n',fldNams{end},ic,catchment(ic).name);
        fprintf('%6s[mm/d]',fldNams{:}); fprintf('\n');
        for iLine=1:Nlines
            for iNm=1:numel(fldNams)
                fprintf('%12.3f',catchment(ic).UserData.(fldNams{iNm})(iLine) * 1000);
            end
            fprintf('\n');
        end
        fprintf('\n');
    end

