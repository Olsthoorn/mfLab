function layer = BCN2Layer(o,BCN,Column)
%% layer = gr.BCN2Layer(BCN,Column)
% change the BCN list Column into a layer
% Instead of Column numbers you may use
% H for head, C for conductance Q vor Q(WEL), B for hBot(RIV)
% I for H2(CHD) J for CHDDENSOPT(CHD)
%
% TO 120605

Idx = cellIndex(BCN(:,[4 3 2]),o.size);

if isnumeric(Column)
    Column=max(1,min(size(BCN,2),Column));
else
    switch upper(Column(1))
        case 'Q', Column=5; % WEL: Q
        case 'H', Column=5; % DRN,GHB,RIV: Head / Stage
        case 'I', Column=6; % CHD: H2
        case 'J', Column=7; % CHD: CHDDENSTOPT
        case 'C', Column=6; % DRN,GHB,RIV: Conductance
        case 'B', Column=min(size(BCN,2),7); % RIV: HBOT
        otherwise
            Column=5;
    end
end

layer = o.const(0);

layer(Idx) = BCN(:,Column);


