function [o, HK] = BCF2LPF(o,LAYCON)
%  [Model, HK] = Model.BCF2LPF(gr) -- converts TRAN and HY to HK and VONT to VKCB using grid
%  the VCONT in BCF become VKCB and therefore, congining beds.
% Model is an array of of class modelObj, each array element contains a model
% array variable with fields name, var, type.
% where name is the name of the variable (HK, TRAN etc),
% var the variable itself
% type the mfLab type of the variable (zlist, 3Dlay, 3Dcbd, stress etc.)
%
% see mmodelObj methods and modelObj.types
%
% TO 120810

% first get the grid contained in the array of model variables (o) using
% method grid:

gr = o.grid(); % get grid from Model arrays

iTRAN = strmatchi('TRAN' ,{o.name},'exact');
iHY   = strmatchi('HY'   ,{o.name},'exact');
iVCONT= strmatchi('VCONT',{o.name},'exact');
iC    = strmatchi('C'    ,{o.name},'exact'); % Just for NHI

% Turn VCONT or C into VKCB
if iVCONT,  % obligatory in BCF if Nlay>1
    o(iVCONT).var  = o(iVCONT).var .* gr.DZcbd;
    o(iVCONT).name ='VKCB';
    o(iVCONY).type ='3Dcbd';
elseif iC % NHI
    o(iC).var = gr.DZcbd ./ o(iC).var;
    o(iC).name = 'VKCB';
    o(iC).type = '3Dcbd';
end

if iTRAN % obligatory in BCF if if any layer is nonconvertible
    nTRAN = size(o(iTRAN).var,3);
    if nTRAN == gr.Nlay  % all layers specified in TRAN, without LAYCON given
        o(iTRAN).var   = o(iTRAN).var ./ gr.DZlay; % compute HK
        o(iTRAN).name  = 'HK';  % change name trom TRAN -> HK
        if iHY % obligatory if any of the layers are convertible
            nHY=size(o(iHY).var,3);
            % put HY in the slots starting on the top as far as HY exists
            % make o(iHY) equal to o(iHK)
            o(iTRAN).var(:,:,1:nHY) = o(iHY).var;
        end
    elseif exist('LAYCON','var') % with LAYCON telling which layers are convertible and which are not
            Inonc = find(LAYCON==0 || LAYCON==2); nNonc = length(Inonc); % non-convertible layers (HY)
            Iconv = find(LAYCON==1 || LAYCON==3); nConv = length(Iconv); % convertible layers     (TRAN)
            if nTRAN>=nNonc, % sufficient number of TRAN layers given ?
                % put them in the TRAN slots
                o(iTRAN).var(:,:,Inonc)   = o(iTRAN).var(:,:,1:nNonc) ./ gr.DZlay(:,:,Inonc);
            else
                error('%s: Insufficient layers %d in TRAN compared to LAYCON=0 or LAYCON==2',mfilename,nTRAN);
            end
            if iHY
                nHY = size(o(iHY).var,3);
                if nHY>= nConv, % sufficient HY layers given?
                    % put them in the slots as prescibed by LAYCONB
                    o(iTRAN).var(:,:,Iconv) = o(iHY).var(:,:,1:nConv);
                else
                    error('%s: Insufficient layers %d in HY compared to LAYCON=0 or LAYCON==2',mfilename,nHY);
                end
            end
    else % assume all layers are given in TRAN
        error('insufficient TRAN layers %d < gr.Nlay =  %d',nTRAN,gr.Nlay);
    end
    % finally repace HY by VK == HK
    if ~iHY, iHY=length(o)+1; end
    o(iHY) = o(iTRAN);
    o(iHY).name = 'VK';
    o(iHY).type = '3Dlay'; % was already the case
else
    warning('BCF2LFP:TRAN:insufficientLayers',...
        '%s: Can'' find TRAN in your model to convert to HK',mfilename);
end

if nargout>1
    iHK = strmatchi('HK',{o.name});
    if (iHK), HK = o(iHK).var; end
end
