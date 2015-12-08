classdef MNW1Obj < wellObj
%MNW1Obj class def for multinode well according to first USGS report on MNW
%
% there will be also MNW2Obj in the near future
%
% TO 120624, Buzios, Brazil
    properties
        %Qdes  --> use Q instead
        %QWval --> use C instead
        h        = [];    % head in the multi-node well (see method hdWell)
        LOSSTYPE = 'skin' % {SKIN | LINEAR | NONLINEAR}
        skin = 0;    % friction loss due owing to screen or formation damage. If LOSSTYPE==LINEAR
                     % then skin is variable B in equation 5.
                     
                     % Note that if the Qdes=0 of the well, then Hlim and
                     % Href have no effect. (non-extracting well)
        Hlim = 0;,   % limiting water level (min for dicharging and max for rechargeing wells)
                     % if DD is set hLim is with respect to reference elevation
                     % Notice that this can be specified per stress period
                     % to take varying well flow (infiltration and
                     % extraction) into account.
                     
                     % hRef                     
        Href = 0;    % hRef reference elevation (explanation in manual is confusing). However,
                     % hRef is only used in case we limit the flow of the
                     % well based on drawdown instead of hlim. In that case
                     % the well will be limited if h_well>Href+Hlim or
                     % h_well<Href-Hilm. If the well is not drawdown
                     % limited, Href has no effect.
                     
        DD   = false;% if DD==false, then Hlim must be set in absolute elevation.
                     % such as elevation of well or bottom of aquifer.
                     % if DD==true, then Hlim is relative to Href.
                     
        QWval = 1,   % water quality given to compute mixing with other wells
        Iqwgrp,      % water-quality group identifyier for reporting (use id for each MNW) 
        PRINT =true, % if set to non zero, print this site (well id)
        lossC = 0;   % C in CWC if LOSSTYPE is LINEAR | NONLINEAR
        lossP = 1;   % P in CWC if LOSSTYPE is NONLINEAR
        QCUT  = '';
        Qfrcmn = 10; % percentage of Q min discharge to switch off pump
        Qfrcmx = 25; % percentage of Q min potential discharge to rectivate pump
    end
    methods
        function o = MNW1Obj(varargin)
            %MNW1OBJ constructor to generate multi node well objects
            %
            % generate multi-node well objects class MNW1Obj
            % USAGE:
            %       wel = MNW1Obj(basename,sheetnm);
            %       wel = MNW1Obj(nr,x,y,z[,rw]);
            % Call is exactly equal to that of oridinary wells.
            % anything in the wells sheet not in the standard variables of
            % the ordinary well will be in UserData of the well objects.
            % This inclues all fields required or optional for the MNW
            % Wen writing MNW output, we can extract this info from the
            % UserData.
            % TO 121116

            o = o@wellObj(varargin{:}); % call wellObj constructor 

            for iw=1:length(o)
                o(iw).ITYPE = 27; % Multinode ITYPE for mt3dms
                o(iw).DD    = logical(o(iw).DD);
                o(iw).PRINT = logical(o(iw).PRINT);
            end

            if nargin==0; return; end
            
        end
        
        function write(o,fid,iPer)
        %% MNW1.write(fid,iPer) -- write this MNW to file for this stress period
        %  used in writeMNW1.m called by mf_setup
        %  TO 120629

           for iw=1:length(o)
                              
               for iCell=1:length(o(iw).idx)
                   if iCell==1,
                       % we have to write out all data pertinent per well
                       % and those pertinent for cells:
                       % L R C Q QWval Rw Skin Hlim Href lqwgrp C Q-%Cut Qfrcmn Qfrcmx Site

                       MN=' '; % First line of each MNW MN is empty.

                       % All info from the wells sheet not necessary for
                       % welObj has been read into UserData. Get this data
                       % for as far as required for MNW1
                       if isfield(o(iw).UserData,'DD'),     o(iw).DD     = o(iw).UserData.DD;     end
                       if isfield(o(iw).UserData,'skin'),   o(iw).skin   = o(iw).UserData.skin;   end
                       if isfield(o(iw).UserData,'lossC'),  o(iw).lossC  = o(iw).UserData.lossC;  end
                       if isfield(o(iw).UserData,'lossP'),  o(iw).lossP  = o(iw).UserData.lossP;  end
                       if isfield(o(iw).UserData,'Href'),   o(iw).hRef   = o(iw).UserData.Href;   end
                       if isfield(o(iw).UserData,'QWval'),  o(iw).QWval  = o(iw).UserData.QWval;  end
                       if isfield(o(iw).UserData,'Iwgrp'),  o(iw).Iqwgrp = o(iw).UserData.Iqwgrp; end
                       if isfield(o(iw).UserData,'QCUT'),   o(iw).QCUT   = o(iw).UserData.QCUT;   end
                       if isfield(o(iw).UserData,'Qfrcmn'), o(iw).Qfrcmn = o(iw).UserData.Qfrcmn; end
                       if isfield(o(iw).UserData,'Qfrcmx'); o(iw).Qfrcmx = o(iw).UserData.Qfrcmx; end

                       % Default hLim: either 1 m above the bottom of
                       % the screen of an extraction well or at ground
                       % surface for an injection well. And DD=false;
                       if isempty(o(iw).Hlim)
                           if o(iw).Q(iPer)<=0, % min level 1 m above bottom of screen
                                o(iw).Hlim = -1e16; % min(o(iw).z)+1;
                                o(iw).DD   = false;
                           else % if injecting maximum level is ground surface
                                o(iw).Hlim = 1e16; % max([o(iw).z(:); o(iw).ztop]);
                                o(iw).DD   = false;
                           end
                       end
                       
                       % skin and C and QWval might depend on stress period, in that
                       % case change to
                       % o(iw).skin(min(length(o(iw).skin),iPer)
                       if o(iw).DD, dd='DD'; else dd='  '; end

                       fprintf(fid,'%10d%10d%10d%10.0f', o(iw).LRC(iCell,:),o(iw).Q(iPer));
                       fprintf(fid,'%5s %2s', MN, dd);  
                       fprintf(fid,'%10g%10g%10g',o(iw).QWval, o(iw).rw, o(iw).skin);
                       fprintf(fid,'%10g%10g',o(iw).Hlim,o(iw).Href);
                       fprintf(fid,' %4d',o(iw).Iqwgrp);
                       if strcmpi(o(iw).LOSSTYPE,'NONLINEAR')
                            fprintf(fid,' Cp:%-10g',o(iw).lossC);
                       end
                       if ~isempty(o(iw).QCUT)
                           fprintf(fid,' Q-%%CUT: %g %g DEFAULT',o(iw).Qfrcmn,o(iw).Qfrcmx);
                       end
                       if o(iw).PRINT
                           fprintf(fid,' SITE:%-s',o(iw).name);
                       end
                       fprintf(fid,'\n');
                   else
                        % we have to write out all data pertinent per cell + Qdes default = 0
                        % L R C Q Rw Skin C
                       fprintf(fid,'%10d%10d%10d%10.0f', o(iw).LRC(iCell,:),0);
                       fprintf(fid,'%5s %2s','MN',dd);  
                       fprintf(fid,'%10g%10g%10g',o(iw).QWval, o(iw).rw, o(iw).skin);
                                              fprintf(fid,'%10g%10g',o(iw).Hlim,o(iw).Href);
                       fprintf(fid,' %4d',o(iw).Iqwgrp);
                       if strcmpi(o(iw).LOSSTYPE,'NONLINEAR')
                            fprintf(fid,' Cp:%-10g',o(iw).lossC);
                       end
                       if ~isempty(o(iw).QCUT)
                           fprintf(fid,' Q-%%CUT: %g %g DEFAULT',o(iw).Qfrcmn,o(iw).Qfrcmx);
                       end
                       if o(iw).PRINT
                           fprintf(fid,' SITE:%-s',o(iw).name);
                       end
                       fprintf(fid,'\n');

                       
                       
                       
                       if strcmpi(o(iw).LOSSTYPE,'NONLINEAR')
                           fprintf(fid,' Cp:%g',o(iw).lossC);
                       end
                       fprintf(fid,'\n');                       
                   end
               end

           end
        end
        function o = hWell(o,gr,H,B)
            %MNW1.HWELL computes the head in the well given the head in the
            % model cells and the cell-by-cell flows. The computation uses
            % the actual skin loss, either linear or non-linear.
            % This method assumes the well is vertical and the layers are
            % horizontaly isotropic
            %
            % There is only one head in a multi-node well. So we can
            % compute this head from the head in any model cell penetrated
            % by the well, using the resistance of that cell and the flow
            % to/from the well in the same cell.
            %
            % The head in the cell we get from H (read by readDat)
            % The flow we get from the budget (i.e. its  field 'MNW').
            %
            % We use the first cell of the well, i.e. mnw(iw).idx(1),
            % for these compuations.
            % TO 130610
            
            for iw = numel(o):-1:1 % for all wells
                
                if numel(H)~=numel(B) % check lengths of H and B are equal
                    error('%s: size of B must equal size of H for this function',...
                        mfilename,numel(B),numel(H));
                end
                
                % Get head and mnw flow in first cell of current mnw
                for it = numel(B):-1:1
                    iLbl = strmatchi('MNW',B(it).label);
                    if ~iLbl % check of field MNW exists
                        Q1(it,1)=0;
                    end
                    Q1(1,it) = B(it).term{iLbl}(o(iw).idx(1)); % Flow in this cell                   
                    h (1,it) = H(it).values(o(iw).idx(1));     % head in this cell
                end                
                
                % effective cell radius (see manual)
                r0 = 0.28*(sqrt( (gr.dx(o(iw).ix(1))).^2 + (gr.dy(o(iw).iy(1))).^2))/2;

                % than add skin of appropriate LOSSTYPE
                switch lower(o(iw).LOSSTYPE)
                    case 'skin'
                        Res = o(iw).skin/(2*pi);
                    case 'linear'
                        Res = o(iw).B;
                    case 'non-linear'
                        Res = o(iw).B + o(iw).C.*(Q1.^(o(iw).P-1));
                    otherwise
                        error('%s: illegal LOSSTYPE <<%s>> use one of <<''skin'',''linear'',''non-linear''>>',...
                            mfilename,o(iw).LOSSTYPE);
                end
                
                % add radial resisance between cell and well
                Res = 1/(2*pi*sum(o(iw).T))*log(r0/o(iw).rw) + Res;
                
                o(iw).h = h + Q1.*Res; %res is a scalar or a vector
            end
        end
       function writePNTSRC(o,fid,iPer)
           %% MNW1.writePNTSRC(fid,iPer) -- write this MNW to file for this stress period
           %  using in writeMNW1.m called by mf_setup
           %  TO 120629
           
           iCell=1;  % only for one cell of idx

           for iw=1:length(o)
               fmt = sprintf('%%10d%%10d%%10d%%14g%%3d%s\n',repmat('%14g',[1,size(o(iw).C,1)])); 
               fprintf(fid,fmt,...
                   o(iw).LRC(iCell,:),...
                   o(iw).C(1,iPer),...
                   o(iw).ITYPE,...
                   o(iw).C(:,iPer).');
           end
        end
    end
end