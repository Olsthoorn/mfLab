function o = setProps(o,varargin)
    % o = setProps(o,varargin)
    % where varargin is ...,varNm,value,varNm,value,varNm,value...
    % with ... non string args and the remainer in varNm,value,pairs.
    % When varNm matches a field of o, the field is set with the
    % corresponding value.
    % This should work for any object of struct type
    % TO 130201

    if nargin<2 || isempty(varargin), return; end

    fldnames = fieldnames(o);
    for i=1:numel(varargin)
        if ischar(varargin{i})
            % first argument of type character found
            % the rest must be [varNm,value] pairs
            for iv= i:2:numel(varargin)
                if ~ischar(varargin{iv})
                    error('%s: varargin{%d} must be the name of a variable not of class <<%s>>',...
                    mfilename,iv,class(varargin{iv}));
                end
                I = strmatchi(varargin{iv},fldnames);
                if ~I, continue; end  % skip missing or misspelled field
                % always pick the last hit I(end), which is most specific
                o(i).(fldnames{I(end)}) = varargin{iv+1};
            end
            break;
        end
    end                       
end
