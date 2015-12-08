function URL=mf_GM2PNG(url,varargin)
%MF_GM2PNG retrieves a GM image as png file by calling mf_GM2PIC
%
% Example:
%   see mg_GM2PIC for details, mf_GM2PNG is just a wrapper
%
% TO 110501

if isempty(varargin)
    if nargout==0
            mf_GM2PIC(url);
    else
        URL=mf_GM2PIC(url);
    end
else
    if nargout==0
            mf_GM2PIC(url,varargin);
    else
        URL=mf_GM2PIC(url,varargin);
    end
end