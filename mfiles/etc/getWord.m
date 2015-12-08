function [boolean,varargin] = getWord(varargin,word)
%GETWORD looks of word is one of the values in varargin
%
% USAGE:
%    [boolean,varargin] = getWord(varargin,word)
%
% useful to make input processing robust. If the word is found in varargin,
% then boolean is set to true and the word is removed form varargin.
% If not, varargin is left untouched and boolean is set to false.
%
% SEE ALSO: getNext, getProp, getType, getWord
%
% TO 130419

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if isempty(varargin)
    boolean = false;
    return; 
end

    I = strncmp(word,varargin,numel(word));
    boolean = any(I);
    
    if any(I)
        varargin(I) = [];
    end
end
