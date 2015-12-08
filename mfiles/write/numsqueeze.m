function sx=numsqueeze(x,fmt,sz)
% sx=numsqueeze(x,fmt,sz)
% squeezes string format of num x to size sz
% sometimes needed to force numbers to stay within format limits
% as Matlab formats can not always be forced to stay within limits
% if these are less than 12 or 13.
%
% TO 100406
%
sx=sprintf(fmt,x);
trunk=length(sx)-sz;
if trunk>0
    i=findstr('e',sx);  % e format ??
    if ~isempty(i)
        sx=sx([1:i-trunk-1, i:end]);
    else
        sx=sx(1:end-trunk);
    end
end
        