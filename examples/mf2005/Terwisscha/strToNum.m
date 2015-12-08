function z = strToNum(str)
    if all(str==' ')
        z = NaN;
    else
        z = str2num(deblank(str));
    end
end