system('uname -a > test.log');

fid = fopen('test.log');

s = fscanf(fid,'%s');

if isempty(regexp(s,'86_64','once'))
    fprintf('This is a 32 bit version of OSX\n');
else
    fprintf('This is a 64 bit version of OSX\n');
end