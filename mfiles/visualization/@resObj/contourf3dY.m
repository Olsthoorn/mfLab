function hout = contourf3dY(cax,x,y,z,c,n,varargin)
%Creates a contourf plane and places it in the right position

[~,handle]=contourf(cax,z,x,c,n);
    Child=get(handle,'Children');
    for i=1:length(Child)
        A=get(Child(i),'XData');
        B=get(Child(i),'YData');
        C=y*ones(size(A));
        set(Child(i),'XData',B,'YData',C,'ZData',A,varargin{:});
    end
    set(handle,'Children',Child)
    hout=handle;
end