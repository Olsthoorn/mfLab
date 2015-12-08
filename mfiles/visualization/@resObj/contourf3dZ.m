function hout = contourf3dZ(cax,x,y,z,c,n,varargin)
%Creates a contourf plane and places it in the right position

[~,handle]=contourf(cax,x,y,c,n);
    Child=get(handle,'Children');
    for i=1:length(Child)
        A=get(Child(i),'XData');
        B=get(Child(i),'YData');
        C=z*ones(size(A));
        set(Child(i),'XData',A,'YData',B,'ZData',C,varargin{:});
    end
    set(handle,'Children',Child)
    hout=handle;
end