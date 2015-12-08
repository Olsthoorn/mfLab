function hout = contourf3dX(cax,x,y,z,c,n,varargin)
%Creates a contourf plane and places it in the right position

[~,handle]=contourf(cax,z,y,c,n);
    Child=get(handle,'Children');
    for i=1:length(Child)
        A=get(Child(i),'XData');
        B=get(Child(i),'YData');
        C=x*ones(size(A));
        set(Child(i),'XData',C,'YData',B,'ZData',A,varargin{:});
    end
    set(handle,'Children',Child)
    hout=handle;
end