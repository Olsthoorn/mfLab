function D=mf_observe(gr,type,Points,varargin)
%MF_OBSERVE gets data from given obvervation points
%
% Example:
%    D=mf_observe(gr,type,Points,lbl,data,lbl,data,....)
%
% Inputs:
%    D=struct with fields {'label','time',values}
%    label   is the string you used in the call of mf_observe(....,lbl,data,...)
%    time    is time obtained form data struct
%    values  are value obtained form datastruct, one column per point given
%    gr      is gridObj, generated by gr=gridObj(xGr,yGr,zGr)
%
% when type is'Points', 'Points' are intepreted as a list of points [x y] in
%    model coordinates.
% when type is 'Cells', 'Points' are interpreted as cell indices either as 'Idx or as [L R C]
%
% data is struct obtained from readDat or readMT3D or readBud
% or any struct having field time and values(Ny,Nx,Nz)
%
% Used in:
%    mflab/examples/mt3dms/Benchmarks/OneD-Lin_noneq
%    mflab/examples/mt3dms/Benchmarks/OneD-Nonlin
%    mflab/examples/mt3dms/Benchmarks/TwoD-Diagonal
%    mflab/examples/mt3dms/Benchmarks/TwoD-Radial
%    mflab/examples/mt3dms/Benchmarks/TwoD-Uniform
%
% ToDo: replace by obsObj (TO 130428)
%
% See also: gridObj obsObj observationObj
%
% TO 120415

switch lower(type)
    case 'points'
          [ix,iy]=xyzindex(Points,gr.xGr,gr.yGr);
          Idx= cellIndex(ix,iy,ones(size(ix)),gr.size);
          
    case 'cells'
         if size(Points,2)==1
             Idx=Points;
         elseif size(Points,2)==2
             Idx= gr.Ny*(Points(:,2)-1)+Points(:,1);
         else
             Idx=(gr.Ny*gr.Nx)*(Points(:,3)-1)+gr.Ny*(Points(:,2)-1)+Points(:,1);
         end
    otherwise
        
        error('mf_observe: unknown type <<%s>>',type);
end

Ndata=size(varargin,2)/2;
if floor(Ndata)~=Ndata, error('Input arguments must come in pairs ,...,Lbl,Data,Lbl,Data,..\m'); end
 

D(Ndata,1)=struct('label','','time',[],'values',[]);

for iD=1:Ndata
    D(iD).label=varargin{2*iD-1};
    Data=varargin{2*iD};
    try
        D(iD).time=vertcat(Data.time); % when concentrations
    catch ME %#ok
        D(iD).time=vertcat(Data.totim); % when heads
    end
    
    NT=size(D(iD).time,1);
    NP=size(Points,1);
    
    D(iD).values=NaN(NT,NP);
        
    for ip=1:NP
        for it=1:NT
            D(iD).values(it,ip)=Data(it).values(Idx(ip));
        end
    end
end