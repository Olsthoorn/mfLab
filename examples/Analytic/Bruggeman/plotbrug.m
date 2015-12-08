function plotbrug(x,Z,Phix,Qx,s)
%PLOTBRUG plot multilayer solution Bruggeman(1999)
%   PLotBRUG(x,Z,Phix,Qx,s) plots the multilayer solution obtrainde by
%   the script Bruggeman720series.m
%   x,Z coordianates of x and elevation of layers that is the
%   top_Atard1 bot_Atard1 bot_Aquif1 bot_Atart2 bot_Aquif2 .... bot_Aquifn
%   length(Z) must be 2*Nz+1 where Nz is number of aquifers=numer of
%   aquitards.
%   Phix and Qx computed by multilyer solution are the head and
%   the discharge in the aquifers of the multilayer system.
%   s is string for plot title.
%
%EXAMPLE
%   plotbrug(x,Z,Phix,Qx)
%
% see also Bruggeman720series.m
%
% revision TO110424

global transmissivity
persistent isOctave;
persistent isOct;

if (isempty (isOctave)); isOctave = exist ('OCTAVE_VERSION', 'builtin') ~= 0; end
if (isOctave)
  OV = strsplit (OCTAVE_VERSION, '.');
  % Octave 3.6.1 and up has multiline titles natively (fltk & qt)
  isOct = (str2num (OV{2}) >= 6);
  % Check if gnuplot = current graphics backend. It handles multi-line titles
  % if cellstr -> \n separated but there's no room in the plotbrug plots
  gnup = strcmp (graphics_toolkit(), tolower('gnuplot'));
end

if (isOctave && iscellstr(s))
  if (gnup);
    % Only show first line in multiplots
    s = s{1}; 
  end
end

yellow = [1   1   0.5];  % soft yellow
grey   = [0.8 0.8 0.8];  % soft grey

Nx=size(Phix,2); Nz=size(Phix,1); Iz=sort([1:Nz 1:Nz]);

figure; 

%% First plot just head lines
subplot(3,1,1); hold on; grid on;
plot(x,Phix); grid on;
title(s); xlabel('x [m]'); ylabel('head [m]');
 
%% Second plot the iso head lines and streamline in cross section
if (isOct) 
  subplot(3,1,2); 
else 
  subplot(3,1,2,'color',yellow); 
end; 
hold on;
Psi = real(flipud(cumsum(flipud(Qx))));

contour(x,Z,[zeros(1,Nx);real(Phix(Iz,:))],'b');  % heads
contour(x,Z,[Psi(Iz,:);zeros(1,Nx)],'r');   % streamlines
xlabel('x [m]'); ylabel('elevation [m]');
legend('headlines','streamlines');

%% plot aquitards
for i=1:Nz
    if (isOctave)
       fill(x([1 end end 1]),Z([2*i-1 2*i-1 2*i 2*i]),grey);
		else
       fill(x([1 end end 1]),Z([2*i-1 2*i-1 2*i 2*i]),grey,'facealpha',0.7);
    end
end

axis('tight');

%% Check the heads, compare with Q
subplot(3,1,3); hold on
xlabel('x [m]'); ylabel('Q [m^2/d] or [m^3/d]');
Qx2=-diff(Phix,1,2).*(transmissivity*(1./diff(x,1,2)));

xm=0.5*(x(1:end-1)+x(2:end));

if (isOctave); s = {s}; end
if ~isempty(strfind('720',s{1})), Qx2=(ones(size(transmissivity))*2*pi*xm) .* Qx2; end

plot(x ,Qx ,'+');
plot(xm,Qx2,'x');
grid on

end
