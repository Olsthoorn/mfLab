%% Van Genuchten relaties for theta en K voor de Staringreeks

[parLbls,soilParam,~,soilNm] = getExcelData('ClapHornberger78SoilData','Staringreeks','H');

h   = -logspace(-1,log10(16000),40)';
psi = -h;
%h = [1 10 20 31 50 100 250 500 1000 2500 5000 10000 16000]';
theta = NaN(numel(psi),numel(soilNm));
K     = NaN(numel(psi),numel(soilNm));
Se    = NaN(numel(psi),numel(soilNm));
Se_a  = NaN(numel(psi),numel(soilNm));
KvG   = NaN(numel(psi),numel(soilNm));
K2    = NaN(numel(psi),numel(soilNm));
K3    = NaN(numel(psi),numel(soilNm));
Ss    = NaN(numel(psi),numel(soilNm));
Ss2   = NaN(numel(psi),numel(soilNm));

for i = 1:numel(parLbls)
    eval([parLbls{i} '= soilParam(:,' sprintf('%d',i), ');']);
end
mvG = 1-1./nvG;

for ib = 1:numel(soilNm)
    
    theta(:,ib) = thetar(ib) + (thetas(ib) - thetar(ib)) ./ ((1+(alpha(ib) * psi).^(nvG(ib))).^mvG(ib));
    
%    Se(:,ib)   = (theta(:,ib)-thetar(ib))/(thetas(ib)-thetar(ib));
    Se(:,ib)   = (theta(:,ib) - thetar(ib))./(thetas(ib)-thetar(ib));
    
    Se_a(:,ib) = (1+(alpha(ib)*psi).^nvG(ib)).^(-mvG(ib));

    K(:,ib)  = KsvG(ib) * ( (1+(alpha(ib)*psi).^nvG(ib)).^mvG(ib) -(alpha(ib)*psi).^(nvG(ib)-1) ).^2 ./...
                  (1+(alpha(ib)*psi).^nvG(ib)).^ (mvG(ib)*(LvG(ib)+2));
    K2(:,ib) = KsvG(ib) * (Se(:,ib).^(LvG(ib)+2)) .* ...
        ( Se(:,1).^(-1) - (Se(:,ib).^(-1./mvG(ib)) - 1) ./ (alpha(ib)*psi) ).^2;
    K3(:,ib) = KsvG(ib) * Se(:,ib).^LvG(ib)  .* ...
        ( 1 - Se(:,ib) .* (Se(:,ib).^(-1./mvG(ib)) - 1) ./ (alpha(ib)*psi) ).^2;

    KvG(:,ib) = KsvG(ib) * Se(:,ib).^LvG(ib) .* (1 - (1-Se(:,ib).^(1/mvG(ib))).^mvG(ib)).^2; 

    Ss(:,ib) = -mvG(ib).*nvG(ib).*alpha(ib).*(thetas(ib)-thetar(ib)).* ...
        (1+(alpha(ib).*psi).^nvG(ib)).^(-mvG(ib)-1).*(alpha(ib).*psi).^(nvG(ib)-1);
   ss = diff(theta(:,ib))./diff(psi);  % numerical approx
   s2 = [ss(1); ss; ss(end)];
   Ss2(:,ib) = 0.5 * (s2(1:end-1)+s2(2:end));
     
end

%% Check de beringscoefficient
u = ones(size(h));

u     = @(u) ones(size(u));
theta = @(psi) u(psi) * thetar' + u(psi) * (thetas - thetar)' .* ((1+(psi * alpha').^(u(psi) * nvG')).^(u(psi) * -mvG'));
SS    = @(psi) (u(psi) * mvG') .* (u(psi) * (thetas-thetar)') .*  ((1 + (psi * alpha').^(u(psi) * nvG')) .^ (u(psi) * (-mvG'-1))) .* ...
        (u(psi) * nvG') .* (psi * alpha').^(u(psi) * (nvG'-1)) .* (u(psi) * alpha');

defaults = {'nextPlot','add','xGrid','on','yGrid','on'};
figure;
ax1 = subplot(2,1,1,defaults{:},'xlim',[0 1]);
xlabel('theta'); ylabel('pF'); title('theta(h), all Staring soils');

ax2 = subplot(2,1,2,defaults{:},'xScale','log'); ylabel('h cm'); title('Ss, all Staring soils');
xlabel('Ss [1/cm]'); ylabel('pF');

plot(ax1, theta(-h),log10(-h));
plot(ax2, SS(-h),log10(-h),'r-o'); set(gca,'xScale','log');

deltah=1;
plot(ax2, (theta(-h) - theta(-(h-deltah)))./deltah,log10(-h),'b');

% Nu nemen we een aantal profielen in het model met randvoorwaarde h =0 aan
%% de onderzijde
D = 150; dD=15;
z = (0.5*dD:-dD:-D-0.5*dD)'; zm=0.5*(z(1:end-1)+z(2:end));  dz= abs(diff(z));
h0 = -(0:dD:D); % drukhoogte aan maaiveld
h      = bsxfun( @plus,h0,zm(1)-zm);
h(h>0) = 0;

Storage = NaN(size(h,2),numel(thetas));
for iCase = 1:size(h,2)
    psi = -h(:,iCase);
    Storage(iCase,:) = sum(bsxfun(@times,dz,SS(psi)));
end

fprintf('%s','gsd[cm]');
fprintf('\t%s',soilNm{2:end}); fprintf('\n');
for i=1:numel(h0)
    fprintf('%10.0f',h0(i));
    fprintf('\t%10.3f',Storage(i,2:end));
    fprintf('\n');
end
fprintf('%10.0f',Inf);
fprintf('\t%10.3f',thetas-thetar);
fprintf('\n');
    
fprintf('%10s','grwd[cm]');
fprintf('\t%10.0f',h0');
fprintf('\t%10.0f\n',Inf);
for ib = 2:numel(soilNm)
    fprintf('%10s',soilNm{ib});
    fprintf('\t%10.3f',Storage(:,ib));
    fprintf('\t%10.3f\n',thetas(ib)-thetar(ib));
end

%beta = alpha(ib).*h;
%[1-Se(:,ib).*(Se(:,ib).^(-1./mvG(ib))-1)./beta  1-(1-Se(:,ib).^(1/mvG(ib))).^mvG(ib)]
%[Se(:,ib).*(Se(:,ib).^(-1./mvG(ib))-1)./beta  (1-Se(:,ib).^(1/mvG(ib))).^mvG(ib)]
%[Se(:,ib).* beta.^(nvG(ib)-1)  (1-Se(:,ib).^(1/mvG(ib))).^mvG(ib)]
%[Se(:,ib).^(1./mvG(ib)).* beta.^((nvG(ib)-1)./mvG(ib))  1-Se(:,ib).^(1./mvG(ib))]
%[ beta.^((nvG(ib)-1)./mvG(ib))  Se(:,ib).^(-1./mvG(ib))-1]
%[1+beta.^((nvG(ib)-1)./mvG(ib))  Se(:,ib).^(-1./mvG(ib))]
%[(1+beta.^((nvG(ib)-1)./mvG(ib))).^(-mvG(ib))  Se(:,ib)]
%[(1+beta.^((nvG(ib)-1)./mvG(ib))).^(-mvG(ib))  Se(:,ib)]

%% Test equality
%bCor.^nvG(ib)./(1+bCor.^nvG(ib))-(bCor.^nvG(ib)./(1+bCor.^nvG(ib))).^mvG(ib)
%(1./((1+bCor.^nvG(ib)).^(2*mvG(ib)))-1./(1+bCor.^nvG(ib)))/bCor

fprintf('%13s','psi','theta','Se','Se_a','K','KvG','Ss','Ss2'); fprintf \nvG

ib = 1;
display([psi theta(:,ib) Se(:,ib) Se_a(:,ib) K(:,ib) KvG(:,ib) Ss(:,ib) Ss2(:,ib)]);

%% Plot the curves voor de bodems
defaults = {'nextPlot','add','xGrid','on','yGrid','on','yScale','log','fontSize',12,'ylim',[1,1e5]};

figure('position',[37   325   560   420]);
ax1 = axes(defaults{:},'xlim',[0,1]);
xlabel(ax1,'theta(h) [-]'); ylabel(ax1,'|h| [cm]'); title(ax1,'waterretentiecurves Staringreeks');

figure('position',[665   318   560   420]);
ax2 = axes(defaults{:},'xScale','log','xlim',[1e-9,1e3]);
xlabel(ax2,'K(h) [cm/d]'); ylabel(ax2,'|h| [cm]'); title(ax2,'doorlatendheidskarateristiek Staringreeks');

for ib=1:numel(soilNm)
    if soilNm{ib}(1)=='O',
        lineWidth=2;
        lineType = '--';
    else
        lineWidth =1;
        lineType = '-';
    end
    plot(ax1,theta(:,ib), h,[mf_color(ib) lineType],'lineWidth',lineWidth);
    plot(ax2,KvG(  :,ib), h,[mf_color(ib) lineType],'lineWidth',lineWidth);

end
legend(ax1,soilNm);
legend(ax2,soilNm);

%% Corey table

[CorLbls,CorParam,~,CorNm] = getExcelData('ClapHornberger78SoilData','usaCorey','H');

%h = [1 10 20 31 50 100 250 500 1000 2500 5000 10000 16000]';
SCor     = NaN(numel(h),numel(CorNm));
thetaCor = NaN(numel(h),numel(CorNm));
KCor     = NaN(numel(h),numel(CorNm));

for i = 1:numel(CorLbls)
    eval([CorLbls{i} '= CorParam(:,' sprintf('%d',i), ');']);
end
c = 2*bCor+3;

for ib = 1:numel(CorNm)    
    SCor(:,ib)     = (h./psiAE(ib)).^(-1/bCor(ib));
    thetaCor(:,ib) = peff(ib) .* SCor(:,ib);
    KCor(:,ib)     = KsCor(ib).*SCor(:,ib).^c(ib); 
end


fprintf('%13s','h','thetaCor','SCore','KCor'); fprintf \nvG

ib = 1;
display([h thetaCor(:,ib) SCor(:,ib) KCor(:,ib)]);

%% Plot the curves voor de bodems

figure('position',[37   325   560   420]);
ax3 = axes(defaults{:},'xlim',[0,1]);
xlabel(ax3,'theta(h) [-]'); ylabel(ax3,'|h| [cm]'); title(ax3,'waterretentiecurves Staringreeks');

figure('position',[665   318   560   420]);
ax4 = axes(defaults{:},'xScale','log','xlim',[1e-9,1e3]);
xlabel(ax4,'K(h) [cm/d]'); ylabel(ax4,'|h| [cm]'); title(ax4,'doorlatendheidskarateristiek Staringreeks');

for ib=1:numel(CorNm)
    if CorNm{ib}(1)=='O',
        lineWidth=2;
        lineType = '--';
    else
        lineWidth =1;
        lineType = '-';
    end
    plot(ax3,thetaCor(:,ib), h,[mf_color(ib) lineType],'lineWidth',lineWidth);
    plot(ax4,KCor(    :,ib), h,[mf_color(ib) lineType],'lineWidth',lineWidth);
    
    plot(ax3,interp1(h,thetaCor(:,ib),psiAE(ib)),psiAE(ib),[mf_color(ib),'o']);
    plot(ax4,interp1(h,KCor(    :,ib),psiAE(ib)),psiAE(ib),[mf_color(ib),'o']);
    

end
legend(ax3,CorNm);
legend(ax4,CorNm);

