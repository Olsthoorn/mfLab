function Khet_intersection(H,gr,khettaras)
% Here you can choose your own cross section, it will find if there is an
% intersection with a khettara, and it will also plot the location of the
% khettara.

figure; hold on;
gr.contourf(H(end).values(:,:,1),min(H(end).values(:)):2:max(H(end).values(:)));
title('Choose your crossection'); axis equal;
xlabel('x UTM [m]'); ylabel('y UTM [m]');
hb = colorbar; set(get(hb,'title'),'string','Phi [m]');

[xx,yy] = ginput;
n = 50;

% xx = 1.0e+05 *[3.631007077450461 3.575095225005760]';
% yy = 1.0e+06 *[3.493265585157602 3.485378418403801]';

L = zeros([1 length(xx)]);
s = zeros([1 length(xx)]);
xn = zeros([length(xx) 1]);
yn = zeros([length(xx) 1]);

for iL = 1:length(xx)
    if iL == 1;
        L(iL) = 0;
        s(iL) = 0;
    else
        L(iL) = sqrt((xx(iL-1)-xx(iL))^2 + (yy(iL-1)-yy(iL))^2);
        l = 0:n:L(iL);
        s(iL) = s(iL-1) + l(end);
        for iS = 1:length(l)
            xn(iL,iS) = (xx(iL) - xx(iL-1)) / numel(l) * iS + xx(iL-1);
            yn(iL,iS) = (yy(iL) - yy(iL-1)) / numel(l) * iS + yy(iL-1);
        end
    end
end
        
WaterH      = interp2(gr.Xm,gr.Ym,H(end).values(:,:,1),xn,yn);
Silt        = interp2(gr.Xm,gr.Ym,gr.Z(:,:,1),xn,yn);
Gravel      = interp2(gr.Xm,gr.Ym,gr.Z(:,:,2),xn,yn);
Limestone   = interp2(gr.Xm,gr.Ym,gr.Z(:,:,3),xn,yn);
Conglomerate= interp2(gr.Xm,gr.Ym,gr.Z(:,:,4),xn,yn);
BottomL     = interp2(gr.Xm,gr.Ym,gr.Z(:,:,5),xn,yn);

Dx = zeros([length(xx) 2*numel(0:n:round(max(L)))]);
for iD = 1:length(xx)
    if iD == 1
        continue;
    else
        Len = s(iD-1):n:(round(max(L))+s(iD-1));
        Dx(iD,:) = [Len Len(end:-1:1)];
    end
end

Silt = [Silt Gravel(:,end:-1:1)];
Gravel = [Gravel Limestone(:,end:-1:1)];
Limestone = [Limestone Conglomerate(:,end:-1:1)];
Conglomerate = [Conglomerate BottomL(:,end:-1:1)];

Silt(isnan(Silt)) = 0;
Gravel(isnan(Gravel)) = 0;
Limestone(isnan(Limestone)) = 0;
Conglomerate(isnan(Conglomerate)) = 0;

figure; hold on;
for iPl = 2:length(xx)
    [ASilt] = area(Dx(iPl,:),Silt(iPl,:));
    [AGrav] = area(Dx(iPl,:),Gravel(iPl,:));
    [ALime] = area(Dx(iPl,:),Limestone(iPl,:));
    [ACong] = area(Dx(iPl,:),Conglomerate(iPl,:));
    [Water] = plot(Dx(iPl,1:length(WaterH)),WaterH(iPl,:),'r','linewidth',2);
    
    set(ASilt,'FaceColor','b');
    set(AGrav,'FaceColor','m');
    set(ALime,'FaceColor','y');
    set(ACong,'FaceColor','g');
    alpha(.3);
end

ylim([770 max(Silt(:))]);
xlabel('Length [m]'); ylabel('DEM [m]'); title('Crossection');
legend([ASilt AGrav ALime ACong Water],'Silt','Gravel','Limestone','Conglomerate','GWL')

%% Ginput lines and slopes
SlopeLine = zeros([length(xx)-1 1]);
CLine     = zeros([length(xx)-1 1]);
yLine     = zeros([length(xx)-1 1]);

for iLx = 1:length(xx)
    if iLx == 1 
        continue;
    else
        SlopeLine(iLx-1) = (yy(iLx) - yy(iLx-1)) / (xx(iLx) - xx(iLx-1));
        CLine(iLx-1) = yy(iLx-1) - SlopeLine(iLx-1) * xx(iLx-1);
        yLine(iLx-1) = SlopeLine(iLx-1) * xx(iLx-1) + CLine(iLx-1);
        yLine(iLx) = SlopeLine(iLx-1) * xx(iLx) + CLine(iLx-1);
    end
end


x         = zeros([length(khettaras) 1]); y = zeros([length(khettaras) 1]);
SlopeKhet = zeros([length(khettaras) 1]); CKhet= zeros([length(khettaras) 1]);

for iK = 1:length(khettaras)
    for iP = 1:length(khettaras(iK).P)
        if iP == 1
            x(iK,iP) = khettaras(iK).P(iP).xm;
            y(iK,iP) = khettaras(iK).P(iP).ym;
        else
            x(iK,iP) = khettaras(iK).P(iP).xm;
            y(iK,iP) = khettaras(iK).P(iP).ym; 
            SlopeKhet(iK,iP-1) = (y(iK,iP) - y(iK,iP-1)) / (x(iK,iP) - x(iK,iP-1));
            
            if abs(SlopeKhet(iK,iP-1)) == Inf || SlopeKhet(iK,iP-1) == 0
                CKhet(iK,iP-1) = NaN;
            else
                CKhet(iK,iP-1) = y(iK,iP-1) - SlopeKhet(iK,iP-1) * x(iK,iP-1);
            end
        end
    end  
end

CKhet(CKhet == 0) = NaN;
SlopeKhet(SlopeKhet == 0) = NaN;

figure; hold on;
xlabel('Length [m]'); ylabel('DEM [m]'); title('Crossection with khettaras');
Intersect  = zeros([length(xx)-1 length(khettaras) 1]);  
xintersect = zeros([length(xx)-1 length(khettaras) 1]);  
yintersect = zeros([length(xx)-1 length(khettaras) 1]); 
zintersect = zeros([length(xx)-1 length(khettaras) 1]); 
Icx = zeros([length(xx)-1 length(khettaras) 1]); 
Icy = zeros([length(xx)-1 length(khettaras) 1]); 
for iCross = 1:length(xx)-1
    for iK = 1:length(khettaras)
        for iC = 1:length(CKhet(iCross,:))
            xintersect(iCross,iK,iC) = (CKhet(iK,iC) - CLine(iCross)) / (SlopeLine(iCross) - SlopeKhet(iK,iC));
            yintersect(iCross,iK,iC) = SlopeLine(iCross) * xintersect(iCross,iK,iC) + CLine(iCross);
            zintersect(iCross,iK,iC) = 0;
            Intersect(iCross,iK,iC)  = false;
            Icx(iCross,iK,iC) = 0;
            Icy(iCross,iK,iC) = 0;
            
            if xx(iCross) > xx(iCross+1)
                if (xintersect(iCross,iK,iC) <= xx(iCross) && xintersect(iCross,iK,iC) >= xx(iCross+1)) &&...
                        (xintersect(iCross,iK,iC) >= x(iK,iC) && xintersect(iCross,iK,iC) <= x(iK,iC+1))
                    Intersect(iCross,iK,iC) = true;
                    zintersect(iCross,iK,iC) = khettaras(iK).P(iC).zm;
                    Icx(iCross,iK,iC) = khettaras(iK).P(iC).ix;
                    Icy(iCross,iK,iC) = khettaras(iK).P(iC).iy;
                    gr.xm(Icx(iCross,iK,iC));
                    gr.ym(Icy(iCross,iK,iC));
            
                else
                    Intersect(iCross,iK,iC) = false;
                end
            else
                if (xintersect(iCross,iK,iC) >= xx(iCross) && xintersect(iCross,iK,iC) <= xx(iCross+1)) &&...
                        (xintersect(iCross,iK,iC) >= x(iK,iC) && xintersect(iCross,iK,iC) <= x(iK,iC+1))
                    Intersect(iCross,iK,iC) = true;
                    zintersect(iCross,iK,iC) = khettaras(iK).P(iC).zm;
                    Icx(iCross,iK,iC) = khettaras(iK).P(iC).ix;
                    Icy(iCross,iK,iC) = khettaras(iK).P(iC).iy;
                else
                    Intersect(iCross,iK,iC) = false;
                end
            end
        end
    end
    
    Xintersect = xintersect(iCross,(find(Intersect(iCross,:,:)>0)));
    Yintersect = yintersect(iCross,(find(Intersect(iCross,:,:)>0)));
    Zintersect = zintersect(iCross,(find(Intersect(iCross,:,:)>0)));
    ICX = Icx(iCross,(find(Intersect(iCross,:,:)>0)));
    ICY = Icy(iCross,(find(Intersect(iCross,:,:)>0)));
 
    R = zeros([length(Xintersect) 1]);
    RM= zeros([length(Xintersect) 1]);
    for iX = 1:length(Xintersect)
        R(iX) = sqrt((xx(iCross) - Xintersect(iX))^2 + (yy(iCross) - Yintersect(iX))^2) + s(iCross);
        RM(iX)= sqrt((xx(iCross) - gr.xm(ICX(iX)))^2 + (yy(iCross) - gr.ym(ICY(iX)))^2) + s(iCross);
    end
   
    for iF = 1:length(xx)
        [ASilt] = area(Dx(iF,:),Silt(iF,:));
        [AGrav] = area(Dx(iF,:),Gravel(iF,:));
        [ALime] = area(Dx(iF,:),Limestone(iF,:));
        [ACong] = area(Dx(iF,:),Conglomerate(iF,:));
        [Water] = plot(Dx(iF,1:length(WaterH)),WaterH(iF,:),'r','linewidth',2);

        set(ASilt,'FaceColor','b');
        set(AGrav,'FaceColor','m');
        set(ALime,'FaceColor','y');
        set(ACong,'FaceColor','g');
        alpha(.3);
    end

    ylim([770 max(Silt(:))]);
    try
        [Khett] = plot(R,Zintersect,'ko','linewidth',3);
%         plot(RM,Zintersect,'mo');
        legend([ASilt AGrav ALime ACong Water Khett],'Silt','Gravel','Limestone','Conglomerate','GWL','Khettara');
    catch
        fprintf('no khettaras at crosssection\n');
    end
end

