%move particles
REV.x=[   -0.3894   -0.5323   -0.6521   -0.6521   -0.6429   -0.5369,...
          -0.4585   -0.4078   -0.3664   -0.3018   -0.2051   -0.0161,...
           0.1820    0.3203    0.3756    0.4171    0.5323    0.6244,...
           0.7627    0.7903    0.7120    0.6475    0.4724    0.3664,...
           0.3433    0.3387    0.2834    0.1866    0.0161   -0.1636   -0.3894];
REV.y=[   -0.6754   -0.5234   -0.3538   -0.1140    0.0906    0.1901,...
           0.2778    0.3480    0.5643    0.7573    0.8684    0.9561,...
           0.9444    0.9211    0.8509    0.7398    0.6111    0.5117,...
           0.3480    0.0673   -0.1316   -0.1901   -0.2427   -0.3655,...
          -0.4766   -0.6287   -0.7749   -0.8392   -0.8626   -0.8567   -0.6754];
field.a = 0.0002;
field.ux =  0.002;
field.uy =  0.002;
field.scale =0.05;

Np=30; BP1.x=REV.x( 3)*ones(Np,1); BP1.y=REV.y( 3)*ones(Np,1);
       BP2.x=REV.x(15)*ones(Np,1); BP2.y=REV.y(15)*ones(Np,1);
Np=100;

xlim=[-2 5];
ylim=[-2 5];
%%
NP=5000;
C=[];       
C.x=xlim(1)+0.25*diff(xlim)*rand(NP,1)-0.21;
C.y=ylim(1)+1.50*diff(ylim)*rand(NP,1)-3.0;

%%       
NT=1500;

figure; hold on; set(gca,'xlim',xlim,'ylim',ylim,'xgrid','on','ygrid','on');

cloud(1)  =particleCloud(REV.x ,REV.y,'r',field);
cloud(2)  =particleCloud(BP1.x ,BP1.y,'g',field);
cloud(3)  =particleCloud(BP2.x ,BP2.y,'m',field);
cloud(4)  =particleCloud(C.x   ,C.y  ,'b',field);

basename='cloud';

xlabel('x'); ylabel('y');

vidObj=VideoWriter(basename);
vidObj.open;

for it=1:NT
    ttl=sprintf('particle cloud step %4d',it);
    if it==1
        ht=title(ttl);
        
        cloud(4)=cloud(4).plot;    
        cloud(3)=cloud(3).plot;    
        cloud(2)=cloud(2).plot;    
        cloud(1)=cloud(1).draw;
    else
        set(ht,'string',ttl);
        
        if rem(it,100)==0
        % reset particles on circumference of particle
            cloud(2).x(:)=mean(cloud(2).x); cloud(2).y(:)=mean(cloud(2).y);
            cloud(3).x(:)=mean(cloud(3).x); cloud(3).y(:)=mean(cloud(3).y);
            cloud(2).show;
            cloud(3).show;
        end

        cloud(4)=cloud(4).flow;
        cloud(3)=cloud(3).flow;
        cloud(2)=cloud(2).flow;
        cloud(1)=cloud(1).flow;
        
        cloud(4)=cloud(4).brown;
        cloud(3)=cloud(3).brown;
        cloud(2)=cloud(2).brown;
    end
    writeVideo(vidObj,getframe(gcf));
end

vidObj.close;

