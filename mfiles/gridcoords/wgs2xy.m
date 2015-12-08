function obs=wgs2xy(obs,e0,n0,angle)
%WGS2XY adds x,y coordinates to obs struct with fields E,N usign e0,n0 as center
%
% Example:
%    obs=wgs2xy(obs,e0,n0,angle)
%
%    adds x,y coordinates to obs struct with fields E,N usign e0,n0 as center
%    and rotation angle angle
%
% TO 111218

global Rearth

    for i=1:length(obs)
        obs(i).v = pi/180*(obs(i).N-n0)*Rearth;
        obs(i).u = pi/180*(obs(i).E-e0)*Rearth*cos(pi/180*n0);
        if exist('angle','var')
            [obs(i).x,obs(i).y]=mf_rotate(obs(i).u,obs(i).v,0,0,-angle);
        else
            obs(i).x=[];
            obs(i).y=[];
        end
    end
