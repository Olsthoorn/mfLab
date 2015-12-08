function obs=xy2wgs(obs,e0,n0,angle)
%XY2WGS adds x,y cordinates to obs struct
%
% Example:
%    obs=xy2wgs(obs,e0,n0,angle)
%
%  add x,y coordinates to struct with fields x,y.  e0,n0 is center and
%  angle = angle of xaxis with respect to horizontal, anti clockwise
%  positive
%
% TO 111218

global Rearth

    for i=1:length(obs)
        [obs(i).u,obs(i).v]=mf_rotate(obs.x,obs.y,0,0,angle);

        obs.N = 180/pi*obs(i).v/ Rearth                 + n0;
        obs.E = 180/pi*obs(i).u/(Rearth*cos(pi/180*n0)) + e0;
    end
end
