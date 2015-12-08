function animator_loglog(action)
%ANIMATOR_LOGLOG allows picking data series on screen with mouse and moving them (on logscale).
%
% USAGE:
%    testAnimator();
%    Pumptest
%
% Suitable for hand-calbrating pumping tests on double log scale.
% The only requirement being that points to move have propertie
% tag myData
%
% See also testAnimator Pumptest
%
% testAinimator uses animator_loglog for demonstration. It shows
% how data can be dragged across a log log scaled figure.
%
% TO 120114

h=findobj('tag','myData');
    
switch action
    case 'start'
        set(gcbf,'WindowButtonMotionFcn','animator_loglog move');
        set(gcbf,'WindowButtonUpFcn',    'animator_loglog stop');
        set(gcbf,'WindowButtonDownFcn',  '');

        xyP=get(gca,'CurrentPoint');
        for i=1:length(h)
            set(h(i),'UserData',xyP(1,1:2));
        end
    case 'move'
        xyP = get(gca,'CurrentPoint');
        xyP = xyP(1,1:2);
        for i=1:length(h)
            xy0 = get(h(i)  ,'UserData');
            set(h(i),'xdata',get(h(i),'xdata')*xyP(1)/xy0(1));
            set(h(i),'ydata',get(h(i),'ydata')*xyP(2)/xy0(2));
            set(h(i),'UserData',xyP);
        end
    case 'stop'
        set(gcbf,'WindowButtonMotionFcn','');
        set(gcbf,'WindowButtonUpFcn','');
        getMatchPoint();
end
        