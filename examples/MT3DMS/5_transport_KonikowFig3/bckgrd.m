%% Plot a background in the loop at end of t(1)

load lines % contains xLake yLake xMdl yMdl xIsl1 yIsl1 xIsl2 yIsl2 xRiv yRiv

plot(xMdl,yMdl,'m','linewidth',2);
plot(xRiv,yRiv,'c','linewidth',5);
fill(xLake,yLake,'b');
fill(xIsl1,yIsl1,'g');
fill(xIsl2,yIsl2,'g');
