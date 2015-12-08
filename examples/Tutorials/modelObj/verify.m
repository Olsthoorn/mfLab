%verify
% display transmissivity HK, resistance and VK at given position

gr = Model.grid();

Ix = hit(gr.xGr,mean(gr.xGr));
Iy = hit(gr.yGr,mean(gr.yGr));

[kD,HK] = Model.transm(    Ix(1),Iy(1));
[c ,VK] = Model.resistance(Ix(1),Iy(1));

DZ = gr.DZ(Iy(1),Ix(1),:);

fprintf('\n\n%s\n',Model(1).description{end});

fprintf('%10s %10s %10s %10s %10s %10s\n','layer','DZ','kD','HK','c','VK');
fprintf('%10g %10g %10g %10g %10g %10g\n',[(1:gr.Nz)' XS(DZ) XS(kD)   XS(HK) XS(c) XS(VK)]'); fprintf \n
fprintf('%10s %10g %10g %10s %10g %10s\n','total:',sum(DZ,3), sum(kD,3),' ', sum(c), ' ');
fprintf('\n');
