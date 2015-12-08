function gr = setAxial(gr,AXIAL)
% adding AXIAL to gridObj
% must be done be regenerating the grid
% TO 120415
    gr=gridObj(gr.xGr, gr.yGr, gr.Z, gr.LAYCBD, gr.MINDZ, AXIAL);
