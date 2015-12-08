% LAYERELEVATIONS: Script which extracts layer elevations from spreadsheet
%    made by Ruben to put them
%    in the well (screen) database for direct and generic use
%
% DEPRECATED
%
% TO 091120

xlsf='registraties.xls';
sht ='layers';

[LayLabels,LayVals]=getExcelData(xlsf,sht,'H');

I=strmatchi('filter',LayLabels);
top=NaN(size(I));
bot=NaN(size(I));
for i=1:length(I)
    K=find(LayVals(:,I(i))==1);
    top(i)=LayVals(K(  1),1);
    bot(i)=LayVals(K(end),2);
end

n=(1:length(I))';
display([n top bot]);

fprintf('%d\t%e\t%e\n',[n, top, bot]');



