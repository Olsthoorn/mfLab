% read CSV files

P = '/Users/Theo/Instituten-Groepen-Overleggen/Delfland-A4/PeilbuizenWareco/Meetreeksen/';
P = './';
pbuis = dir([P '*.csv']);
for i=numel(pbuis):-1:1
    values = csvread([P pbuis(i).name],1);
    pbuis(i).TH = [datenum(values(:,[3,2,1,4,5,6])),values(:,end)];
end

names = {pbuis.name};
%%
figure; axes('nextplot','add','xGrid','on','yGrid','on','fontsize',12);
xlabel('tijd'); ylabel('mNAP');
title(sprintf('Verloop peilbuizen WARECO raai %s HHD',sprintf('%s,',names{:})));
for i=1:numel(pbuis)
    I=pbuis(i).TH(:,1)>=datenum(2010,4,1) & pbuis(i).TH(:,1)<=datenum(2010,8,1);
    plot(pbuis(i).TH(I,1),pbuis(i).TH(I,2),mf_color(i,'gbr'));
    
    mean(pbuis(i).TH(I,2)-(-4.80))
end
datetick;
legend(names{:});
