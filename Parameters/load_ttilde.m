function [] = load_ttilde()

filename = 'RevenuesPctGDP_ToJorge.xlsx';

rev_base    = xlsread(filename,1,'B5:B28')'/100;
rev_trump   = xlsread(filename,1,'D5:D28')'/100;
rev_clinton = xlsread(filename,1,'C5:C28')'/100;
rev_ryan    = xlsread(filename,1,'E5:E28')'/100;

revenue_percent = [ rev_base; rev_trump; rev_clinton; rev_ryan]; %#ok<NASGU>


save('param_gtilde.mat','revenue_percent','-append')

end