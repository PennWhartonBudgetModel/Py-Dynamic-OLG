function [] = solve_trans()

jobdir = 'Testing';
load(fullfile(jobdir, 'imm_polparams.mat'), 'pop_trans')

load('params.mat')


T_life   = T;
T_model  = Tss;

% Initialize aggregates
series = {'assets', 'beqs', 'labeffs', 'labs', 'lfprs', 'pits', 'ssts', 'bens'};
for o = series, Dynamic.(o{1}) = zeros(1,T_model); end


for startyear = (-T_life+1):(T_model-1)
    
    fprintf('Cohort %+03d\n', startyear);
    
    T_shift  = max(+startyear, 0);
    T_active = min(startyear+T_life, T_model) - T_shift;
    
    for idem = 1:ndem
        Cohort = solve_cohort(startyear, idem);
        for o = series, Dynamic.(o{1})(T_shift+(1:T_active)) = Dynamic.(o{1})(T_shift+(1:T_active)) + Cohort.(o{1}); end
    end
    
end


save(fullfile(jobdir, 'trans_aggregates.mat'), '-struct', 'Dynamic');




%% Testing

trans_aggregates        = load(fullfile(jobdir  , 'trans_aggregates.mat'));
trans_aggregates_freeze = load(fullfile('Freeze', 'trans_aggregates.mat'));

fprintf('trans_aggregates\n');
valuenames = fields(trans_aggregates);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = trans_aggregates.(valuename)(:) - trans_aggregates_freeze.(valuename)(:);
    if any(delta)
        pdev = abs(nanmean(delta*2 ./ (trans_aggregates.(valuename)(:) + trans_aggregates_freeze.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');



end