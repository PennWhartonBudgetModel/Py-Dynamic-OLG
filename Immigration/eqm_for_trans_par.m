function [] = eqm_for_trans_par()

jobdir = 'Testing';
load(fullfile(jobdir, 'imm_polparams_1.mat'), 'pop_trans')

load('params.mat')
load('polparams_1.mat')


T_life   = T;
T_model  = Tss;

Dynamic.assets  = zeros(1,T_model);
Dynamic.beqs    = zeros(1,T_model);
Dynamic.labeffs = zeros(1,T_model);
Dynamic.labs    = zeros(1,T_model);
Dynamic.lfprs   = zeros(1,T_model);
Dynamic.pits    = zeros(1,T_model);
Dynamic.ssts    = zeros(1,T_model);
Dynamic.bens    = zeros(1,T_model);


for startyear = (-T_life+1):(T_model-1)
    
    fprintf('Cohort %+03d\n', startyear);
    
    T_shift  = max(+startyear, 0);
    T_active = min(startyear+T_life, T_model) - T_shift;
    
    for idem = 1:ndem
        
        Cohort = solve_cohort(startyear, idem);
        
        Dynamic.assets (T_shift+(1:T_active)) = Dynamic.assets (T_shift+(1:T_active)) + Cohort.assets ;
        Dynamic.beqs   (T_shift+(1:T_active)) = Dynamic.beqs   (T_shift+(1:T_active)) + Cohort.beqs   ;
        Dynamic.labeffs(T_shift+(1:T_active)) = Dynamic.labeffs(T_shift+(1:T_active)) + Cohort.labeffs;
        Dynamic.labs   (T_shift+(1:T_active)) = Dynamic.labs   (T_shift+(1:T_active)) + Cohort.labs   ;
        Dynamic.lfprs  (T_shift+(1:T_active)) = Dynamic.lfprs  (T_shift+(1:T_active)) + Cohort.lfprs  ;
        Dynamic.pits   (T_shift+(1:T_active)) = Dynamic.pits   (T_shift+(1:T_active)) + Cohort.pits   ;
        Dynamic.ssts   (T_shift+(1:T_active)) = Dynamic.ssts   (T_shift+(1:T_active)) + Cohort.ssts   ;
        Dynamic.bens   (T_shift+(1:T_active)) = Dynamic.bens   (T_shift+(1:T_active)) + Cohort.bens   ;
        
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