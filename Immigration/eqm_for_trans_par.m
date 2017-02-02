function [] = eqm_for_trans_par()

polno   = 1;
impolno = 1;

jobdir = 'Testing';
load(fullfile(jobdir, 'imm_polparams_1.mat'), 'pop_trans')

load('SSVALS.mat', 'DEBTss')
DEBT = (0.74/0.7) * DEBTss * pop_trans' / pop_trans(1);

load('params.mat')
load('polparams_1.mat')


KPR    = zeros(T+Tss,Tss);
BEQ    = zeros(T+Tss,Tss);
ELAB   = zeros(T+Tss,Tss);
LAB    = zeros(T+Tss,Tss);
DIST   = zeros(T+Tss,Tss);
FEDIT  = zeros(T+Tss,Tss);
SSREV  = zeros(T+Tss,Tss);
SSEXP  = zeros(T+Tss,Tss);
LFP    = zeros(T+Tss,Tss);
SSBASE = zeros(T+Tss,Tss);


% Newborns during the transition
for startyear = 0:Tss-1
    
    fprintf('Head %2u\n', startyear);
    trans_thread_head(startyear, polno, impolno);
    
    for idem = 1:ndem
        
        load(fullfile(jobdir, sprintf('transvars_%u_%u_head_%u.mat', idem, startyear+1, polno)), ...
             'dist_1', 'dist_r', 'Kalive', 'Kdead', 'ELab', 'Lab', 'Dist', 'Fedit', 'SSrev', 'SSexp', 'Lfp', 'SS_base');
        
        KPR   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = KPR   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Kalive (1:min(Tss-startyear,T)) + Kdead(1:min(Tss-startyear,T));
        BEQ   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = BEQ   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Kdead  (1:min(Tss-startyear,T));
        ELAB  (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = ELAB  (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + ELab   (1:min(Tss-startyear,T));
        LAB   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = LAB   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Lab    (1:min(Tss-startyear,T));
        DIST  (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = DIST  (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Dist   (1:min(Tss-startyear,T));
        FEDIT (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = FEDIT (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Fedit  (1:min(Tss-startyear,T));
        SSREV (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = SSREV (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + SSrev  (1:min(Tss-startyear,T));
        SSEXP (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = SSEXP (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + SSexp  (1:min(Tss-startyear,T));
        LFP   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = LFP   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Lfp    (1:min(Tss-startyear,T));
        SSBASE(Tss-startyear, startyear+(1:min(Tss-startyear,T))) = SSBASE(Tss-startyear, startyear+(1:min(Tss-startyear,T))) + SS_base(1:min(Tss-startyear,T));
        
    end
    
end


for tt1 = 2:T
    fprintf('Tail %2u\n', tt1);
    trans_thread_tail(tt1, polno, impolno);
end

for t = 1:T-1
    for idem = 1:ndem
        
        load(fullfile(jobdir, sprintf('transvars_%u_%u_tail_%u.mat', idem, t+1, polno)), ...
             'dist_1', 'dist_r', 'Kalive', 'Kdead', 'ELab', 'Lab', 'Dist', 'Fedit', 'SSrev', 'SSexp', 'Lfp', 'SS_base');
        
        KPR   (Tss+t+1, 1:min(Tss,T-t)) = KPR   (Tss+t+1, 1:min(Tss,T-t)) + Kalive (1:min(Tss,T-t)) + Kdead(1:min(Tss,T-t));
        BEQ   (Tss+t+1, 1:min(Tss,T-t)) = BEQ   (Tss+t+1, 1:min(Tss,T-t)) + Kdead  (1:min(Tss,T-t));
        ELAB  (Tss+t+1, 1:min(Tss,T-t)) = ELAB  (Tss+t+1, 1:min(Tss,T-t)) + ELab   (1:min(Tss,T-t));
        LAB   (Tss+t+1, 1:min(Tss,T-t)) = LAB   (Tss+t+1, 1:min(Tss,T-t)) + Lab    (1:min(Tss,T-t));
        DIST  (Tss+t+1, 1:min(Tss,T-t)) = DIST  (Tss+t+1, 1:min(Tss,T-t)) + Dist   (1:min(Tss,T-t));
        FEDIT (Tss+t+1, 1:min(Tss,T-t)) = FEDIT (Tss+t+1, 1:min(Tss,T-t)) + Fedit  (1:min(Tss,T-t));
        SSREV (Tss+t+1, 1:min(Tss,T-t)) = SSREV (Tss+t+1, 1:min(Tss,T-t)) + SSrev  (1:min(Tss,T-t));
        SSEXP (Tss+t+1, 1:min(Tss,T-t)) = SSEXP (Tss+t+1, 1:min(Tss,T-t)) + SSexp  (1:min(Tss,T-t));
        LFP   (Tss+t+1, 1:min(Tss,T-t)) = LFP   (Tss+t+1, 1:min(Tss,T-t)) + Lfp    (1:min(Tss,T-t));
        SSBASE(Tss+t+1, 1:min(Tss,T-t)) = SSBASE(Tss+t+1, 1:min(Tss,T-t)) + SS_base(1:min(Tss,T-t));
        
    end
end


OUTPUT = A*((sum(KPR)-DEBT).^(alp)).*(sum(ELAB).^(1-alp)); %#ok<NASGU>

save(fullfile(jobdir, sprintf('trans_distvars_%u_%u.mat', polno, impolno)), ...
     'dist_1', 'dist_r', 'KPR', 'BEQ', 'ELAB', 'LAB', 'DIST', 'FEDIT', 'SSREV', 'SSEXP', 'LFP', 'SSBASE', 'OUTPUT');




%% Testing

trans_distvars_1_1   = load(fullfile(jobdir  , 'trans_distvars_1_1.mat'));
trans_distvars_1_1_0 = load(fullfile('Freeze', 'trans_distvars_1_1.mat'));

fprintf('trans_distvars_1_1\n');
valuenames = fields(trans_distvars_1_1);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = trans_distvars_1_1.(valuename)(:) - trans_distvars_1_1_0.(valuename)(:);
    if any(delta)
        pdev = abs(nanmean(delta*2 ./ (trans_distvars_1_1.(valuename)(:) + trans_distvars_1_1_0.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');



end