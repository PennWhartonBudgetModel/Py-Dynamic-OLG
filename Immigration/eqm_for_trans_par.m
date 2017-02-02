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
for tt1 = 1:Tss
    fprintf('Head %2u\n', tt1);
    trans_thread_head(tt1, polno, impolno);
end

for t = 0:Tss-1
    for demtype = 1:ndem
        
        load(fullfile(jobdir, sprintf('transvars_%u_%u_head_%u.mat', demtype, t+1, polno)), ...
             'dist_1', 'dist_r', 'Kalive', 'Kdead', 'ELab', 'Lab', 'Dist', 'Fedit', 'SSrev', 'SSexp', 'Lfp', 'SS_base');
        
        KPR   (Tss-t, t+(1:min(Tss-t,T))) = KPR   (Tss-t, t+(1:min(Tss-t,T))) + Kalive (1:min(Tss-t,T)) + Kdead(1:min(Tss-t,T));
        BEQ   (Tss-t, t+(1:min(Tss-t,T))) = BEQ   (Tss-t, t+(1:min(Tss-t,T))) + Kdead  (1:min(Tss-t,T));
        ELAB  (Tss-t, t+(1:min(Tss-t,T))) = ELAB  (Tss-t, t+(1:min(Tss-t,T))) + ELab   (1:min(Tss-t,T));
        LAB   (Tss-t, t+(1:min(Tss-t,T))) = LAB   (Tss-t, t+(1:min(Tss-t,T))) + Lab    (1:min(Tss-t,T));
        DIST  (Tss-t, t+(1:min(Tss-t,T))) = DIST  (Tss-t, t+(1:min(Tss-t,T))) + Dist   (1:min(Tss-t,T));
        FEDIT (Tss-t, t+(1:min(Tss-t,T))) = FEDIT (Tss-t, t+(1:min(Tss-t,T))) + Fedit  (1:min(Tss-t,T));
        SSREV (Tss-t, t+(1:min(Tss-t,T))) = SSREV (Tss-t, t+(1:min(Tss-t,T))) + SSrev  (1:min(Tss-t,T));
        SSEXP (Tss-t, t+(1:min(Tss-t,T))) = SSEXP (Tss-t, t+(1:min(Tss-t,T))) + SSexp  (1:min(Tss-t,T));
        LFP   (Tss-t, t+(1:min(Tss-t,T))) = LFP   (Tss-t, t+(1:min(Tss-t,T))) + Lfp    (1:min(Tss-t,T));
        SSBASE(Tss-t, t+(1:min(Tss-t,T))) = SSBASE(Tss-t, t+(1:min(Tss-t,T))) + SS_base(1:min(Tss-t,T));
        
    end
end


for tt1 = 2:T
    fprintf('Tail %2u\n', tt1);
    trans_thread_tail(tt1, polno, impolno);
end

for t = 1:T-1
    for demtype = 1:ndem
        
        load(fullfile(jobdir, sprintf('transvars_%u_%u_tail_%u.mat', demtype, t+1, polno)), ...
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