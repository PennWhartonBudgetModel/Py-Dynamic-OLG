function [] = eqm_for_trans_par()

polno   = 1;
impolno = 1;

jobdir = 'Testing';
load(fullfile(jobdir, 'imm_polparams_1.mat'), 'pop_trans')

load('SSVALS.mat', 'DEBTss')
DEBT = (0.74/0.7) * DEBTss * pop_trans' / pop_trans(1);

load('params.mat')
load('polparams_1.mat')


% (First dimension is redundant; summation to be taken over columns)
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


for startyear = -T+1:-1
    
    fprintf('Tail %+3d\n', startyear);
    trans_thread_tail(startyear, polno, impolno);
    
    for idem = 1:ndem
        
        load(fullfile(jobdir, sprintf('transvars_%u_%u_tail_%u.mat', idem, -startyear+1, polno)), ...
             'Kalive', 'Kdead', 'ELab', 'Lab', 'Dist', 'Fedit', 'SSrev', 'SSexp', 'Lfp', 'SS_base');
        
        % (First dimension is redundant; summation to be taken over columns)
        KPR   (Tss-startyear+1, 1:min(Tss,T+startyear)) = KPR   (Tss-startyear+1, 1:min(Tss,T+startyear)) + Kalive + Kdead;
        BEQ   (Tss-startyear+1, 1:min(Tss,T+startyear)) = BEQ   (Tss-startyear+1, 1:min(Tss,T+startyear)) + Kdead  ;
        ELAB  (Tss-startyear+1, 1:min(Tss,T+startyear)) = ELAB  (Tss-startyear+1, 1:min(Tss,T+startyear)) + ELab   ;
        LAB   (Tss-startyear+1, 1:min(Tss,T+startyear)) = LAB   (Tss-startyear+1, 1:min(Tss,T+startyear)) + Lab    ;
        DIST  (Tss-startyear+1, 1:min(Tss,T+startyear)) = DIST  (Tss-startyear+1, 1:min(Tss,T+startyear)) + Dist   ;
        FEDIT (Tss-startyear+1, 1:min(Tss,T+startyear)) = FEDIT (Tss-startyear+1, 1:min(Tss,T+startyear)) + Fedit  ;
        SSREV (Tss-startyear+1, 1:min(Tss,T+startyear)) = SSREV (Tss-startyear+1, 1:min(Tss,T+startyear)) + SSrev  ;
        SSEXP (Tss-startyear+1, 1:min(Tss,T+startyear)) = SSEXP (Tss-startyear+1, 1:min(Tss,T+startyear)) + SSexp  ;
        LFP   (Tss-startyear+1, 1:min(Tss,T+startyear)) = LFP   (Tss-startyear+1, 1:min(Tss,T+startyear)) + Lfp    ;
        SSBASE(Tss-startyear+1, 1:min(Tss,T+startyear)) = SSBASE(Tss-startyear+1, 1:min(Tss,T+startyear)) + SS_base;
        
    end
    
end


% Newborns during the transition
for startyear = 0:Tss-1
    
    fprintf('Head %+3d\n', startyear);
    trans_thread_head(startyear, polno, impolno);
    
    for idem = 1:ndem
        
        load(fullfile(jobdir, sprintf('transvars_%u_%u_head_%u.mat', idem, startyear+1, polno)), ...
             'Kalive', 'Kdead', 'ELab', 'Lab', 'Dist', 'Fedit', 'SSrev', 'SSexp', 'Lfp', 'SS_base');
        
        % (First dimension is redundant; summation to be taken over columns)
        KPR   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = KPR   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Kalive + Kdead;
        BEQ   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = BEQ   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Kdead  ;
        ELAB  (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = ELAB  (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + ELab   ;
        LAB   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = LAB   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Lab    ;
        DIST  (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = DIST  (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Dist   ;
        FEDIT (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = FEDIT (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Fedit  ;
        SSREV (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = SSREV (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + SSrev  ;
        SSEXP (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = SSEXP (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + SSexp  ;
        LFP   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) = LFP   (Tss-startyear, startyear+(1:min(Tss-startyear,T))) + Lfp    ;
        SSBASE(Tss-startyear, startyear+(1:min(Tss-startyear,T))) = SSBASE(Tss-startyear, startyear+(1:min(Tss-startyear,T))) + SS_base;
        
    end
    
end


OUTPUT = A*((sum(KPR)-DEBT).^(alp)).*(sum(ELAB).^(1-alp)); %#ok<NASGU>

save(fullfile(jobdir, sprintf('trans_distvars_%u_%u.mat', polno, impolno)), ...
     'KPR', 'BEQ', 'ELAB', 'LAB', 'DIST', 'FEDIT', 'SSREV', 'SSEXP', 'LFP', 'SSBASE', 'OUTPUT');




%% Testing

trans_distvars_1_1        = load(fullfile(jobdir  , 'trans_distvars_1_1.mat'));
trans_distvars_1_1_freeze = load(fullfile('Freeze', 'trans_distvars_1_1.mat'));

fprintf('trans_distvars_1_1\n');
valuenames = fields(trans_distvars_1_1);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = trans_distvars_1_1.(valuename)(:) - trans_distvars_1_1_freeze.(valuename)(:);
    if any(delta)
        pdev = abs(nanmean(delta*2 ./ (trans_distvars_1_1.(valuename)(:) + trans_distvars_1_1_freeze.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');



end