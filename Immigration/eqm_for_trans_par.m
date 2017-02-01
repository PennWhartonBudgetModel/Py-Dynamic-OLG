function [] = eqm_for_trans_par()

%-------------------------------------
% Place main_HPC2 material here
%-------------------------------------

clear; clc;
polno=1;
impolno=1;
warning('off')
% options = optimset('MaxFunEvals',200000,'MaxIter',4000,'Algorithm','sqp','Display','off');  % run interior-point algorithm
options = optimset('Display','off','TolFun',1.e-3);  % run interior-point algorithm

% mypool = parpool(15);
%------------------------------------------------------------------- TRANSITION PATH ----------------------------------------------------------------------------------
load params.mat
pref_params = [sigma4,chi3];
rhotol=.001;
rhoeps=100;
iter=0;

rate=0;


% determining debt 
load FedDebt.mat

load SSVALS Y wg rate beq kpr DEBTss rho
% load results_1_1 beq
beqss = beq;

% beq = beq.*ones(1,Tss);
rho1 = rho.*ones(1,Tss);    % taking the steady-state value as the initial guess

jobdir = 'Testing';
load(fullfile(jobdir, 'imm_polparams_1.mat'), 'pop_trans')

% DD = (DEBTss/pop_trans(1)).*pop_trans;    % assume SS debt grows at population rate
nd = length(FederalDebtHeldbythePublic);
Tdebt = 100;    % debt to gdp ratio stabilizes at period Tdebt

load('CBODATA.mat', 'r_cbo')
r_g = 1+r_cbo;


load('gtilde.mat', 'Gtilde')


KK = (kpr/pop_trans(1)).*pop_trans';
beq = (beq/pop_trans(1)).*pop_trans';
DEBT = ((.74/.7)*DEBTss/pop_trans(1)).*pop_trans';


load params.mat
clear surv
load Surv_Probs surv_proj
surv = 1-surv_proj(1,1:T);
%             surv = 1-surv_proj(tt1,:);
load inctax_data.mat
avg_deduc = deduc_scale*avg_deduc;  % calibrating to match FIT/Y = 8% in baseline

filename = ['polparams' '_' num2str(polno) '.mat'];
load(filename)


surv(T) =0;

while rhoeps>rhotol
    tic
    iter = iter+1
    
    DD = DEBT;
    
    
    KPR = zeros(T+Tss,Tss); % (cross section for year, year)
    BEQ = zeros(T+Tss,Tss); % (cross section for year, year)
    HRS = zeros(T+Tss,Tss);
    ELAB = zeros(T+Tss,Tss);
    LAB = zeros(T+Tss,Tss);
    DIST = zeros(T+Tss,Tss);
    FEDIT = zeros(T+Tss,Tss);
    SSREV = zeros(T+Tss,Tss);
    SSEXP = zeros(T+Tss,Tss);
    LFP = zeros(T+Tss,Tss);
    SSBASE = zeros(T+Tss,Tss);


    % Newborns during the transition
%     for tt1=1:Tss
    for tt1=1:Tss
%         tt1
    %     T+tt1
        fprintf('Head %u\n', tt1);
        trans_thread_head(tt1,polno,impolno);
    end

    for tt1 =1:Tss
        for demtype = 1:ndem
            filename = ['transvars_' num2str(demtype) '_' num2str(tt1) '_head_' num2str(polno) '.mat'];
            totfile = fullfile(jobdir,filename);
            load(totfile,'dist_1','dist_r','Kalive','Kdead','ELab','Lab','Dist','Fedit','SSrev','SSexp','Lfp','SS_base');

            LAB(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1) = LAB(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1)+Lab(1:min(Tss-(tt1-1),T));
            ELAB(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1) = ELAB(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1)+ELab(1:min(Tss-(tt1-1),T));
            DIST(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1) = DIST(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1)+Dist(1:min(Tss-(tt1-1),T));

            KPR(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1) =   KPR(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1)+Kalive(1:min(Tss-(tt1-1),T))+Kdead(1:min(Tss-(tt1-1),T));
            BEQ(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1) = BEQ(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1)+Kdead(1:min(Tss-(tt1-1),T)); 
            FEDIT(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1) = FEDIT(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1)+Fedit(1:min(Tss-(tt1-1),T));
            SSREV(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1) = SSREV(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1)+SSrev(1:min(Tss-(tt1-1),T));
            SSEXP(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1) = SSEXP(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1)+SSexp(1:min(Tss-(tt1-1),T));
            LFP(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1) = LFP(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1)+Lfp(1:min(Tss-(tt1-1),T));
            SSBASE(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1) = SSBASE(Tss-(tt1-1),tt1:tt1+min(Tss-(tt1-1),T)-1)+SS_base(1:min(Tss-(tt1-1),T));
        end
    end
    % 
    % 
    % % % Alive at start of transition
    % 
    % 
%     for tt1 = 2:T
    for tt1 = 2:T
    %     tt1
        fprintf('Tail %u\n', tt1);
        trans_thread_tail(tt1,polno,impolno);
    end

    for tt1 = 2:T  

        for demtype = 1:ndem
            filename = ['transvars_' num2str(demtype) '_' num2str(tt1) '_tail_' num2str(polno) '.mat'];
            totfile = fullfile(jobdir,filename);
            load(totfile,'dist_1','dist_r','Kalive','Kdead','ELab','Lab','Dist','Fedit','SSrev','SSexp','Lfp','SS_base');

            LAB(tt1+Tss,1:min(Tss,T-(tt1-1))) = LAB(tt1+Tss,1:min(Tss,T-(tt1-1)))+Lab(1:min(Tss,T-(tt1-1)));
            ELAB(tt1+Tss,1:min(Tss,T-(tt1-1))) = ELAB(tt1+Tss,1:min(Tss,T-(tt1-1)))+ELab(1:min(Tss,T-(tt1-1)));
            DIST(tt1+Tss,1:min(Tss,T-(tt1-1))) = DIST(tt1+Tss,1:min(Tss,T-(tt1-1)))+Dist(1:min(Tss,T-(tt1-1)));
            KPR(tt1+Tss,1:min(Tss,T-(tt1-1))) = KPR(tt1+Tss,1:min(Tss,T-(tt1-1)))+Kalive(1:min(Tss,T-(tt1-1)))+Kdead(1:min(Tss,T-(tt1-1)));
            BEQ(tt1+Tss,1:min(Tss,T-(tt1-1))) = BEQ(tt1+Tss,1:min(Tss,T-(tt1-1)))+Kdead(1:min(Tss,T-(tt1-1)));
            FEDIT(tt1+Tss,1:min(Tss,T-(tt1-1))) = FEDIT(tt1+Tss,1:min(Tss,T-(tt1-1)))+Fedit(1:min(Tss,T-(tt1-1)));
            SSREV(tt1+Tss,1:min(Tss,T-(tt1-1))) = SSREV(tt1+Tss,1:min(Tss,T-(tt1-1)))+SSrev(1:min(Tss,T-(tt1-1)));
            SSEXP(tt1+Tss,1:min(Tss,T-(tt1-1))) = SSEXP(tt1+Tss,1:min(Tss,T-(tt1-1)))+SSexp(1:min(Tss,T-(tt1-1)));
            LFP(tt1+Tss,1:min(Tss,T-(tt1-1))) = LFP(tt1+Tss,1:min(Tss,T-(tt1-1)))+Lfp(1:min(Tss,T-(tt1-1)));
            SSBASE(tt1+Tss,1:min(Tss,T-(tt1-1))) = SSBASE(tt1+Tss,1:min(Tss,T-(tt1-1)))+SS_base(1:min(Tss,T-(tt1-1)));

        end

    end

    Beq = sum(BEQ)./pop_trans';  % this will give a sequence of bequests
%     Beq = [beqss Beq(1:end-1)];
    

    OUTPUT = A*((sum(KPR)-DEBT).^(alp)).*(sum(ELAB).^(1-alp));
    
%     outputeps = max(abs(OUTPUT-OUTPUTpr));
%     OUTPUTpr = OUTPUT;

    filename = ['trans_distvars_' num2str(polno) '_' num2str(impolno)  '.mat'];
    save(fullfile(jobdir, filename),'dist_1','dist_r','KPR','LAB','ELAB','OUTPUT','DIST','FEDIT','SSEXP','SSREV','LFP','SSBASE');
    
    
    
    %*********************************************************************************************************
    %*********************************************************************************************************
    %*********************DIST_TRANS***********************************************************************
    %*********************************************************************************************************
    
    
    GREV = sum(FEDIT) + sum(SSREV);
%     rate = 1 + (A*alp).*(rho1.^(alp-1)) - d;
    rate_portfolio = (rate.*sum(KPR) + (r_g).*DD)./(sum(KPR)+DD);
%     rate_g = (1+r_cbo);
    
    ssexp = sum(SSEXP);
    for t1 = 1:Tss-1
        DEBT(t1+1) = Gtilde(t1) - GREV(t1) + DEBT(t1)*r_g(t1) + ssexp(t1);
    end
%     d_y = DEBT(min(Tdebt,Tss))/OUTPUT(min(Tss,Tdebt));
%     DEBT(min(Tdebt+1,Tss):Tss) = d_y.*OUTPUT(min(Tss,Tdebt+1):Tss);   % holds the debt to gdp ratio constant past Tdebt
    
    KPRtot = sum(KPR);
%     KK = KPRtot;
%     DD = DEBT;
    ELABtot = sum(ELAB);
    
%     rhopr = ([kpr KPRtot(1:end-1)] - DEBT)./ELABtot;
    rhopr = ([kpr KPRtot(1:end-1)] - DEBT)./ELABtot;
    
    rhoeps = max(abs(rho1 - rhopr))
    beqeps = max(abs(beq-Beq))
%     debteps = max(abs(DEBT-DEBTpr))
%     DEBTpr = DEBT;
    phi_1 = .35;
    if (iter>8)&&(iter<15)
        phi_1 = .225;
    elseif iter>=15
        finaltol = rhoeps;
        rhoeps = 0;
    end
        
        
    
    rho1 = (1-phi_1).*rho1 + phi_1.*rhopr;
%     DEBTpr = DEBT;
    
    beq = Beq;
%     DEBT = DebtRatio.*OUTPUT;
%     OUTPUTpr = OUTPUT;
    
    toc
    
    rhoeps = 0;
end

% for t1 = 1:Tss-1
%     DEBT(t1+1) = GEXP(t1) - GREV(t1) + DEBT(t1)*r_g(t1) + ssexp(t1);
% end


done = 1;

% delete(mypool)




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