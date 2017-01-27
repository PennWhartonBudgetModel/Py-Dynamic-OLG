% this function solves the steady state

clear;clc;
polno=1;
impolno=1;
warning('off')
% options = optimset('MaxFunEvals',200000,'MaxIter',4000,'Algorithm','sqp','Display','off');  % run interior-point algorithm
options = optimset('Display','off','TolFun',1.e-3);%,'Algorithm','trust-region-reflective');  % run interior-point algorithm

jobdir = ['Steady State Values'];  % gives name
if exist(jobdir, 'dir')==7          % 7 = directory.
    rmdir(jobdir, 's')                     % remove all files and that directory
end
mkdir(jobdir);

load('params.mat')%,'mpci','ndem','Trmax','Tss')
Tr = NRA_baseline;
pref_params = [sigma4,chi3];

load inctax_data.mat
avg_deduc = deduc_scale*avg_deduc;

load('Imm_Data.mat')

load('CBODATA.mat', 'r_cbo')
meanr_cbo = mean(r_cbo);


tic
iter = 1;

% solving the steady-state distribution of each
% DEBTss = 2.360429083647738;    % based on 70% debt/output
% kpr = 11.4;%3*mpci;   % initial guess.  logic here is that the model is calibrated at K/Y = 3 and mpci = Y
% DD = DEBTss;
% rho =  7; 
rhosseps = 10;
% rhosstol = .1;  % tolerance on K/L ratio (rho)
rhosstol = .0001;
% beq = 0.164394453625637;
ssloc=1;
% parpool('local',ndem);

im_scale = 1;    % scaling the immigration flow to target 13% immigrant population

load SSVALS.mat   % loads the last set of outputs as the first guess


% starting steady state equilbrium calculation
while rhosseps>rhosstol
    DD = DEBTss;  % debt used in the calculation of prices.  zero for open economy
    tic
%     for demtype = 1:ndem
    %*******************************************************************************************************************************************************************
    %*******************************************************************************************************************************************************************
    %*******************************************************************************************************************************************************************
    %*******************************************************************************************************************************************************************
    %*******************************************************************************************************************************************************************
    % 

    iteration_counter = [0,0];

    for demtype = 1:ndem
        

%         meas_pop_age = MU2(demtype,:);    % measure of population by age (for this demographic type)
%         mu3 = MU3(demtype,:);
        



        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        

        filename = ['polparams' '_' num2str(polno) '.mat'];
        load(filename)
        
        
        

        % ssloc gives 1 if beginning of transition and 2 if end of transition
        if ssloc==1
            ben = ss_benefit(:,1);
        %     ss_disc = ss_disc(:,1);
            beq = beq(1);
%             surv = 1-surv_proj(1,:);
        elseif ssloc==2
            ben = ss_benefit(:,Tss);
        %     ss_disc = ss_disc(:,Tss);
            beq = beq(end);
%             surv = 1-surv_proj(end,:);
        end

        taxmax = taxmax(1);


        wage = (A*(1-alp))*(rho^(alp));
        rate1 = 1 + (A*alp)*(rho^(alp-1)) - d;
%         rate=rate1;
        rate = (rate1*kpr + (1+meanr_cbo)*DD)/(kpr+DD);



        % Creating arrays
        V = zeros(nk,nz,nb,Tr+1);
        labopt = zeros(nk,nz,nb,Tr);
        kopt = zeros(nk,nz,nb,Tr);
        copt = zeros(nk,nz,nb,Tr);
        bopt = zeros(nk,nz,nb,Tr);
        fedincome = zeros(nk,nz,nb,Tr);
        fitax = zeros(nk,nz,nb,Tr);
        fsstax = zeros(nk,nz,nb,Tr);

        % cont_util = zeros(nk,nz,nb,T);

        Vss = zeros(nk,nb,T+1-Tr);
        laboptss = zeros(nk,nb,T-Tr);
        koptss = zeros(nk,nb,T-Tr);
        boptss = zeros(nk,nb,T-Tr);
        coptss = zeros(nk,nb,T-Tr);
        fedincomess = zeros(nk,nb,T-Tr);
        fitaxss = zeros(nk,nb,T-Tr);
        fsstaxss = zeros(nk,nb,T-Tr);
        % cont_utilss = zeros(nk,nz,nb,T-Tr);

        % bequest motive
        Vbeq = phi1.*((1+kgrid./phi2).^(1-phi3));



        surv(T) =0;
        clear argmax val
        for t1 = T-Tr:-1:1
%             t1+Tr
            totben = ben;
            EV = ((1-surv(t1+Tr))*Vbeq)*ones(1,nb); % gives continuation value dimensions [k',b']
            EV = EV + (surv(t1+Tr)*b)*squeeze(Vss(:,:,t1+1));
            for i1 = 1:nb
                for k1 = 1:nk
                    
                    income = rate.*kgrid(k1) + totben(i1);
                    fincome = (rpci/mpci)*((rate-1).*kgrid(k1)  + (1-ss_tax_cred)*totben(i1));
                    deduc = avg_deduc + coefs(1).*fincome + coefs(2).*fincome.^2 + coefs(3).*fincome.^3 + coefs(4).*fincome.^4;
                    ftax = (mpci/rpci)*limit.*((max(fincome-deduc,0)) - ((max(fincome-deduc,0)).^(-X(1)) + (X(2))).^(-1/X(1)));
                    sstax = 0;

                    eff_inc = income - ftax - sstax + beq;    % effective income after tax and bequest.  just subtract savings to get consumption.
                    [argmax, val] = fminsearch(@(x) valfuncold(x,kgrid,eff_inc,EV(:,i1),pref_params),kgrid(k1),options);
                    
                    koptss(k1,i1,t1) = argmax;
                    Vss(k1,i1,t1) = -1*val;
                    coptss(k1,i1,t1) = eff_inc-argmax;
%                     fedincomess(k1,i1,t1) = ;
                    fitaxss(k1,i1,t1) = ftax;
        %             fsstaxss(:,i1,t1) = sstax(IND1);
                end
            end
        end


        % Pre-retirement age
        
        E2 =  (1-surv(Tr))*(Vbeq*ones(1,nb));
        for j1 = 1:nz
            for j2 = 1:nz
                E2 = E2 + (surv(Tr)*b)*tr_z(j1,j2)*squeeze(Vss(:,:,1));
            end
        end
                            

        for t1 = Tr:-1:1
%             tic
%             t1
%             tax_coefs = [avg_deduc; coefs; limit; X; mpci; rpci; tau_ss; v_ss_max];
            cap_gain = (rate-1).*kgrid;
            cap_inc = rate.*kgrid;
            for i1 = 1:nb
                for j1 = 1:nz
                    eff_wage = wage*max(z(j1,t1,demtype),0);%  + soc_sec*ge(t1,Tr);
                    EV = (1-surv(t1))*(Vbeq*ones(1,nb)); % gives continuation value dimensions [k',b']
                    for j2 = 1:nz
                        EV = EV + (surv(t1)*b)*tr_z(j1,j2)*squeeze(V(:,j2,:,t1+1));%
                    end
                    
                    EV = lt(t1,Tr)*EV + eq(t1,Tr)*E2;    %
                    
                    for k1 = 1:nk
%                         [argmax] = fsolve(@(x) valzero(x,kgrid,bgrid,cap_inc(k1),cap_gain(k1),eff_wage,beq,EV,pref_params,avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max,t1,i1,nb,nk),[kgrid(k1), .2],options);
%                         argmax(1) = max(argmax(1),kgrid(1));
%                         val = valfunc(argmax,kgrid,bgrid,cap_inc(k1),cap_gain(k1),eff_wage,beq,EV,pref_params,avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max,t1,i1,nb,nk);
                        
                        [argmax,val] = fminsearch(@(x) valfunc(x,kgrid,bgrid,cap_inc(k1),cap_gain(k1),eff_wage,beq,EV,pref_params,avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max,t1,i1,nb,nk),[kgrid(k1), .2],options);
                        income = cap_inc(k1) + eff_wage*argmax(2);
                        fincome = (rpci/mpci)*(cap_gain(k1) + eff_wage*argmax(2));%  + soc_sec*ge(t1,Tr);
                        deduc = avg_deduc + coefs(1).*fincome + coefs(2).*fincome.^2 + coefs(3).*fincome.^3 + coefs(4).*fincome.^4;
                        ftax = (mpci/rpci)*limit*((max(fincome-deduc,0)) - ((max(fincome-deduc,0)).^(-X(1)) + (X(2))).^(-1/X(1)));
                        sstax = tau_ss*min(eff_wage*argmax(2),v_ss_max);
                        cons1 = income - ftax - sstax - argmax(1) + beq;

                        V(k1,j1,i1,t1) = -1*val;
                        kopt(k1,j1,i1,t1) = argmax(1);  %(nl,nk,nk,4,T-Tr,nep,nh,ng)
                        labopt(k1,j1,i1,t1) = argmax(2);
                        copt(k1,j1,i1,t1) = cons1;
                        fedincome(k1,j1,i1,t1) = fincome;
                        fitax(k1,j1,i1,t1) = ftax;
                        fsstax(k1,j1,i1,t1) = sstax;
                        bopt(k1,j1,i1,t1) =(1/t1).*(bgrid(i1)*(t1-1) + min(v_ss_max,eff_wage*argmax(2))); % could also do bprime(Il);
                        
                    end
                    

                end
            end
        end
    

        % filename = ['agent_pol' num2str(demtype) '.mat'];
        filename = ['sspol' num2str(demtype) '.mat'];
        jobdir = ['Steady State Values'];  % gives name
        totfile = fullfile(jobdir,filename);
        
        save(totfile,'V','Vss','copt','coptss','kopt','koptss','labopt','laboptss','fedincome','fedincomess','fitax','fitaxss','fsstax','fsstaxss','bopt','boptss','ben','bgrid');
    end



        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %----------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    for demtype = 1:ndem
        Kalive=zeros(1,T);
        Kdead=zeros(1,T);
        Lab = zeros(1,T);
        ELab = zeros(1,T);
        
        filename = ['sspol' num2str(demtype) '.mat'];
        jobdir = ['Steady State Values'];  % gives name
        totfile = fullfile(jobdir,filename);
        load(totfile);
        iter = 0;

        dist1ss_previous = zeros(nk,nz,nb,T,3);
        dist1ss_r_previous = zeros(nk,nb,T,3);
        
        disteps = 100;
        disttol = .0001;
        pop_prev = 0;
        pop = 1;   % initial assertion of population
        
        dist_age_previous = ones(T,1);
        
        while disteps>disttol
            dist1ss = zeros(nk,nz,nb,T,3);   % last dimension is [native, legal, illegal]
            dist1ss_r = zeros(nk,nb,T,3);
            iter = iter + 1
            iteration_counter(demtype) = iter;

            % Initiating Natives
%             k0 = 0;    % initial capital  **********CHANGE IMMIGRANT CAPITAL IF THIS IS DIFFERENT FROM ZERO**********
            
            % pgr is population growth rate of existing population.  only grows youngest cohort.
            im_flow = [pop*(pgr); im_scale*pop*imm_age(1)*legal_rate(1); im_scale*pop*imm_age(1)*illegal_rate(1)];   % using period 1 imm rate values for steady state

            loc1 = 1; loc2 = 1;    % initiating 
            for i1=1:nz
                dist1ss(loc1,i1,loc2,1,:) = squeeze(proddist_age(i1,1,:)).*(im_flow);
            end
            

            
            for t1 = 1:Tr
                im_flow = [0; im_scale*pop*imm_age(t1)*legal_rate(1); im_scale*pop*imm_age(t1)*illegal_rate(1)];   % using period 1 imm rate values for steady state
                for j2 = 1:nz
                    for i1 = 1:nk
                        for i2 = 1:nb
                            
                            point_k = max(kopt(i1,j2,i2,t1),kgrid(1));    % placing floor at lowest gridpoint.
                            loc1 = find(kgrid(1:nk-1)<=point_k,1,'last');
                            w1 = (point_k - kgrid(loc1))/(kgrid(loc1+1)-kgrid(loc1));  % amount allocated to higher gridpoint
                            w1 = min(w1,1);

                            point_b = max(bopt(i1,j2,i2,t1),bgrid(1));    % placing floor at lowest gridpoint.
                            loc2 = find(bgrid(1:nb-1)<=point_b,1,'last');    % lower gridpoint
                            w2 = (point_b - bgrid(loc2))/(bgrid(loc2+1)-bgrid(loc2));  % amount allocated to higher gridpoint
                            w2 = min(w2,1);

%                             dist_hold = squeeze(dist1ss(i1,j2,i2,t1,:)).*(1+squeeze(proddist_age(j2,t1,:)).*im_flow);
                            dist_hold = squeeze(dist1ss_previous(i1,j2,i2,t1,:)) + (eq(i1,1)*eq(i2,1))*squeeze(proddist_age(j2,t1,:)).*im_flow;

                            for j4 = 1:nz
                                
                                if t1<Tr
                                    dist1ss(loc1,j4,loc2,t1+1,:) = squeeze(dist1ss(loc1,j4,loc2,t1+1,:)) + (surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                    dist1ss(loc1+1,j4,loc2,t1+1,:) = squeeze(dist1ss(loc1+1,j4,loc2,t1+1,:)) + (surv(t1)*(1-w2)*(w1)*tr_z(j2,j4)).*dist_hold;
                                    dist1ss(loc1,j4,loc2+1,t1+1,:) = squeeze(dist1ss(loc1,j4,loc2+1,t1+1,:)) + (surv(t1)*(w2)*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                    dist1ss(loc1+1,j4,loc2+1,t1+1,:) = squeeze(dist1ss(loc1+1,j4,loc2+1,t1+1,:)) + (surv(t1)*(w2)*(w1)*tr_z(j2,j4)).*dist_hold;
                                elseif t1==Tr
                                    dist1ss_r(loc1,loc2,t1+1,:) = squeeze(dist1ss_r(loc1,loc2,t1+1,:)) + (surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                    dist1ss_r(loc1+1,loc2,t1+1,:) = squeeze(dist1ss_r(loc1+1,loc2,t1+1,:)) + (surv(t1)*(1-w2)*(w1)*tr_z(j2,j4)).*dist_hold;
                                    dist1ss_r(loc1,loc2+1,t1+1,:) = squeeze(dist1ss_r(loc1,loc2+1,t1+1,:)) + (surv(t1)*(w2)*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                    dist1ss_r(loc1+1,loc2+1,t1+1,:) = squeeze(dist1ss_r(loc1+1,loc2+1,t1+1,:)) + (surv(t1)*(w2)*(w1)*tr_z(j2,j4)).*dist_hold;
                                end
                            end
                        end
                    end
                end
            end
            
            
            

            for t1 = Tr+1:T-1
                im_flow = [0; im_scale*pop*imm_age(t1)*legal_rate(1); im_scale*pop*imm_age(t1)*illegal_rate(1)];   % using period 1 imm rate values for steady state
                for i1 = 1:nk
                    for i2 = 1:nb
                        point_k = max(koptss(i1,i2,t1-Tr),kgrid(1));    % placing floor at lowest gridpoint.
                        loc1 = find(kgrid(1:nk-1)<=point_k,1,'last');
                        w1 = (point_k - kgrid(loc1))/(kgrid(loc1+1)-kgrid(loc1));  % amount allocated to higher gridpoint
                        w1 = min(w1,1);
        
                        dist_hold = squeeze(dist1ss_r_previous(i1,i2,t1,:)) + (eq(i1,1)*eq(i2,1)).*im_flow;
                        dist1ss_r(loc1,i2,t1+1,:) = squeeze(dist1ss_r(loc1,i2,t1+1,:)) + (surv(t1)*(1-w1)).*dist_hold;
                        dist1ss_r(loc1+1,i2,t1+1,:) = squeeze(dist1ss_r(loc1+1,i2,t1+1,:)) + (surv(t1).*w1)*dist_hold;

                    end
                end
            end
            
            pop_prev = pop;
            pop = sum(sum(sum(sum(sum(dist1ss))))) + sum(sum(sum(sum(dist1ss_r))));    % summing population
            
%             dist1ss = dist1ss.*(T/dist_tot);    % normalizing the distribution so that it has size T
%             dist1ss_r = dist1ss_r.*(T/dist_tot);
           
            dist_age = [squeeze(sum(sum(sum(sum(permute(dist1ss,[1,2,3,5,4]))))))+squeeze(sum(sum(sum(permute(dist1ss_r,[1,2,4,3])))))];
            death_rate = sum(dist_age'.*(1-surv(1:T)))/pop;  % Just for debug
            % Below, convergence is relative to
            disteps = max(abs((1/dist_age(1)).*dist_age(2:end) - (1/dist_age_previous(1)).*dist_age_previous(2:end)))
            dist_age_previous = dist_age;
            dist1ss_previous = dist1ss;   % updating working-age distribution array
            dist1ss_r_previous = dist1ss_r;   % updating retiree distribution array
            
        end
            
        immigrant_population = 1-(sum(sum(sum(sum(dist1ss(:,:,:,:,1)))))+sum(sum(sum(dist1ss_r(:,:,:,1)))))/(sum(sum(sum(sum(sum(dist1ss)))))+sum(sum(sum(sum(dist1ss_r)))))
        % solving for aggregates by age
        for i3 = 1:3
            for k1 = 1:Tr
    %             node1 = 1;
                for k2 = 1:nz
                    for k4 = 1:nk
                        for k5 = 1:nb
                            Kalive(k1) = Kalive(k1) + kopt(k4,k2,k5,k1)*dist1ss(k4,k2,k5,k1,i3);
                            Kdead(k1) = Kdead(k1) + (1-surv(k1))*kopt(k4,k2,k5,k1)*dist1ss(k4,k2,k5,k1,i3);
                            Lab(k1) =Lab(k1)+ labopt(k4,k2,k5,k1)*dist1ss(k4,k2,k5,k1,i3);
                            ELab(k1) =ELab(k1)+ z(k2,k1,demtype)*labopt(k4,k2,k5,k1)*dist1ss(k4,k2,k5,k1,i3);
                        end
%                         node1 = node1+1;
                    end
                end
            end
%         end
            for k1 = Tr+1:T
                for k4 = 1:nk
                    for k5 = 1:nb
                        Kalive(k1) = Kalive(k1) + koptss(k4,k5,k1-Tr)*dist1ss_r(k4,k5,k1,i3);
                        Kdead(k1) = Kdead(k1) + (1-surv(k1))*koptss(k4,k5,k1-Tr)*dist1ss_r(k4,k5,k1,i3);
                    end
                end
            end
        end

        KPR = sum(Kalive+Kdead);
        ELAB = sum(ELab);

        
        jobdir = ['Steady State Values'];  % gives name
        filename = ['distvars_' num2str(demtype) '.mat'];
        totfile = fullfile(jobdir,filename);
        save(totfile,'dist1ss','dist1ss_r','Kalive','Kdead','KPR','ELab','Lab','ELAB','-double');


    end

    %*******************************************************************************************************************************************************************
    %*******************************************************************************************************************************************************************
    %*******************************************************************************************************************************************************************
    %*******************************************************************************************************************************************************************
%     for demtype = 1:ndem
%         [done] = dist_ss_moments(demtype,polno,impolno);
%     end
%     end
    toc
    kpr = 0;
    elab = 0;
    lab = 0;
    lab2 = 0;
    kdeadpr = 0;
    
%     plot([squeeze(sum(sum(sum(sum(permute(dist1ss,[1,2,3,5,4]))))))+squeeze(sum(sum(sum(permute(dist1ss_r,[1,2,4,3])))))])    % age distribution
    for demtype = 1:ndem
%         jobdir = ['job_' num2str(polno) '_' num2str(impolno)];  % gives name
        filename = ['distvars_' num2str(demtype) '.mat'];
        totfile = fullfile(jobdir,filename);
        load(totfile,'KPR','ELAB','Lab','Kdead');
        lab = lab + sum(Lab);
%         lab2 = lab2 + sum(Lab(1:Trmax));
        kpr = kpr+KPR;
        elab = elab+ELAB;
        kdeadpr = kdeadpr + sum(Kdead);
        
    end
    
    rhopr = (kpr-DD)/elab
    
    beqeps = beq - kdeadpr
    beq = kdeadpr/pop

    rhopr1 = (max(.5,kpr-DD))/elab;
    rhosseps = abs(rhopr1-rho)
    rho = .25*rhopr1 + .75*rho
    Y = A*((kpr - DD)^alp)*(elab^(1-alp))
%     rate = 1 + (A*alp).*(rho.^(alp-1)) - d
    wg = (A.*(1-alp)).*(rhopr.^(alp))
%     debteps = abs(DEBTss - .7*Y)
    DEBTss = .7*Y;
    
    K_Y = (kpr-DD)/Y
    I_Y = d*(kpr-DD)/Y
%     KK = kpr
%     Grev_to_GDP = sum(SSrev1+SSrev2+FEDREV)/Y
%     
%     agg_ss = sum(SSexp)/Y
%     avg_ss = sum((ben*(rpci/mpci)).*BENEFITS)/(sum(sum(sum(sum(distr(:,:,:,Tr:T))))))
%     avg_ss2 = sum(ben*(rpci/mpci).*BENEFITS)/sum(BENEFITS)
%     avg_ss3 = (rpci/mpci)*(sum(BENEFITS)/(sum(sum(sum(sum(distr(:,:,:,Tr:T)))))))
%     
%     full_time = sum(FT1+FT2)
%     part_time = sum(PT1+PT2)
%     not_working = sum(NW1+NW2)
        
end



display('done with ss')

% load params.mat

% flse = ((1-lab)/lab)*((1-chi3*(1-sigma4))/sigma4)
% flse2 = ((1-lab2)/lab2)*((1-chi3*(1-sigma4))/sigma4)
% moments

dist1 = zeros(nk,nz,nb,T,3,ndem);
distr = zeros(nk,nb,T,3,ndem);
% K_lifecycle = zeros(1,T);
% EL_lifecycle = zeros(1,T);

% this patches the distribution to add demographic type
for demtype = 1:ndem
    filename = ['distvars_' num2str(demtype) '.mat'];
    totfile = fullfile(jobdir,filename);
    load(totfile,'dist1ss','dist1ss_r','Kalive','Kdead');
    dist1(:,:,:,:,:,demtype) = dist1ss;
    distr(:,:,:,:,demtype) = dist1ss_r;
end

% jobdir = ['Steady State Values'];  % gives name
filename = ['eqmdist'];
totfile = fullfile(jobdir,filename);
save(totfile,'dist1','distr','-double');

save SSVALS Y wg rate beq kpr DEBTss rho pop pop_prev% sexp srev pit



dist1 = permute(dist1,[1,2,4,3]);  % [nk*nz,nb,ndem,T]
distr = permute(distr,[1,2,4,3]);
% load Surv_Probs
% deadpop = sum(squeeze(sum(sum(sum(dist1+distr)))).*(surv_proj(1,1:T))')
% avg_beq = kdead/deadpop
% b_mpci = avg_beq/mpci   % should be 1/3

% % distribution over benefits
% 
% dist1 = permute(dist1,[1,4,3,2]);
% distr = permute(distr,[1,4,3,2]);
% 
% plot(bgrid,squeeze(sum(sum(sum(dist1+distr)))))


% 
% for demtype = 1:ndem
%     [done] = dist_ss_moments(demtype)
% end
% 
% 
% NW1 = zeros(T,1);
% NW2 = zeros(T,1);
% PT1 = zeros(T,1);
% PT2 = zeros(T,1);
% FT1 = zeros(T,1);
% FT2 = zeros(T,1);
% SSexp = zeros(T,1);
% SSrev1 = zeros(T,1);
% SSrev2 = zeros(T,1);
% Klife = zeros(T,1);
% Llife = zeros(Tr,1);
% K0 = zeros(T,1);
% KUB = zeros(T,1);
% FEDREV = zeros(T,1);
% FIREV = zeros(T,1);
% BENEFITS=0;
% FRISCH = 0;
% WORKMASS = 0;
% for demtype = 1:ndem
%     filename = ['distmoments_' num2str(demtype) '.mat'];
%     totfile = fullfile(jobdir,filename);
%     load(totfile,'Benefits','SSrev','FIrev','Finc','SSexpss','SSrevss','FIrevss','Fincss','FT','PT','NW','FTss','PTss','NWss','FTz','PTz','NWz','FTzss','PTzss','NWzss','K_zero','FLSE','Working_mass','k_life','l_life','K_ub');
%     Klife = Klife + k_life;
%     Llife = Llife + l_life;
%     NW1 = NW1 + NW;
%     NW2 = NW2 + NWss;
%     FT1 = FT1 + FT;
%     FT2 = FT2 + FTss;
%     PT1 = PT1 + PT;
%     PT2 = PT2 + PTss;
%     SSexp = SSexp + SSexpss;
%     SSrev1 = SSrev1 + SSrev;
%     SSrev2 = SSrev2 + SSrevss;
%     BENEFITS = BENEFITS + Benefits;
%     K0 = K0 + K_zero;
%     KUB = KUB + K_ub;
%     FIREV = FIREV + FIrev;
%     FEDREV = FEDREV+FIrev+FIrevss;
%     FRISCH = FRISCH + sum(sum(FLSE));
%     WORKMASS = WORKMASS + sum(sum(Working_mass));
%     
% end
% Llife = Llife'./(Mu2(1:Tr)./sum(Mu2));
% Klife = Klife'./(Mu2./sum(Mu2));
% FITREV = FIREV'./(Mu2./sum(Mu2));
% save lifecycle_properties Llife Klife FEDREV 
% 
% sexp = sum(BENEFITS);  % SS expenditures
% srev = sum(SSrev1); % payroll tax revenue
% pit = sum(FIREV);  % personal income tax revenue
% 
% save SSVALS Y wg rate beq kpr DEBTss rho % sexp srev pit
% 
% 
% gexp = sum(FEDREV) + sum(SSrev1) - meanr_cbo*DEBTss;
% 
% save baselineG gexp DEBTss
% 
% frisch = FRISCH/WORKMASS
% 
% Y = A*((kpr-DD)^alp)*(elab^(1-alp))
% 
% K_Y = (kpr-DD)/Y
% 
% Grev_to_GDP = sum(SSrev1+SSrev2+FEDREV)/Y
% 
% agg_ss = sum(SSexp)/Y
% avg_ss = sum((ben(:,1)*(rpci/mpci)).*BENEFITS)/(sum(sum(sum(sum(distr(:,:,:,Tr:T))))))
% avg_ss2 = sum(ben(:,1)*(rpci/mpci).*BENEFITS)/sum(BENEFITS)
% avg_ss3 = (rpci/mpci)*(sum(BENEFITS)/(sum(sum(sum(sum(distr(:,:,:,Tr:T)))))))
% 
% full_time = sum(FT1+FT2)
% part_time = sum(PT1+PT2)
% not_working = sum(NW1+NW2)
% 
% % full_time_percent = (FT1+FT2)./squeeze(sum(sum(sum(dist1)))+sum(sum(sum(distr))));
% % part_time_percent = (PT1+PT2)./squeeze(sum(sum(sum(dist1)))+sum(sum(sum(distr))));
% % not_working_percent = (NW1+NW2)./squeeze(sum(sum(sum(dist1)))+sum(sum(sum(distr))));
% 
% toc
% 
% 
% 
% 

