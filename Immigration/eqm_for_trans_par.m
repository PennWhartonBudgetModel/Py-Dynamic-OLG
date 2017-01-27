% function [done] = eqm_for_trans(polno,impolno)

%-------------------------------------
% Place main_HPC2 material here
%-------------------------------------

% clear;clc;
% polno=1;
% impolno=1;
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

filename = ['imm_polparams' '_' num2str(polno) '.mat'];
load(filename,'pop_trans')

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
    
    %*********************************************************************************************************
    %*********************************************************************************************************
    %*********************************************************************************************************
    %**********************DIST_TRANS***********************************************************************
    
    % function reads: function [KPR,ELAB,OUTPUT,Beq,FEDIT,SSEXP,SSREV,LFP] = dist_trans(rho1,beq,polno,impolno)
    for demtype = 1:ndem
%     demtype
%         for tt1 = 1:Tss
        for tt1 = 1:Tss
            cohort = tt1+(T-1);
            
            %-----------------------------------------------------------------------------------------------------------------------------
            %-----------------------------------------------------------------------------------------------------------------------------
            %-----------------------------------------------------------------------------------------------------------------------------
            %-----------------------------HEAD-----------------------------------------------------------------------------------------

            % function reads: function [hello] = head4(rho3,beq1,agetrans,tt1,Tss,demtype,polno,impolno)
            


            Tr = NRA(cohort);
            ben = ss_benefit;   % renamed for parallelization non-conflict
            
            %-------------------------------------------------------------------------------------------------------------------------------------------------------------
            % agetrans will give the age of individual when the transition ended (i.e., the age in the LAST period of the transition, not the 1st pd of ss).
            % variables loaded are ss policy and value functions in their entirety.
            % when needed, the values are replaced by this algorithm
            %-------------------------------------------------------------------------------------------------------------------------------------------------------------
            surv1 = surv;
            rho2 = rho1;
            beq2 = beq;
%             r_portfolio1 = r_portfolio;
            KK2 = KK;
            DD2 = DD;
            r_cbo1 = r_cbo;
            taxmax1 = taxmax;
            ss_tax1 = ss_tax;
            rho = [rho2(tt1:Tss), rho2(end).*ones(1,T-(Tss-(tt1-1)))];
            beq1 = [beq2(tt1:end), beq2(end).*ones(1,T-(Tss-(tt1-1)))];
            DD1 = [DD2(tt1:Tss), DD2(end).*ones(1,T-(Tss-(tt1-1)))];
            KK1 = [KK2(tt1:Tss), KK2(end).*ones(1,T-(Tss-(tt1-1)))];
            r_cbo1 = [r_cbo1(tt1:Tss), r_cbo1(end).*ones(1,T-(Tss-(tt1-1)))];


            wage = (A*(1-alp)).*(rho.^(alp));
            rate1 = 1 + (A*alp).*(rho.^(alp-1)) - d;
            r_portfolio1 = (rate1.*KK1 + (1+r_cbo1).*DD1)./(KK1+DD1);


            % %     % Creating arrays
                V = zeros(nk,nz,nb,Tr+1);
                labopt = zeros(nk,nz,nb,Tr);
                kopt = zeros(nk,nz,nb,Tr);
                copt = zeros(nk,nz,nb,Tr);
                bopt = zeros(nk,nz,nb,Tr);
                fedincome = zeros(nk,nz,nb,Tr);
                fitax = zeros(nk,nz,nb,Tr);
                fsstax = zeros(nk,nz,nb,Tr);
                ss_select = zeros(nk,nz,nb,Tr);
                benopt = zeros(nk,nz,nb,Tr);
                ss_base = zeros(nk,nz,nb,Tr);

                Vss = zeros(nk,nb,T+1-Tr);
                laboptss = zeros(nk,nb,T-Tr);
                koptss = zeros(nk,nb,T-Tr);
                boptss = zeros(nk,nb,T-Tr);
                coptss = zeros(nk,nb,T-Tr);
                fedincomess = zeros(nk,nb,T-Tr);
                fitaxss = zeros(nk,nb,T-Tr);
                fsstaxss = zeros(nk,nb,T-Tr);
                benoptss = zeros(nk,nb,T-Tr);
                ss_basess = zeros(nk,nb,T-Tr);
            
            % bequest motive
            Vbeq = phi1.*((1+kgrid./phi2).^(1-phi3));
            % disutility of working

            
%             clear argmax val
            for t1 = T-Tr:-1:1
                age = t1+Tr;
                year = max(1,min(cohort +age - T,Tss));   % year = tt1 + age - 1 = cohort + age - T = tt1+ T - 1 + age - T = tt1 + age - 1; check.
                EV = (1-surv1(t1+Tr))*Vbeq*ones(1,nb); % gives continuation value dimensions [k',b']
                EV = EV + (surv1(t1+Tr)*b)*squeeze(Vss(:,:,t1+1));
                for i1 = 1:nb
                    for k1 = 1:nk
                        totben = ben(i1,year);% (ben_change(t1+Tr)*eq(yr_start(t1+Tr),1) + (1- eq(yr_start(t1+Tr),1))*1).*(ben(i1).*(1 + b_increase*b_inc1(t1+Tr)*ge(t1+Tr+19,inc_age)));
                        income = r_portfolio1(age).*kgrid(k1) +  totben;
                        fincome = (rpci/mpci)*((max(0,r_portfolio1(age)-1)).*kgrid(k1) + (1-ss_tax_cred)*totben);
                        deduc = avg_deduc + coefs(1).*fincome + coefs(2).*fincome.^2 + coefs(3).*fincome.^3 + coefs(4).*fincome.^4;
                        ftax = (mpci/rpci)*limit.*((max(fincome-deduc,0)) - ((max(fincome-deduc,0)).^(-X(1)) + (X(2))).^(-1/X(1)));
                        eff_inc = income - ftax + beq1(t1+Tr);
                        [argmax, val] = fminsearch(@(x) valfuncold(x,kgrid,eff_inc,EV(:,i1),pref_params),kgrid(k1),options);
                        
                        Vss(k1,i1,t1) = -1*val;
                        koptss(k1,i1,t1) = argmax;
                        coptss(k1,i1,t1) = eff_inc - argmax;
                        fedincomess(k1,i1,t1) = fincome*(mpci/rpci);
                        fitaxss(k1,i1,t1) = ftax;
                        benoptss(k1,i1,t1) = totben;
                    end
                end
            end
            
%             [V,kopt,labopt,copt,fedincome,fitax,fsstax,bopt,ss_base] = workAgent(Vss,Vbeq,nb,nk,nz,Tr,tr_z,b,surv1,wage,z,rate,tau_ss,kgrid,bgrid,beq1,pref_params,avg_deduc,coefs,limit,X,mpci,rpci,v_ss_max)
            
%             coder.ceval('



            % Pre-retirement age
            E2 = 0;
            for j1 = 1:nz
                for j2 = 1:nz
                    E2 = E2 + surv1(Tr)*b*tr_z(j1,j2)*squeeze(Vss(:,:,1));
                end
            end
            
            
            for t1 = Tr:-1:1
                age = t1;
                year = max(1,min(cohort +age - T,Tss));
                cap_gain = (r_portfolio1(t1)-1).*kgrid;
                cap_inc = r_portfolio1(t1).*kgrid;
                v_ss_max1 = taxmax1(year);
                tau_ss1 = ss_tax1(year);
            %     t1
                for i1 = 1:nb
                    for j1 = 1:nz
                        eff_wage = wage(t1)*max(z(j1,t1,demtype),0);    % effective income after tax and bequest.  just subtract savings to get consumption.
                        EV = (1-surv1(t1))*(Vbeq*ones(1,nb)); % gives continuation value dimensions [k',b']
                        for j2 = 1:nz
                            EV = EV + (surv1(t1)*b)*tr_z(j1,j2)*squeeze(V(:,j2,:,t1+1));%
                        end
                        
                        EV = lt(t1,Tr)*EV + eq(t1,Tr)*E2;    %
                        for k1 = 1:nk
                            
                            [argmax,val] = fminsearch(@(x) valfunc(x,kgrid,bgrid,cap_inc(k1),cap_gain(k1),eff_wage,beq1(t1),EV,pref_params,avg_deduc, coefs, limit, X, mpci, rpci, tau_ss1, v_ss_max,t1,i1,nb,nk),[kgrid(k1), .2],options);
                            income = cap_inc(k1) + eff_wage*argmax(2);
                            fincome = (rpci/mpci)*(cap_gain(k1) + eff_wage*argmax(2));%  + soc_sec*ge(t1,Tr);
                            deduc = avg_deduc + coefs(1).*fincome + coefs(2).*fincome.^2 + coefs(3).*fincome.^3 + coefs(4).*fincome.^4;
                            ftax = (mpci/rpci)*limit*((max(fincome-deduc,0)) - ((max(fincome-deduc,0)).^(-X(1)) + (X(2))).^(-1/X(1)));
                            sstax = tau_ss1*min(eff_wage*argmax(2),v_ss_max1);
                            ssbase = min(eff_wage*argmax(2),v_ss_max1);
                            cons1 = income - ftax - sstax - argmax(1) + beq1(t1);
 
                            V(k1,j1,i1,t1) = -1*val;
                            kopt(k1,j1,i1,t1) = argmax(1);
                            labopt(k1,j1,i1,t1) = argmax(2);

                            copt(k1,j1,i1,t1) = cons1;
                            fedincome(k1,j1,i1,t1) = fincome;
                            fitax(k1,j1,i1,t1) = ftax;
                            fsstax(k1,j1,i1,t1) = sstax;
                            bopt(k1,j1,i1,t1) = (1/t1).*(bgrid(i1)*(t1-1) + min(eff_wage*argmax(2),v_ss_max)); % could also do bprime(Il);
                            ss_base(k1,j1,i1,t1) = ssbase;
                        end
                    end
                end
            end

            filename = ['head' num2str(tt1) '_' num2str(demtype) '_' num2str(polno) '.mat'];  % abs(agetrans - Tss - 1) = tt1.  just minimizing function inputs here.
            jobdir = ['job_' num2str(polno) '_' num2str(impolno)];  % gives name
            totfile = fullfile(jobdir,filename);
            save_par(totfile,V,Vss,copt,coptss,kopt,koptss,labopt,laboptss,fedincome,fedincomess,fitax,fitaxss,fsstax,fsstaxss,ss_select,bopt,boptss,Tr,ben,benopt,benoptss,ss_base,ss_basess);
%             save(totfile,'V','Vss','copt','coptss','kopt','koptss','labopt','laboptss','fedincome','fedincomess','fitax','fitaxss','fsstax','fsstaxss','ss_select','bopt','boptss','Tr','ben','benopt','benoptss','ss_base','ss_basess');





            %-----------------------------HEAD-----------------------------------------------------------------------------------------
            %-----------------------------------------------------------------------------------------------------------------------------
            %-----------------------------------------------------------------------------------------------------------------------------
            %-----------------------------------------------------------------------------------------------------------------------------
        end
        
%         Ttail=2;



        for Ttail = 2:T
            cohort = T-(Ttail-1);
            
            %-----------------------------------------------------------------------------------------------------------------------------
            %-----------------------------------------------------------------------------------------------------------------------------
            %-----------------------------------------------------------------------------------------------------------------------------
            %-----------------------------TAIL-------------------------------------------------------------------------------------------
            

            Tr = NRA(cohort);
            ben = ss_benefit;
            surv1 = surv;

            %-----------------------------------------------------------------------------------------------------------
            %*******************NOTE: Ttail gives the age of person in 1st pd of transition*****************************
            %-----------------------------------------------------------------------------------------------------------
            
            beq2 = beq;
            rho2 = rho1;
            KK2 = KK;
            DD2 = DD;
            r_cbo1 = r_cbo;
            taxmax1 = taxmax;
            ss_tax1 = ss_tax;

            rho = [zeros(1,Ttail-1),rho2(1:Tss)];
            rho = [rho, rho(end).*ones(1,T)];
            KK2 = [zeros(1,Ttail-1),KK2(1:Tss)];
            KK1 = [KK2, KK2(end).*ones(1,T)];
            DD2 = [zeros(1,Ttail-1),DD2(1:Tss)];
            DD1 = [DD2, DD2(end).*ones(1,T)];
            r_cbo1 = [zeros(1,Ttail-1),r_cbo1(1:Tss)];
            r_cbo1 = [r_cbo1, r_cbo1(end).*ones(1,T)];
            beq1 = [zeros(1,Ttail-1),beq2(1:Tss)];
            beq1 = [beq1, beq1(end).*ones(1,T)];
            wage = (A*(1-alp)).*(rho.^(alp));
            rate1 = 1 + (A*alp).*(rho.^(alp-1)) - d;
            r_portfolio1 = (rate1.*KK1 + (1+r_cbo1).*DD1)./(KK1+DD1);

            % Creating arrays

%             if Ttail<Tr+1
                V = zeros(nk,nz,nb,T- (Ttail-1)+1);
                labopt = zeros(nk,nz,nb,T- (Ttail-1));
                kopt = zeros(nk,nz,nb,T- (Ttail-1));
                copt = zeros(nk,nz,nb,T- (Ttail-1));
                bopt = zeros(nk,nz,nb,T- (Ttail-1));
                fedincome = zeros(nk,nz,nb,T- (Ttail-1));
                fitax = zeros(nk,nz,nb,T- (Ttail-1));
                fsstax = zeros(nk,nz,nb,T- (Ttail-1));
                ss_select = zeros(nk,nz,nb,T- (Ttail-1));
                benopt = zeros(nk,nz,nb,T- (Ttail-1));
                ss_base = zeros(nk,nz,nb,T- (Ttail-1));
%             end



            % if Ttail<Tr
            Vss = zeros(nk,nb,min(T+1-Tr,T+1-(Ttail-1)));
            laboptss = zeros(nk,nb,min(T-Tr,T-(Ttail-1)));
            koptss = zeros(nk,nb,min(T-Tr,T-(Ttail-1)));
            boptss = zeros(nk,nb,min(T-Tr,T-(Ttail-1)));
            coptss = zeros(nk,nb,min(T-Tr,T-(Ttail-1)));
            fedincomess = zeros(nk,nb,min(T-Tr,T-(Ttail-1)));
            fitaxss = zeros(nk,nb,min(T-Tr,T-(Ttail-1)));
            fsstaxss = zeros(nk,nb,min(T-Tr,T-(Ttail-1)));
            benoptss = zeros(nk,nb,min(T-Tr,T-(Ttail-1)));
            ss_basess = zeros(nk,nb,min(T-Tr,T-(Ttail-1)));
            % end


            % ytax = zeros(nk,nz,T);
            % ytaxst = zeros(nk,nz,T);
            % ytaxf = zeros(nk,nz,T);
            % transf = zeros(nk,nz,T);

            % bequest motive
            Vbeq = phi1.*((1+kgrid./phi2).^(1-phi3));
            
            % disutility of working



%             surv(T) =0;
            
            
            %--------------------------------------------------------------------------------    REGION 3    --------------------------------------------------------------------------------------


            if Ttail<=Tr % region 3
%                 clear argmax val

                % Retirement age, accepted ss
                for t1 = T-Tr:-1:1
                    age = t1+Tr;    % cohort = T-(Ttail-1);
                    year = max(1,min(age+cohort-T,Tss));  %  age - (Ttail-1) = age + cohort - T = age - T + (Tail-1)
                    totben = ben(:,year);
                    EV = (1-surv1(age))*Vbeq*ones(1,nb); % gives continuation value dimensions [k',b']
                    EV = EV + (surv1(age)*b)*squeeze(Vss(:,:,t1+1));
%                     E = repmat(EV,[1,nk]); 
                    for i1 = 1:nb
                        for k1 = 1:nk
                            income = r_portfolio1(age).*kgrid(k1)  +  totben(i1);
                            fincome = (rpci/mpci)*((max(0,r_portfolio1(age)-1)).*kgrid(k1)  + (1-ss_tax_cred)*totben(i1));
                            deduc = avg_deduc + coefs(1).*fincome + coefs(2).*fincome.^2 + coefs(3).*fincome.^3 + coefs(4).*fincome.^4;
                            ftax = (mpci/rpci)*limit.*((max(fincome-deduc,0)) - ((max(fincome-deduc,0)).^(-X(1)) + (X(2))).^(-1/X(1)));
                            sstax = 0;
                            eff_inc = income - ftax - sstax + beq1(age);    % effective income after tax and bequest.  just subtract savings to get consumption.
                            [argmax, val] = fminsearch(@(x) valfuncold(x,kgrid,eff_inc,EV(:,i1),pref_params),kgrid(k1),options);
                            
                            Vss(k1,i1,t1) = -1*val;
                            koptss(k1,i1,t1) = argmax;  
                            
                            coptss(k1,i1,t1) = eff_inc - argmax;
                            fedincomess(k1,i1,t1) = fincome;
                            fitaxss(k1,i1,t1) = ftax;
                            benoptss(k1,i1,t1) = totben(i1);
                        end
                    end
                end
            
                % Pre-retirement age
                E2 = 0;
                for j1 = 1:nz
                    for j2 = 1:nz
                        E2 = E2 + surv1(Tr)*b*tr_z(j1,j2)*squeeze(Vss(:,:,1));
                    end
                end


                for t1 = Tr-(Ttail-1):-1:1
                    age = t1+Ttail-1;
                    year = max(1,min(cohort +age - T,Tss));
                    cap_gain = (r_portfolio1(age)-1).*kgrid;
                    cap_inc = r_portfolio1(age).*kgrid;
                    v_ss_max1 = taxmax1(year);
                    tau_ss1 = ss_tax1(year);
                    if (v_ss_max1==0)||(wage(age)==0)
                        display('problem')
                        tt1
                        t1
                    end
            %         t1
                    for i1 = 1:nb
                        for j1 = 1:nz
                            eff_wage = wage(age)*max(z(j1,age,demtype),0);%  + soc_sec*ge(t1,Tr);
                            EV = (1-surv1(age))*(Vbeq*ones(1,nb)); % gives continuation value dimensions [k',b']
                            for j2 = 1:nz
                                EV = EV + (surv1(age)*b)*tr_z(j1,j2)*squeeze(V(:,j2,:,t1+1));%
                            end
                            EV = lt(age,Tr)*EV + eq(age,Tr)*E2;
                            for k1 = 1:nk
                                
                                [argmax,val] = fminsearch(@(x) valfunc(x,kgrid,bgrid,cap_inc(k1),cap_gain(k1),eff_wage,beq1(age),EV,pref_params,avg_deduc, coefs, limit, X, mpci, rpci, tau_ss1, v_ss_max,age,i1,nb,nk),[kgrid(k1), .2],options);
                                
                                income = cap_inc(k1) + eff_wage*argmax(2);
                                fincome = (rpci/mpci)*(cap_gain(k1) + eff_wage*argmax(2));%  + soc_sec*ge(t1,Tr);
                                deduc = avg_deduc + coefs(1).*fincome + coefs(2).*fincome.^2 + coefs(3).*fincome.^3 + coefs(4).*fincome.^4;
                                ftax = (mpci/rpci)*limit*((max(fincome-deduc,0)) - ((max(fincome-deduc,0)).^(-X(1)) + (X(2))).^(-1/X(1)));
                                sstax = tau_ss1*min(eff_wage*argmax(2),v_ss_max1);
                                ssbase = min(eff_wage*argmax(2),v_ss_max1);
                                cons1 = income - ftax - sstax - argmax(1) + beq1(age);


                                V(k1,j1,i1,t1) = -1*val;

                                kopt(k1,j1,i1,t1) = argmax(1);  %(nl,nk,nk,4,T-Tr,nep,nh,ng)
                                labopt(k1,j1,i1,t1) = argmax(2);
% 
                                copt(k1,j1,i1,t1) = cons1;
                                fedincome(k1,j1,i1,t1) = fincome;
                                fitax(k1,j1,i1,t1) = ftax;
                                fsstax(k1,j1,i1,t1) = sstax;
                                bopt(k1,j1,i1,t1) = (1/(age)).*(bgrid(i1)*(age-1) + min(v_ss_max,eff_wage*argmax(2)));%(1/t1).*(bgrid(i1)*(t1-1) + wage*z(j1,t1,demtype).*lgrid(Il)); % could also do bprime(Il);
                                ss_base(k1,j1,i1,t1) = ssbase;
                            end
                        end
                    end
                end
%                 phi1

            end

            %---------------------------------------------------------------------------------     REGION 1      ------------------------------------------------------------------------------------------

            if Ttail>Tr  % region 1
%                 clear argmax val
                % Retirement age, accepted ss
            %     T-(Ttail-1)
                for t1 = T-(Ttail-1):-1:1
                    age = t1+Ttail-1;
                    year = max(1,min(age+cohort-T,Tss));
                    totben = ben(:,year);
                    EV = (1-surv1(age))*Vbeq*ones(1,nb); % gives continuation value dimensions [k',b']
                    EV = EV + (surv1(age)*b)*squeeze(Vss(:,:,t1+1));
                    for i1 = 1:nb
                        for k1 = 1:nk
                            income = r_portfolio1(age).*kgrid(k1) +  totben(i1);
                            fincome = (rpci/mpci)*((max(0,r_portfolio1(age)-1)).*kgrid(k1) + (1-ss_tax_cred)*totben(i1));
                            deduc = avg_deduc + coefs(1).*fincome + coefs(2).*fincome.^2 + coefs(3).*fincome.^3 + coefs(4).*fincome.^4;
                            ftax = (mpci/rpci)*limit.*((max(fincome-deduc,0)) - ((max(fincome-deduc,0)).^(-X(1)) + (X(2))).^(-1/X(1)));
                            eff_inc = income - ftax + beq1(age);
                            [argmax, val] = fminsearch(@(x) valfuncold(x,kgrid,eff_inc,EV(:,i1),pref_params),kgrid(k1),options);

                            cons1 = eff_inc - argmax;

                            Vss(k1,i1,t1) = -1*val;
                            koptss(k1,i1,t1) = argmax;  %(nl,nk,nk,4,T-Tr,nep,nh,ng)
                            benoptss(k1,i1,t1) = totben(i1);
                            coptss(k1,i1,t1) = cons1;
                            fedincomess(k1,i1,t1) = fincome;
                            fitaxss(k1,i1,t1) = ftax;
                        end
                    end
                end
%                 phi1

            end


            filename = ['tail' num2str(Ttail) '_' num2str(demtype) '_' num2str(polno)  '.mat'];
            jobdir = ['job_' num2str(polno) '_' num2str(impolno)];  % gives name
            totfile = fullfile(jobdir,filename);
            
            
            if Ttail<Tr+1
                save_par(totfile,V,Vss,copt,coptss,kopt,koptss,labopt,laboptss,fedincome,fedincomess,fitax,fitaxss,fsstax,fsstaxss,ss_select,bopt,boptss,Tr,ben,benopt,benoptss,ss_base,ss_basess);
            else 
                save_par2(totfile,Vss,coptss,koptss,laboptss,fedincomess,fitaxss,fsstaxss,boptss,Tr,ben,benoptss,ss_basess,wage,rate);
%                 save(totfile,'Vss','coptss','koptss','laboptss','fedincomess','fitaxss','fsstaxss','boptss','Tr','ben','benoptss','ss_basess','wage','rate');
            end

            
%             display('pause')
            
            %-----------------------------TAIL-------------------------------------------------------------------------------------------
            %-----------------------------------------------------------------------------------------------------------------------------
            %-----------------------------------------------------------------------------------------------------------------------------
            %-----------------------------------------------------------------------------------------------------------------------------
        end
    end
    
    
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
        trans_thread_head(tt1,polno,impolno);
    end

    for tt1 =1:Tss
        for demtype = 1:ndem
            filename = ['transvars_' num2str(demtype) '_' num2str(tt1) '_head_' num2str(polno) '.mat'];
            jobdir = ['job_' num2str(polno) '_' num2str(impolno)];  % gives name
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
        trans_thread_tail(tt1,polno,impolno);
    end

    for tt1 = 2:T  

        for demtype = 1:ndem
            filename = ['transvars_' num2str(demtype) '_' num2str(tt1) '_tail_' num2str(polno) '.mat'];
            jobdir = ['job_' num2str(polno) '_' num2str(impolno)];  % gives name
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
    save(filename,'dist_1','dist_r','KPR','LAB','ELAB','OUTPUT','DIST','FEDIT','SSEXP','SSREV','LFP','SSBASE');
    
    
    
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
    filename = ['current_results_' num2str(polno) '_' num2str(impolno) '.mat'];
    save(filename)
    

end

% for t1 = 1:Tss-1
%     DEBT(t1+1) = GEXP(t1) - GREV(t1) + DEBT(t1)*r_g(t1) + ssexp(t1);
% end

delete(filename)

filename = ['results_' num2str(polno) '_' num2str(impolno)  '.mat'];
save(filename)

done = 1;

% delete(mypool)
    

