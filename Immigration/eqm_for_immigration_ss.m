clear; clc;

polno = 1;
impolno = 1;

options = optimset('Display','off','TolFun',1.e-3);

jobdir = 'Steady State Values';
if (exist(jobdir, 'dir') == 7)
    rmdir(jobdir, 's')
end
mkdir(jobdir);

load('params.mat')
Tr = NRA_baseline;
pref_params = [sigma4,chi3];

load inctax_data.mat
avg_deduc = deduc_scale*avg_deduc;

load('Imm_Data.mat')

load('CBODATA.mat', 'r_cbo')
meanr_cbo = mean(r_cbo);


iter = 1;

rhosseps = 10;
rhosstol = .0001;
ssloc = 1;

im_scale = 1;    % scaling the immigration flow to target

load SSVALS.mat   % loads the last set of outputs as the first guess


% starting steady state equilbrium calculation
while (rhosseps > rhosstol)
    DD = DEBTss;  % debt used in the calculation of prices (zero for open economy)

    iteration_counter = [0,0];

    for demtype = 1:ndem
        
        filename = ['polparams' '_' num2str(polno) '.mat'];
        load(filename)
        
        % ssloc gives 1 if beginning of transition and 2 if end of transition
        if (ssloc == 1)
            ben = ss_benefit(:,1);
            beq = beq(1);
        elseif (ssloc == 2)
            ben = ss_benefit(:,Tss);
            beq = beq(end);
        end

        taxmax = taxmax(1);

        wage = (A*(1-alp))*(rho^(alp));
        rate1 = 1 + (A*alp)*(rho^(alp-1)) - d;
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

        Vss = zeros(nk,nb,T+1-Tr);
        laboptss = zeros(nk,nb,T-Tr);
        koptss = zeros(nk,nb,T-Tr);
        boptss = zeros(nk,nb,T-Tr);
        coptss = zeros(nk,nb,T-Tr);
        fedincomess = zeros(nk,nb,T-Tr);
        fitaxss = zeros(nk,nb,T-Tr);
        fsstaxss = zeros(nk,nb,T-Tr);

        % bequest motive
        Vbeq = phi1.*((1+kgrid./phi2).^(1-phi3));

        
        surv(T) = 0; %#ok<SAGROW>
        clear argmax val
        for t1 = T-Tr:-1:1
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

                    eff_inc = income - ftax - sstax + beq;   % effective income after tax and bequest.  just subtract savings to get consumption.
                    [argmax, val] = fminsearch(@(x) valfuncold(x,kgrid,eff_inc,EV(:,i1),pref_params),kgrid(k1),options);
                    
                    koptss(k1,i1,t1) = argmax;
                    Vss(k1,i1,t1) = -1*val;
                    coptss(k1,i1,t1) = eff_inc-argmax;
                    fitaxss(k1,i1,t1) = ftax;
                    
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
            cap_gain = (rate-1).*kgrid;
            cap_inc = rate.*kgrid;
            for i1 = 1:nb
                for j1 = 1:nz
                    eff_wage = wage*max(z(j1,t1,demtype),0);
                    EV = (1-surv(t1))*(Vbeq*ones(1,nb)); % gives continuation value dimensions [k',b']
                    for j2 = 1:nz
                        EV = EV + (surv(t1)*b)*tr_z(j1,j2)*squeeze(V(:,j2,:,t1+1));%
                    end
                    
                    EV = lt(t1,Tr)*EV + eq(t1,Tr)*E2;
                    
                    for k1 = 1:nk
                      
                        [argmax,val] = fminsearch(@(x) valfunc(x,kgrid,bgrid,cap_inc(k1),cap_gain(k1),eff_wage,beq,EV,pref_params,avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max,t1,i1,nb,nk),[kgrid(k1), .2],options);
                        income = cap_inc(k1) + eff_wage*argmax(2);
                        fincome = (rpci/mpci)*(cap_gain(k1) + eff_wage*argmax(2));
                        deduc = avg_deduc + coefs(1).*fincome + coefs(2).*fincome.^2 + coefs(3).*fincome.^3 + coefs(4).*fincome.^4;
                        ftax = (mpci/rpci)*limit*((max(fincome-deduc,0)) - ((max(fincome-deduc,0)).^(-X(1)) + (X(2))).^(-1/X(1)));
                        sstax = tau_ss*min(eff_wage*argmax(2),v_ss_max);
                        cons1 = income - ftax - sstax - argmax(1) + beq;

                        V(k1,j1,i1,t1) = -1*val;
                        kopt(k1,j1,i1,t1) = argmax(1);
                        labopt(k1,j1,i1,t1) = argmax(2);
                        copt(k1,j1,i1,t1) = cons1;
                        fedincome(k1,j1,i1,t1) = fincome;
                        fitax(k1,j1,i1,t1) = ftax;
                        fsstax(k1,j1,i1,t1) = sstax;
                        bopt(k1,j1,i1,t1) =(1/t1).*(bgrid(i1)*(t1-1) + min(v_ss_max,eff_wage*argmax(2)));
                        
                    end
                    

                end
            end
        end
        
        filename = ['sspol' num2str(demtype) '.mat'];
        totfile = fullfile(jobdir,filename);
        
        save(totfile,'V','Vss','copt','coptss','kopt','koptss','labopt','laboptss','fedincome','fedincomess','fitax','fitaxss','fsstax','fsstaxss','bopt','boptss','ben','bgrid');
        
    end


    
    for demtype = 1:ndem
        
        Kalive=zeros(1,T);
        Kdead=zeros(1,T);
        Lab = zeros(1,T);
        ELab = zeros(1,T);
        
        filename = ['sspol' num2str(demtype) '.mat'];
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
        
        while (disteps > disttol)
            
            dist1ss = zeros(nk,nz,nb,T,3);   % last dimension is [native, legal, illegal]
            dist1ss_r = zeros(nk,nb,T,3);
            iter = iter + 1;
            display(iter)
            iteration_counter(demtype) = iter;
       
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
            
            dist_age = squeeze(sum(sum(sum(sum(permute(dist1ss,[1,2,3,5,4]))))))+squeeze(sum(sum(sum(permute(dist1ss_r,[1,2,4,3])))));
            death_rate = sum(dist_age'.*(1-surv(1:T)))/pop;  % Just for debug
            % Below, convergence is relative to
            disteps = max(abs((1/dist_age(1)).*dist_age(2:end) - (1/dist_age_previous(1)).*dist_age_previous(2:end)));
            display(disteps)
            dist_age_previous = dist_age;
            dist1ss_previous = dist1ss;   % updating working-age distribution array
            dist1ss_r_previous = dist1ss_r;   % updating retiree distribution array
            
        end
            
        immigrant_population = 1-(sum(sum(sum(sum(dist1ss(:,:,:,:,1)))))+sum(sum(sum(dist1ss_r(:,:,:,1)))))/(sum(sum(sum(sum(sum(dist1ss)))))+sum(sum(sum(sum(dist1ss_r)))));
        
        % solving for aggregates by age
        for i3 = 1:3
            for k1 = 1:Tr
                for k2 = 1:nz
                    for k4 = 1:nk
                        for k5 = 1:nb
                            Kalive(k1) = Kalive(k1) + kopt(k4,k2,k5,k1)*dist1ss(k4,k2,k5,k1,i3);
                            Kdead(k1) = Kdead(k1) + (1-surv(k1))*kopt(k4,k2,k5,k1)*dist1ss(k4,k2,k5,k1,i3);
                            Lab(k1) =Lab(k1)+ labopt(k4,k2,k5,k1)*dist1ss(k4,k2,k5,k1,i3);
                            ELab(k1) =ELab(k1)+ z(k2,k1,demtype)*labopt(k4,k2,k5,k1)*dist1ss(k4,k2,k5,k1,i3);
                        end
                    end
                end
            end
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
        
        filename = ['distvars_' num2str(demtype) '.mat'];
        totfile = fullfile(jobdir,filename);
        save(totfile,'dist1ss','dist1ss_r','Kalive','Kdead','KPR','ELab','Lab','ELAB','-double');
        
    end
    
    kpr = 0;
    elab = 0;
    lab = 0;
    lab2 = 0;
    kdeadpr = 0;
    
    for demtype = 1:ndem
        filename = ['distvars_' num2str(demtype) '.mat'];
        totfile = fullfile(jobdir,filename);
        load(totfile,'KPR','ELAB','Lab','Kdead');
        lab = lab + sum(Lab);
        kpr = kpr+KPR;
        elab = elab+ELAB;
        kdeadpr = kdeadpr + sum(Kdead);
    end
    
    rhopr = (kpr-DD)/elab;
    
    beqeps = beq - kdeadpr;
    beq = kdeadpr/pop;

    rhopr1 = (max(.5,kpr-DD))/elab;
    rhosseps = abs(rhopr1-rho);
    rho = .25*rhopr1 + .75*rho;
    Y = A*((kpr - DD)^alp)*(elab^(1-alp));
    wg = (A.*(1-alp)).*(rhopr.^(alp));
    DEBTss = .7*Y;
    
    K_Y = (kpr-DD)/Y;
    I_Y = d*(kpr-DD)/Y;
    
end


dist1 = zeros(nk,nz,nb,T,3,ndem);
distr = zeros(nk,nb,T,3,ndem);

% this patches the distribution to add demographic type
for demtype = 1:ndem
    filename = ['distvars_' num2str(demtype) '.mat'];
    totfile = fullfile(jobdir,filename);
    load(totfile,'dist1ss','dist1ss_r','Kalive','Kdead');
    dist1(:,:,:,:,:,demtype) = dist1ss;
    distr(:,:,:,:,demtype) = dist1ss_r;
end


filename = 'eqmdist';
totfile = fullfile(jobdir,filename);
save(totfile,'dist1','distr','-double');




%% Testing

distvars_1   = load(fullfile(jobdir  , 'distvars_1.mat'));
distvars_1_0 = load(fullfile('Freeze', 'distvars_1.mat'));

fprintf('distvars_1\n');
fieldnames = fields(distvars_1)';
for fieldname = fieldnames
    fprintf('\t%-12s%d\n', fieldname{1}, all(distvars_1.(fieldname{1})(:) == distvars_1_0.(fieldname{1})(:)));
end
fprintf('\n');


distvars_2   = load(fullfile(jobdir  , 'distvars_2.mat'));
distvars_2_0 = load(fullfile('Freeze', 'distvars_2.mat'));

fprintf('distvars_2\n');
fieldnames = fields(distvars_2)';
for fieldname = fieldnames
    fprintf('\t%-12s%d\n', fieldname{1}, all(distvars_2.(fieldname{1})(:) == distvars_2_0.(fieldname{1})(:)));
end
fprintf('\n');

