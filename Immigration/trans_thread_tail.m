function [] = trans_thread_tail(startyear, polno, impolno)

load('params.mat')
load('Surv_Probs.mat')
load('Imm_Data.mat')

jobdir = 'Testing';
load(fullfile(jobdir, sprintf('imm_polparams_%u.mat', impolno)))

dist_1 = zeros(nk,nz,nb,T+startyear+1,3,ndem);
dist_r = zeros(nk,   nb,T+startyear+1,3,ndem);

load(fullfile(jobdir, 'eqmdist.mat'));
load('SSVALS.mat', 'pop_prev')
pop_trans = [pop_prev; pop_trans]; %#ok<NODEF>


for idem = 1:ndem
    
    load(fullfile('Freeze', 'Cohorts', sprintf('tail%u_%u_%u.mat', -startyear+1, idem, polno)));
    
    T_life   = T;
    T_work   = Tr;
    T_model  = Tss;
    T_past   = max(-startyear, 0);
    T_shift  = max(+startyear, 0);
    T_active = min(startyear+T_life, T_model) - T_shift;
    
    
    dist_1(:,:,:,1,:,idem) = dist1(:,:,:,-startyear+1,:,idem); %#ok<NODEF>
    dist_r(:,  :,1,:,idem) = distr(:,  :,-startyear+1,:,idem); %#ok<NODEF>
    
    
    for t = 1:T_work+startyear
        
        age  = t + T_past;
        year = max(1, min(startyear+age, T_model)) + 1; % (Addition of 1 may be erroneous; leads to year starting at 2)
        
        im_flow = [ 0                                                ;
                    pop_trans(year) * imm_age(age) * legal_rate(1)   ;
                    pop_trans(year) * imm_age(age) * illegal_rate(1) ];
        
        for iz = 1:nz
            for ik = 1:nk
                for ib = 1:nb
                    
                    point_k = max(kopt(ik,iz,ib,t), kgrid(1));
                    loc1 = find(kgrid(1:nk-1) <= point_k, 1, 'last');
                    w1 = (point_k - kgrid(loc1)) / (kgrid(loc1+1) - kgrid(loc1));
                    w1 = min(w1, 1);
                    
                    point_b = max(bopt(ik,iz,ib,t), bgrid(1));
                    loc2 = find(bgrid(1:nb-1) <= point_b, 1, 'last');
                    w2 = (point_b - bgrid(loc2)) / (bgrid(loc2+1) - bgrid(loc2));
                    w2 = min(w2, 1);
                    
                    for jz = 1:nz
                        for ipop = 1:3
                            
                            dist_hold = dist_1(ik,iz,ib,t,ipop,idem) + (ik == 1)*(ib == 1)*proddist_age(iz,age,ipop)*im_flow(ipop);
                            
                            if (age < T_work)
                                dist_1(loc1  , jz, loc2  , t+1, ipop, idem) = dist_1(loc1  , jz, loc2  , t+1, ipop, idem) + surv(age)*(1-w2)*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_1(loc1+1, jz, loc2  , t+1, ipop, idem) = dist_1(loc1+1, jz, loc2  , t+1, ipop, idem) + surv(age)*(1-w2)*(w1  )*tr_z(iz,jz)*dist_hold;
                                dist_1(loc1  , jz, loc2+1, t+1, ipop, idem) = dist_1(loc1  , jz, loc2+1, t+1, ipop, idem) + surv(age)*(w2  )*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_1(loc1+1, jz, loc2+1, t+1, ipop, idem) = dist_1(loc1+1, jz, loc2+1, t+1, ipop, idem) + surv(age)*(w2  )*(w1  )*tr_z(iz,jz)*dist_hold;
                            else
                                dist_r(loc1  ,     loc2  , t+1, ipop, idem) = dist_r(loc1  ,     loc2  , t+1, ipop, idem) + surv(age)*(1-w2)*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_r(loc1+1,     loc2  , t+1, ipop, idem) = dist_r(loc1+1,     loc2  , t+1, ipop, idem) + surv(age)*(1-w2)*(w1  )*tr_z(iz,jz)*dist_hold;
                                dist_r(loc1  ,     loc2+1, t+1, ipop, idem) = dist_r(loc1  ,     loc2+1, t+1, ipop, idem) + surv(age)*(w2  )*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_r(loc1+1,     loc2+1, t+1, ipop, idem) = dist_r(loc1+1,     loc2+1, t+1, ipop, idem) + surv(age)*(w2  )*(w1  )*tr_z(iz,jz)*dist_hold;
                            end
                            
                        end
                    end
                    
                end
            end
        end
        
        if amnesty
            
            if (age < T_work)
                
                amnesty_dist = squeeze(amnesty*dist_1(:,:,:,t+1,3,idem)); 
                amnesty_dist = permute(amnesty_dist,[2,1,3]);
                amnesty_dist = squeeze(sum(amnesty_dist)); % amnesty_dist has dimensions nk,nb
                
                prod_legal = squeeze(dist_1(:,:,:,t+1,2,idem));   % has dimensions nk,nz,nb
                prod_legal = squeeze(sum(sum(permute(prod_legal,[1,3,2]))));
                prod_legal = prod_legal / sum(prod_legal);
                
                for ik = 1:nk
                    for ib = 1:nb
                        for iz = 1:nz
                            dist_1(ik,iz,ib,t+1,2,idem) = dist_1(ik,iz,ib,t+1,2,idem) + prod_legal(iz)*amnesty_dist(ik,ib);
                        end
                    end
                end
                dist_1(:,:,:,t+1,3,idem) = (1-amnesty)*dist_1(:,:,:,t+1,3,idem);
                
            else
                
                dist_r(:,:,t+1,2,idem) = dist_r(:,:,t+1,2,idem) + amnesty*dist_r(:,:,t+1,3,idem);
                dist_r(:,:,t+1,3,idem) = (1-amnesty)*dist_r(:,:,t+1,3,idem);
                
            end
            
        end
        
        dist_1(:,:,:,t+1,3,idem) = (1-deportation)*dist_1(:,:,:,t+1,3,idem);  
        dist_r(:,  :,t+1,3,idem) = (1-deportation)*dist_r(:,  :,t+1,3,idem);
        
    end
    
    
    for t = max(T_work+startyear+1, 1):T_life+startyear-1
        
        age  = t + T_past;
        year = max(1, min(startyear+age, T_model)) + 1; % (Addition of 1 may be erroneous; leads to year starting at 2)
        
        im_flow = [ 0                                                ;
                    pop_trans(year) * imm_age(age) * legal_rate(1)   ;
                    pop_trans(year) * imm_age(age) * illegal_rate(1) ];
        
        for ik = 1:nk
            for ib = 1:nb
                
                point_k = max(koptss(ik, ib, age - max(T_work, -startyear)), kgrid(1));
                loc1 = find(kgrid(1:nk-1) <= point_k, 1, 'last');
                w1 = (point_k - kgrid(loc1)) / (kgrid(loc1+1) - kgrid(loc1));
                w1 = min(w1, 1);
                
                for ipop = 1:3
                    
                    dist_hold = dist_r(ik,ib,t,ipop,idem) + (ik == 1)*(ib == 1)*im_flow(ipop);
                    
                    dist_r(loc1  , ib, t+1, ipop, idem) = dist_r(loc1  , ib, t+1, ipop, idem) + surv(age)*(1-w1)*dist_hold;
                    dist_r(loc1+1, ib, t+1, ipop, idem) = dist_r(loc1+1, ib, t+1, ipop, idem) + surv(age)*(w1  )*dist_hold; 
                    
                end
                
            end
        end
        
        dist_r(:,:,t+1,3,idem) = (1-amnesty)*dist_r(:,:,t+1,3,idem);
        dist_r(:,:,t+1,2,idem) = dist_r(:,:,t+1,2,idem) + amnesty*dist_r(:,:,t+1,3,idem);
        
        dist_r(:,:,t+1,3,idem) = (1-deportation)*dist_r(:,:,t+1,3,idem);
        
    end
    
    
    Kalive  = zeros(1,T_life);
    Kdead   = zeros(1,T_life);
    Lab     = zeros(1,T_life);
    ELab    = zeros(1,T_life);
    Dist    = zeros(1,T_life);
    Fedit   = zeros(1,T_life);
    SSrev   = zeros(1,T_life);
    SSexp   = zeros(1,T_life);
    Lfp     = zeros(1,T_life);
    SS_base = zeros(1,T_life);
    
    % solving for aggregates by age: region 1
    for t = 1:T_work+startyear
        for ipop = 1:3
            for iz = 1:nz
                for ik = 1:nk
                    for ib = 1:nb
                        age = t - startyear;
                        Kalive (t) = Kalive (t) + dist_1(ik,iz,ib,t,ipop,idem)*kopt   (ik,iz,ib,t);
                        Kdead  (t) = Kdead  (t) + dist_1(ik,iz,ib,t,ipop,idem)*kopt   (ik,iz,ib,t)*(1-surv(age));
                        ELab   (t) = ELab   (t) + dist_1(ik,iz,ib,t,ipop,idem)*labopt (ik,iz,ib,t)*z(iz,age,idem);
                        Lab    (t) = Lab    (t) + dist_1(ik,iz,ib,t,ipop,idem)*labopt (ik,iz,ib,t);
                        Dist   (t) = Dist   (t) + dist_1(ik,iz,ib,t,ipop,idem);
                        Fedit  (t) = Fedit  (t) + dist_1(ik,iz,ib,t,ipop,idem)*fitax  (ik,iz,ib,t);
                        SSrev  (t) = SSrev  (t) + dist_1(ik,iz,ib,t,ipop,idem)*fsstax (ik,iz,ib,t);
                        SSexp  (t) = SSexp  (t) + dist_1(ik,iz,ib,t,ipop,idem)*benopt (ik,iz,ib,t); % (Extra term in head cohorts missing here)
                        Lfp    (t) = Lfp    (t) + dist_1(ik,iz,ib,t,ipop,idem)*(labopt(ik,iz,ib,t) > 0);
                        SS_base(t) = SS_base(t) + dist_1(ik,iz,ib,t,ipop,idem)*ss_base(ik,iz,ib,t);
                    end
                end
            end
        end
    end
    
    % solving for aggregates by age: region 3
    for t = max(T_work+startyear+1, 1):T_life+startyear
        for ipop = 1:3
            for ik = 1:nk
                for ib = 1:nb
                    age = t - startyear;
                    Kalive(t) = Kalive(t) + dist_r(ik,ib,t,ipop,idem)*koptss  (ik,ib,t-max(T_work+startyear, 0));
                    Kdead (t) = Kdead (t) + dist_r(ik,ib,t,ipop,idem)*koptss  (ik,ib,t-max(T_work+startyear, 0))*(1-surv(age));
                    Dist  (t) = Dist  (t) + dist_r(ik,ib,t,ipop,idem);
                    Fedit (t) = Fedit (t) + dist_r(ik,ib,t,ipop,idem)*fitaxss (ik,ib,t-max(T_work+startyear, 0));
                    SSrev (t) = SSrev (t) + dist_r(ik,ib,t,ipop,idem)*fsstaxss(ik,ib,t-max(T_work+startyear, 0));
                    SSexp (t) = SSexp (t) + dist_r(ik,ib,t,ipop,idem)*benoptss(ik,ib,t-max(T_work+startyear, 0));
                end
            end
        end
    end
    
    
    save(fullfile(jobdir, sprintf('transvars_%u_%u_tail_%u.mat', idem, -startyear+1, polno)), ...
         'dist_1', 'dist_r', 'Kalive', 'Kdead', 'ELab', 'Lab', 'Dist', 'Fedit', 'SSrev', 'SSexp', 'Lfp', 'SS_base');
    
end