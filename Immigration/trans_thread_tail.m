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
    
    dist_1(:,:,:,1,:,idem) = dist1(:,:,:,-startyear+1,:,idem); %#ok<NODEF>
    dist_r(:,  :,1,:,idem) = distr(:,  :,-startyear+1,:,idem); %#ok<NODEF>
    
    
    for t = 1:Tr+startyear
        
        age  = t - startyear;
        year = max(1, min(t, Tss)) + 1; % (Addition of 1 may be erroneous; leads to year starting at 2)
        
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
                            
                            if (age < Tr)
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
            
            if (age < Tr)
                
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
    
    
    for t = max(Tr+startyear+1, 1):T-1+startyear
        
        age  = t - startyear;
        year = max(1, min(t, Tss)) + 1; % (Addition of 1 may be erroneous; leads to year starting at 2)
        
        im_flow = [ 0                                                ;
                    pop_trans(year) * imm_age(age) * legal_rate(1)   ;
                    pop_trans(year) * imm_age(age) * illegal_rate(1) ];
        
        for ik = 1:nk
            for ib = 1:nb
                
                point_k = max(koptss(ik, ib, t - startyear -max(Tr, -startyear)), kgrid(1));
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
    
    
    Kalive  = zeros(1,T);
    Kdead   = zeros(1,T);
    Lab     = zeros(1,T);
    ELab    = zeros(1,T);
    Dist    = zeros(1,T);
    Fedit   = zeros(1,T);
    SSrev   = zeros(1,T);
    SSexp   = zeros(1,T);
    Lfp     = zeros(1,T);
    SS_base = zeros(1,T);
    
    % solving for aggregates by age: region 1
    if (-startyear+1 <= Tr)
        for t = -startyear+1:Tr
            for ipop = 1:3
                for k2 = 1:nz
                    for ik = 1:nk
                        for ib = 1:nb
                            Kalive (t+startyear) = Kalive (t+startyear) + dist_1(ik,k2,ib,t+startyear,ipop,idem)*kopt   (ik,k2,ib,t+startyear);
                            Kdead  (t+startyear) = Kdead  (t+startyear) + dist_1(ik,k2,ib,t+startyear,ipop,idem)*kopt   (ik,k2,ib,t+startyear)*(1-surv(t));
                            ELab   (t+startyear) = ELab   (t+startyear) + dist_1(ik,k2,ib,t+startyear,ipop,idem)*labopt (ik,k2,ib,t+startyear)*z(k2,t,idem);
                            Lab    (t+startyear) = Lab    (t+startyear) + dist_1(ik,k2,ib,t+startyear,ipop,idem)*labopt (ik,k2,ib,t+startyear);
                            Dist   (t+startyear) = Dist   (t+startyear) + dist_1(ik,k2,ib,t+startyear,ipop,idem);
                            Fedit  (t+startyear) = Fedit  (t+startyear) + dist_1(ik,k2,ib,t+startyear,ipop,idem)*fitax  (ik,k2,ib,t+startyear);
                            SSrev  (t+startyear) = SSrev  (t+startyear) + dist_1(ik,k2,ib,t+startyear,ipop,idem)*fsstax (ik,k2,ib,t+startyear);
                            SSexp  (t+startyear) = SSexp  (t+startyear) + dist_1(ik,k2,ib,t+startyear,ipop,idem)*benopt (ik,k2,ib,t+startyear);
                            Lfp    (t+startyear) = Lfp    (t+startyear) + dist_1(ik,k2,ib,t+startyear,ipop,idem)*(labopt(ik,k2,ib,t+startyear) > 0);
                            SS_base(t+startyear) = SS_base(t+startyear) + dist_1(ik,k2,ib,t+startyear,ipop,idem)*ss_base(ik,k2,ib,t+startyear);
                        end
                    end
                end
            end
        end
    end
    
    % solving for aggregates by age: region 3
    for t = max(Tr+1, -startyear+1):T
        for ipop = 1:3
            for ik = 1:nk
                for ib = 1:nb
                    Kalive(t+startyear) = Kalive(t+startyear) + dist_r(ik,ib,t+startyear,ipop,idem)*koptss  (ik,ib,t-max(Tr,-startyear));
                    Kdead (t+startyear) = Kdead (t+startyear) + dist_r(ik,ib,t+startyear,ipop,idem)*koptss  (ik,ib,t-max(Tr,-startyear))*(1-surv(t));
                    Dist  (t+startyear) = Dist  (t+startyear) + dist_r(ik,ib,t+startyear,ipop,idem);
                    Fedit (t+startyear) = Fedit (t+startyear) + dist_r(ik,ib,t+startyear,ipop,idem)*fitaxss (ik,ib,t-max(Tr,-startyear));
                    SSrev (t+startyear) = SSrev (t+startyear) + dist_r(ik,ib,t+startyear,ipop,idem)*fsstaxss(ik,ib,t-max(Tr,-startyear));
                    SSexp (t+startyear) = SSexp (t+startyear) + dist_r(ik,ib,t+startyear,ipop,idem)*benoptss(ik,ib,t-max(Tr,-startyear));
                end
            end
        end
    end
    
    
    save(fullfile(jobdir, sprintf('transvars_%u_%u_tail_%u.mat', idem, -startyear+1, polno)), ...
         'dist_1', 'dist_r', 'Kalive', 'Kdead', 'ELab', 'Lab', 'Dist', 'Fedit', 'SSrev', 'SSexp', 'Lfp', 'SS_base');
    
end