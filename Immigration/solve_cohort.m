function [Cohort] = solve_cohort(startyear, idem)

load('params.mat')
load('Surv_Probs.mat')
load('Imm_Data.mat')

jobdir = 'Testing';
load(fullfile(jobdir, sprintf('imm_polparams_1.mat')))

load('SSVALS.mat', 'pop_prev')
pop_trans = [pop_prev; pop_trans]; %#ok<NODEF>

load(fullfile(jobdir, 'eqmdist.mat'));


T_life   = T;
T_work   = Tr;
T_model  = Tss;
T_past   = max(-startyear, 0);
T_shift  = max(+startyear, 0);
T_active = min(startyear+T_life, T_model) - T_shift;

dist_1 = zeros(nk,nz,nb,T_active,3);
dist_r = zeros(nk,   nb,T_active,3);


opt = load(fullfile('Freeze', 'Cohorts', sprintf('cohort=%+03d_idem=%u.mat', startyear, idem)));

K   = opt.K  ;
LAB = opt.LAB;
B   = opt.B  ;
PIT = opt.PIT;
SST = opt.SST;
BEN = opt.BEN;


if (startyear < 0)
    
    dist_1(:,:,:,1,:) = dist1(:,:,:,T_past+1,:,idem); %#ok<NODEF>
    dist_r(:,  :,1,:) = distr(:,  :,T_past+1,:,idem); %#ok<NODEF>
    
else
    
    age  = 1;
    year = max(1, min(startyear+age, T_model));
    
    im_flow = [ pop_trans(year) * pgr                            ;
                pop_trans(year) * imm_age(age) * legal_rate(1)   ;
                pop_trans(year) * imm_age(age) * illegal_rate(1) ];
    
    for iz = 1:nz
        for ipop = 1:3
            dist_1(1,iz,1,1,ipop) = proddist_age(iz,1,ipop) * im_flow(ipop);
        end
    end
    
end


for t = 1:T_active-1
    
    age  = t + T_past;
    year = max(1, min(startyear+age, T_model)) + 1; % (Addition of 1 may be erroneous; leads to year starting at 2)
    
    im_flow = [ 0                                                ;
                pop_trans(year) * imm_age(age) * legal_rate(1)   ;
                pop_trans(year) * imm_age(age) * illegal_rate(1) ];
    
    if (age <= T_work)
        
        for iz = 1:nz
            for ik = 1:nk
                for ib = 1:nb
                    
                    point_k = max(K(ik,iz,ib,t), kgrid(1));
                    loc1 = find(kgrid(1:nk-1) <= point_k, 1, 'last');
                    w1 = (point_k - kgrid(loc1)) / (kgrid(loc1+1) - kgrid(loc1));
                    w1 = min(w1, 1);
                    
                    point_b = max(B(ik,iz,ib,t), bgrid(1));
                    loc2 = find(bgrid(1:nb-1) <= point_b, 1, 'last');
                    w2 = (point_b - bgrid(loc2)) / (bgrid(loc2+1) - bgrid(loc2));
                    w2 = min(w2, 1);
                    
                    for jz = 1:nz
                        for ipop = 1:3
                            
                            dist_hold = dist_1(ik,iz,ib,t,ipop) + (ik == 1)*(ib == 1)*proddist_age(iz,age,ipop)*im_flow(ipop);
                            
                            if (age < T_work)
                                dist_1(loc1  ,jz,loc2  ,t+1,ipop) = dist_1(loc1  ,jz,loc2  ,t+1,ipop) + surv(age)*(1-w2)*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_1(loc1+1,jz,loc2  ,t+1,ipop) = dist_1(loc1+1,jz,loc2  ,t+1,ipop) + surv(age)*(1-w2)*(w1  )*tr_z(iz,jz)*dist_hold;
                                dist_1(loc1  ,jz,loc2+1,t+1,ipop) = dist_1(loc1  ,jz,loc2+1,t+1,ipop) + surv(age)*(w2  )*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_1(loc1+1,jz,loc2+1,t+1,ipop) = dist_1(loc1+1,jz,loc2+1,t+1,ipop) + surv(age)*(w2  )*(w1  )*tr_z(iz,jz)*dist_hold;
                            else
                                dist_r(loc1  ,   loc2  ,t+1,ipop) = dist_r(loc1  ,   loc2  ,t+1,ipop) + surv(age)*(1-w2)*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_r(loc1+1,   loc2  ,t+1,ipop) = dist_r(loc1+1,   loc2  ,t+1,ipop) + surv(age)*(1-w2)*(w1  )*tr_z(iz,jz)*dist_hold;
                                dist_r(loc1  ,   loc2+1,t+1,ipop) = dist_r(loc1  ,   loc2+1,t+1,ipop) + surv(age)*(w2  )*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_r(loc1+1,   loc2+1,t+1,ipop) = dist_r(loc1+1,   loc2+1,t+1,ipop) + surv(age)*(w2  )*(w1  )*tr_z(iz,jz)*dist_hold;
                            end
                            
                        end
                    end
                    
                end
            end
        end
        
        if (age < T_work)
            
            amnesty_dist = squeeze(amnesty*dist_1(:,:,:,t+1,3)); 
            amnesty_dist = permute(amnesty_dist, [2,1,3]);
            amnesty_dist = squeeze(sum(amnesty_dist));
            
            prod_legal = squeeze(dist_1(:,:,:,t+1,2));
            prod_legal = squeeze(sum(sum(permute(prod_legal, [1,3,2]))));
            prod_legal = prod_legal / sum(prod_legal);
            
            for ik = 1:nk
                for ib = 1:nb
                    for iz = 1:nz
                        dist_1(ik,iz,ib,t+1,2) = dist_1(ik,iz,ib,t+1,2) + prod_legal(iz)*amnesty_dist(ik,ib);
                    end
                end
            end
            dist_1(:,:,:,t+1,3) = (1-amnesty)*dist_1(:,:,:,t+1,3);
            
        else
            
            dist_r(:,:,t+1,2) = dist_r(:,:,t+1,2) + amnesty*dist_r(:,:,t+1,3);
            dist_r(:,:,t+1,3) = (1-amnesty)*dist_r(:,:,t+1,3);
            
        end
        
        dist_1(:,:,:,t+1,3) = (1-deportation)*dist_1(:,:,:,t+1,3);
        dist_r(:,  :,t+1,3) = (1-deportation)*dist_r(:,  :,t+1,3);
        
    else
        
        for ik = 1:nk
            for ib = 1:nb
                
                point_k = max(K(ik,1,ib,t), kgrid(1));
                loc1 = find(kgrid(1:nk-1) <= point_k, 1, 'last');
                w1 = (point_k - kgrid(loc1)) / (kgrid(loc1+1) - kgrid(loc1));
                w1 = min(w1, 1);
                
                for ipop = 1:3
                    
                    dist_hold = dist_r(ik,ib,t,ipop) + (ik == 1)*(ib == 1)*im_flow(ipop);
                    
                    dist_r(loc1  ,ib,t+1,ipop) = dist_r(loc1  ,ib,t+1,ipop) + surv(age)*(1-w1)*dist_hold;
                    dist_r(loc1+1,ib,t+1,ipop) = dist_r(loc1+1,ib,t+1,ipop) + surv(age)*(w1  )*dist_hold;
                    
                end
                
            end
        end
        
        dist_r(:,:,t+1,2) = dist_r(:,:,t+1,2) + amnesty*dist_r(:,:,t+1,3);
        dist_r(:,:,t+1,3) = (1-amnesty)*dist_r(:,:,t+1,3);
        
        dist_r(:,:,t+1,3) = (1-deportation)*dist_r(:,:,t+1,3);
        
    end
    
end




Cohort.assets  = zeros(1,T_active);
Cohort.beqs    = zeros(1,T_active);
Cohort.labeffs = zeros(1,T_active);
Cohort.labs    = zeros(1,T_active);
Cohort.lfprs   = zeros(1,T_active);
Cohort.pits    = zeros(1,T_active);
Cohort.ssts    = zeros(1,T_active);
Cohort.bens    = zeros(1,T_active);

for t = 1:T_active
    for ipop = 1:3
        
        age = t + T_past;
        
        if (age <= T_work)
            
            for iz = 1:nz
                for ik = 1:nk
                    for ib = 1:nb
                        Cohort.assets (t) = Cohort.assets (t) + dist_1(ik,iz,ib,t,ipop)*K  (ik,iz,ib,t)*(2-surv(age));
                        Cohort.beqs   (t) = Cohort.beqs   (t) + dist_1(ik,iz,ib,t,ipop)*K  (ik,iz,ib,t)*(1-surv(age));
                        Cohort.labeffs(t) = Cohort.labeffs(t) + dist_1(ik,iz,ib,t,ipop)*LAB(ik,iz,ib,t)*z(iz,age,idem);
                        Cohort.labs   (t) = Cohort.labs   (t) + dist_1(ik,iz,ib,t,ipop)*LAB(ik,iz,ib,t);
                        Cohort.lfprs  (t) = Cohort.lfprs  (t) + dist_1(ik,iz,ib,t,ipop)*(LAB(ik,iz,ib,t) > 0);
                        Cohort.pits   (t) = Cohort.pits   (t) + dist_1(ik,iz,ib,t,ipop)*PIT(ik,iz,ib,t);
                        Cohort.ssts   (t) = Cohort.ssts   (t) + dist_1(ik,iz,ib,t,ipop)*SST(ik,iz,ib,t);
                        Cohort.bens   (t) = Cohort.bens   (t) + dist_1(ik,iz,ib,t,ipop)*BEN(ik,iz,ib,t);
                    end
                end
            end
            
        else
            
            for ik = 1:nk
                for ib = 1:nb
                    Cohort.assets(t) = Cohort.assets(t) + dist_r(ik,ib,t,ipop)*K  (ik,1,ib,t)*(2-surv(age));
                    Cohort.beqs  (t) = Cohort.beqs  (t) + dist_r(ik,ib,t,ipop)*K  (ik,1,ib,t)*(1-surv(age));
                    Cohort.pits  (t) = Cohort.pits  (t) + dist_r(ik,ib,t,ipop)*PIT(ik,1,ib,t);
                    Cohort.ssts  (t) = Cohort.ssts  (t) + dist_r(ik,ib,t,ipop)*SST(ik,1,ib,t);
                    Cohort.bens  (t) = Cohort.bens  (t) + dist_r(ik,ib,t,ipop)*BEN(ik,1,ib,t);
                end
            end
            
        end
        
    end
end


end