function [] = trans_thread_head(startyear, polno, impolno)

load('params.mat')
load('Surv_Probs.mat')
load('Imm_Data.mat')

jobdir = 'Testing';
load(fullfile(jobdir, sprintf('imm_polparams_%u.mat', impolno)))

dist_1 = zeros(nk,nz,nb,min(T+1,Tss+1),3,ndem);
dist_r = zeros(nk,   nb,min(T+1,Tss+1),3,ndem);

load('SSVALS.mat', 'pop_prev')
pop_trans = [pop_prev; pop_trans]; %#ok<NODEF>


for idem = 1:ndem
    
    load(fullfile('Freeze', 'Cohorts', sprintf('head%u_%u_%u.mat', startyear+1, idem, polno)));
    
    year = max(1, min(startyear+1, Tss));
    
    % using period 1 imm rate values for steady state
    im_flow = [ pop_trans(year) * pgr                          ;
                pop_trans(year) * imm_age(1) * legal_rate(1)   ;
                pop_trans(year) * imm_age(1) * illegal_rate(1) ];
    
    for iz = 1:nz
        for ipop = 1:3
            dist_1(1,iz,1,1,ipop,idem) = proddist_age(iz,1,ipop) * im_flow(ipop);
        end
    end
    
    for t = 1:min(Tr, Tss-startyear-1)
        
        age  = t;
        year = max(1, min(startyear+age, Tss)) + 1;
        
        % using period 1 imm rate values for steady state
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
                                
                            dist_hold = dist_1(ik,iz,ib,t,ipop,idem) + (ik == 1)*(ib == 1)*proddist_age(iz,t,ipop)*im_flow(ipop);
                        
                            if (t < Tr)
                                dist_1(loc1  , jz, loc2  , t+1, ipop, idem) = dist_1(loc1  , jz, loc2  , t+1, ipop, idem) + surv(t)*(1-w2)*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_1(loc1+1, jz, loc2  , t+1, ipop, idem) = dist_1(loc1+1, jz, loc2  , t+1, ipop, idem) + surv(t)*(1-w2)*(w1  )*tr_z(iz,jz)*dist_hold;
                                dist_1(loc1  , jz, loc2+1, t+1, ipop, idem) = dist_1(loc1  , jz, loc2+1, t+1, ipop, idem) + surv(t)*(w2  )*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_1(loc1+1, jz, loc2+1, t+1, ipop, idem) = dist_1(loc1+1, jz, loc2+1, t+1, ipop, idem) + surv(t)*(w2  )*(w1  )*tr_z(iz,jz)*dist_hold;
                            else
                                dist_r(loc1  ,     loc2  , t+1, ipop, idem) = dist_r(loc1  ,     loc2  , t+1, ipop, idem) + surv(t)*(1-w2)*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_r(loc1+1,     loc2  , t+1, ipop, idem) = dist_r(loc1+1,     loc2  , t+1, ipop, idem) + surv(t)*(1-w2)*(w1  )*tr_z(iz,jz)*dist_hold;
                                dist_r(loc1  ,     loc2+1, t+1, ipop, idem) = dist_r(loc1  ,     loc2+1, t+1, ipop, idem) + surv(t)*(w2  )*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist_r(loc1+1,     loc2+1, t+1, ipop, idem) = dist_r(loc1+1,     loc2+1, t+1, ipop, idem) + surv(t)*(w2  )*(w1  )*tr_z(iz,jz)*dist_hold;
                            end
                            
                        end
                    end
                    
                end
            end
        end
        
        if (amnesty) && (t < Tr)
            
            amnesty_dist = squeeze(amnesty*dist_1(:,:,:,t+1,3,idem)); 
            amnesty_dist = permute(amnesty_dist, [2,1,3]);
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
            
        elseif (amnesty) && (t == Tr)
            
            dist_r(:,:,t+1,2,idem) = dist_r(:,:,t+1,2,idem) + amnesty*dist_r(:,:,t+1,3,idem);
            dist_r(:,:,t+1,3,idem) = (1-amnesty)*dist_r(:,:,t+1,3,idem);
            
        end
        
        dist_1(:,:,:,t+1,3,idem) = (1-deportation)*dist_1(:,:,:,t+1,3,idem);  
        dist_r(:,  :,t+1,3,idem) = (1-deportation)*dist_r(:,  :,t+1,3,idem);
        
    end
    
    
    % looping over dist_r
    if (Tss-startyear > Tr)
        
        for t = Tr+1:min(T-1, Tss-startyear)
            
            age  = t;
            year = max(1, min(startyear+age, Tss)) + 1;
            
            % using period 1 imm rate values for steady state
            im_flow = [ 0                                                ;
                        pop_trans(year) * imm_age(age) * legal_rate(1)   ;
                        pop_trans(year) * imm_age(age) * illegal_rate(1) ];
            
            for ik = 1:nk
                for ib = 1:nb
                    
                    point_k = max(koptss(ik,ib,t-Tr), kgrid(1));
                    loc1 = find(kgrid(1:nk-1) <= point_k, 1, 'last');
                    w1 = (point_k - kgrid(loc1)) / (kgrid(loc1+1) - kgrid(loc1));
                    w1 = min(w1, 1);
                    
                    for ipop = 1:3
                        
                        dist_hold = dist_r(ik,ib,t,ipop,idem) + (ik == 1)*(ib == 1)*im_flow(ipop);
                        
                        dist_r(loc1  , ib, t+1, ipop, idem) = dist_r(loc1  , ib, t+1, ipop, idem) + surv(t)*(1-w1)*dist_hold;
                        dist_r(loc1+1, ib, t+1, ipop, idem) = dist_r(loc1+1, ib, t+1, ipop, idem) + surv(t)*(w1  )*dist_hold;
                        
                    end
                    
                end
            end
            
            dist_r(:,:,t+1,2,idem) = dist_r(:,:,t+1,2,idem) + amnesty*dist_r(:,:,t+1,3,idem);
            dist_r(:,:,t+1,3,idem) = (1-amnesty)*dist_r(:,:,t+1,3,idem);
            
            dist_r(:,:,t+1,3,idem) = (1-deportation)*dist_r(:,:,t+1,3,idem);
            
        end
        
    end
    
    
    Kalive  = zeros(1,min(Tss-startyear,T));
    Kdead   = zeros(1,min(Tss-startyear,T));
    Lab     = zeros(1,min(Tss-startyear,T));
    ELab    = zeros(1,min(Tss-startyear,T));
    Dist    = zeros(1,min(Tss-startyear,T));
    Fedit   = zeros(1,min(Tss-startyear,T));
    SSrev   = zeros(1,min(Tss-startyear,T));
    SSexp   = zeros(1,min(Tss-startyear,T));
    Lfp     = zeros(1,min(Tss-startyear,T));
    SS_base = zeros(1,min(Tss-startyear,T));
    
    % solving for aggregates by age
    for t = 1:min(Tss-startyear, Tr)
        for ipop = 1:3
            for iz = 1:nz
                for ik = 1:nk
                    for ib = 1:nb
                        Kalive (t) = Kalive (t) + dist_1(ik,iz,ib,t,ipop,idem)*kopt   (ik,iz,ib,t);
                        Kdead  (t) = Kdead  (t) + dist_1(ik,iz,ib,t,ipop,idem)*kopt   (ik,iz,ib,t)*(1-surv(t));
                        ELab   (t) = ELab   (t) + dist_1(ik,iz,ib,t,ipop,idem)*labopt (ik,iz,ib,t)*z(iz,t,idem);
                        Lab    (t) = Lab    (t) + dist_1(ik,iz,ib,t,ipop,idem)*labopt (ik,iz,ib,t);
                        Dist   (t) = Dist   (t) + dist_1(ik,iz,ib,t,ipop,idem);
                        Fedit  (t) = Fedit  (t) + dist_1(ik,iz,ib,t,ipop,idem)*fitax  (ik,iz,ib,t);
                        SSrev  (t) = SSrev  (t) + dist_1(ik,iz,ib,t,ipop,idem)*fsstax (ik,iz,ib,t);
                        SSexp  (t) = SSexp  (t) + dist_1(ik,iz,ib,t,ipop,idem)*benopt (ik,iz,ib,t)*max(ss_select(ik,iz,ib,t)-1, 0);
                        Lfp    (t) = Lfp    (t) + dist_1(ik,iz,ib,t,ipop,idem)*(labopt(ik,iz,ib,t) > 0);
                        SS_base(t) = SS_base(t) + dist_1(ik,iz,ib,t,ipop,idem)*ss_base(ik,iz,ib,t);
                    end
                end
            end
        end
    end

    for t = Tr+1:min(T, Tss-startyear)
        for ipop = 1:3
            for ik = 1:nk
                for ib = 1:nb
                    Kalive(t) = Kalive(t) + dist_r(ik,ib,t,ipop,idem)*koptss  (ik,ib,t-Tr);
                    Kdead (t) = Kdead (t) + dist_r(ik,ib,t,ipop,idem)*koptss  (ik,ib,t-Tr)*(1-surv(t));
                    Dist  (t) = Dist  (t) + dist_r(ik,ib,t,ipop,idem);
                    Fedit (t) = Fedit (t) + dist_r(ik,ib,t,ipop,idem)*fitaxss (ik,ib,t-Tr);
                    SSrev (t) = SSrev (t) + dist_r(ik,ib,t,ipop,idem)*fsstaxss(ik,ib,t-Tr);
                    SSexp (t) = SSexp (t) + dist_r(ik,ib,t,ipop,idem)*benoptss(ik,ib,t-Tr);
                end
            end
        end
    end
    
    
    save(fullfile(jobdir, sprintf('transvars_%u_%u_head_%u.mat', idem, startyear+1, polno)), ...
         'dist_1', 'dist_r', 'Kalive', 'Kdead', 'ELab', 'Lab', 'Dist', 'Fedit', 'SSrev', 'SSexp', 'Lfp', 'SS_base');
    
end