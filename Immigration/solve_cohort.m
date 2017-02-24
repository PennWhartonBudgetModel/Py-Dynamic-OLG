function [Cohort] = solve_cohort(startyear, idem)

% --- Argument verification ---

load('params.mat')
load('Surv_Probs.mat')
load('Imm_Data.mat')

jobdir = 'Testing';
load(fullfile(jobdir, 'imm_polparams.mat'))

load('SSVALS.mat', 'pop_prev')
pops = [pop_prev, pops]; %#ok<NODEF>

load(fullfile(jobdir, 'DIST.mat'));


T_life   = T;
T_work   = Tr;
T_model  = Tss;
T_past   = max(-startyear, 0);
T_shift  = max(+startyear, 0);
T_active = min(startyear+T_life, T_model) - T_shift;

dist = zeros(nk,nz,nb,T_active,3);



% --- Dynamic optimization ---

opt = load(fullfile('Freeze', 'Cohorts', sprintf('cohort=%+03d_idem=%u.mat', startyear, idem)));

K   = opt.K  ;
LAB = opt.LAB;
B   = opt.B  ;
PIT = opt.PIT;
SST = opt.SST;
BEN = opt.BEN;



% --- Distribution generation ---

% --- Initialize distributions ---

if (startyear < 0)
    
    dist(:,:,:,1,:) = DIST(:,:,:,T_past+1,:,idem); %#ok<NODEF>
    
else
    
    age  = 1;
    year = max(1, min(startyear+age, T_model));
    
    im_flow = [ pops(year) * pgr                            ;
                pops(year) * imm_age(age) * legal_rate(1)   ;
                pops(year) * imm_age(age) * illegal_rate(1) ];
    
    for ipop = 1:3
        dist(1,:,1,1,ipop) = reshape(proddist_age(:,1,ipop)*im_flow(ipop), [1,nz,1,1,1]);
    end
    
end



% --- Generate distributions through forward propagation ---

for t = 1:T_active-1
    
    age  = t + T_past;
    year = max(1, min(startyear+age, T_model)) + 1; % (Addition of 1 may be erroneous; also possible that max should be taken with 0 instead)
    
    im_flow = [ 0                                           ;
                pops(year) * imm_age(age) * legal_rate(1)   ;
                pops(year) * imm_age(age) * illegal_rate(1) ];
    
    % Extract optimal k and b values
    k_t = max(K(:,:,:,t), kgrid(1));
    b_t = max(B(:,:,:,t), bgrid(1));
    
    % Find indices of nearest values in ks and bs series
    jk_lt = ones(size(k_t));
    for elem = 1:length(k_t(:))
        jk_lt(elem) = find(kgrid(1:end-1) <= k_t(elem), 1, 'last');
    end
    jk_gt = jk_lt + 1;
    
    jb_lt = ones(size(b_t));
    for elem = 1:length(b_t(:))
        jb_lt(elem) = find(bgrid(1:end-1) <= b_t(elem), 1, 'last');
    end
    jb_gt = jb_lt + 1;
    
    % Calculate linear weights for nearest values
    wk_lt = max((kgrid(jk_gt) - k_t) ./ (kgrid(jk_gt) - kgrid(jk_lt)), 0);
    wk_gt = 1 - wk_lt;
    
    wb_lt = max((bgrid(jb_gt) - b_t) ./ (bgrid(jb_gt) - bgrid(jb_lt)), 0);
    wb_gt = 1 - wb_lt;
    
    for jz = 1:nz
        for ipop = 1:3
            
            DIST_inflow = zeros(nk,nz,nb);
            DIST_inflow(1,:,1) = reshape(proddist_age(:,age,ipop)*im_flow(ipop), [1,nz,1]);
            
            DIST_step = (dist(:,:,:,t,ipop) + DIST_inflow) .* repmat(reshape(tr_z(:,jz), [1,nz,1]), [nk,1,nb]) * surv(age);
            
            for elem = 1:numel(DIST_step)
                dist(jk_lt(elem), jz, jb_lt(elem), t+1, ipop) = dist(jk_lt(elem), jz, jb_lt(elem), t+1, ipop) + wk_lt(elem)*wb_lt(elem)*DIST_step(elem);
                dist(jk_gt(elem), jz, jb_lt(elem), t+1, ipop) = dist(jk_gt(elem), jz, jb_lt(elem), t+1, ipop) + wk_gt(elem)*wb_lt(elem)*DIST_step(elem);
                dist(jk_lt(elem), jz, jb_gt(elem), t+1, ipop) = dist(jk_lt(elem), jz, jb_gt(elem), t+1, ipop) + wk_lt(elem)*wb_gt(elem)*DIST_step(elem);
                dist(jk_gt(elem), jz, jb_gt(elem), t+1, ipop) = dist(jk_gt(elem), jz, jb_gt(elem), t+1, ipop) + wk_gt(elem)*wb_gt(elem)*DIST_step(elem);
            end
            
        end
    end
    
    % Increase legal immigrant population for amnesty, maintaining overall distribution over productivity levels
    distz_legal = sum(sum(dist(:,:,:,t+1,2), 1), 3);
    dist(:,:,:,t+1,2) = dist(:,:,:,t+1,2) + repmat(sum(amnesty*dist(:,:,:,t+1,3), 2), [1,nz,1]).*repmat(distz_legal, [nk,1,nb])/sum(distz_legal);
    
    % Reduce illegal immigrant population for amnesty and deportation
    dist(:,:,:,t+1,3) = (1-amnesty-deportation)*dist(:,:,:,t+1,3);
    
    % if dist_previous(t+1) does not exist, set dist_previous(t+1) to dist(t+1)
    
end



% --- Aggregate generation ---

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
        
        for iz = 1:nz
            for ik = 1:nk
                for ib = 1:nb
                    Cohort.assets (t) = Cohort.assets (t) + dist(ik,iz,ib,t,ipop)*K  (ik,iz,ib,t)*(2-surv(age));
                    Cohort.beqs   (t) = Cohort.beqs   (t) + dist(ik,iz,ib,t,ipop)*K  (ik,iz,ib,t)*(1-surv(age));
                    Cohort.labeffs(t) = Cohort.labeffs(t) + dist(ik,iz,ib,t,ipop)*LAB(ik,iz,ib,t)*z(iz,age,idem);
                    Cohort.labs   (t) = Cohort.labs   (t) + dist(ik,iz,ib,t,ipop)*LAB(ik,iz,ib,t);
                    Cohort.lfprs  (t) = Cohort.lfprs  (t) + dist(ik,iz,ib,t,ipop)*(LAB(ik,iz,ib,t) > 0);
                    Cohort.pits   (t) = Cohort.pits   (t) + dist(ik,iz,ib,t,ipop)*PIT(ik,iz,ib,t);
                    Cohort.ssts   (t) = Cohort.ssts   (t) + dist(ik,iz,ib,t,ipop)*SST(ik,iz,ib,t);
                    Cohort.bens   (t) = Cohort.bens   (t) + dist(ik,iz,ib,t,ipop)*BEN(ik,iz,ib,t);
                end
            end
        end
        
    end
end


end