function [] = solve_trans()



% --- Initialization ---

jobdir = 'Testing';



% -- Parameter loading ---

s = load('params.mat');

T_life  = s.T;
T_model = s.Tss;

nk   = s.nk;
nz   = s.nz;
nb   = s.nb;
ndem = s.ndem;
ks   = s.kgrid;
bs   = s.bgrid;

zs     = s.z;
transz = s.tr_z;

surv = s.surv;
surv(T_life) = 0;

pgr = s.pgr;


s = load('Imm_Data.mat');

illegal_rate = s.illegal_rate;
imm_age      = s.imm_age;


s = load(fullfile(jobdir, 'imm_polparams.mat'));

DISTz_age    = s.proddist_age;
legal_rate   = s.legal_rate;
amnesty      = s.amnesty;
deportation  = s.deportation;



% --- generate_aggregates ---

% Initialize aggregates
series = {'assets', 'beqs', 'labeffs', 'labs', 'lfprs', 'pits', 'ssts', 'bens'};
for o = series, Dynamic.(o{1}) = zeros(1,T_model); end



% --- Solve transition path cohorts ---

for idem = 1:ndem
    
    % Consolidate cohort optimal decision value arrays
    K   = zeros(nk,nz,nb,T_life,T_model);
    LAB = zeros(nk,nz,nb,T_life,T_model);
    B   = zeros(nk,nz,nb,T_life,T_model);
    PIT = zeros(nk,nz,nb,T_life,T_model);
    SST = zeros(nk,nz,nb,T_life,T_model);
    BEN = zeros(nk,nz,nb,T_life,T_model);
    
    for startyear = (-T_life+1):(T_model-1)
        
        T_past   = max(-startyear, 0);
        T_shift  = max(+startyear, 0);
        T_active = min(startyear+T_life, T_model) - T_shift;
        
        opt = load(fullfile('Freeze', 'Cohorts', sprintf('cohort=%+03d_idem=%u.mat', startyear, idem)));
        
        for t = 1:T_active
            
            age  = t + T_past ;
            year = t + T_shift;
            
            K  (:,:,:,age,year) = opt.K  (:,:,:,t);
            LAB(:,:,:,age,year) = opt.LAB(:,:,:,t);
            B  (:,:,:,age,year) = opt.B  (:,:,:,t);
            PIT(:,:,:,age,year) = opt.PIT(:,:,:,t);
            SST(:,:,:,age,year) = opt.SST(:,:,:,t);
            BEN(:,:,:,age,year) = opt.BEN(:,:,:,t);
            
        end
        
    end
    
    
    
    s = load(fullfile(jobdir, sprintf('distvars_%u.mat', idem)));
    
    DIST     = s.dist;
    DIST_age = ones(1,T_life);
    
    DISTeps = Inf;
    DISTtol = -Inf;
    
    year = 1;
    lastyear = T_model;
    
    
    while (DISTeps > DISTtol)
        
        fprintf('Year %3u\n', year);
        
        
        for ik = 1:nk
            for iz = 1:nz
                for ib = 1:nb
                    for age = 1:T_life
                        for ipop = 1:3
                            Dynamic.assets (year) = Dynamic.assets (year) + DIST(ik,iz,ib,age,ipop)*K  (ik,iz,ib,age,year)*(2-surv(age));
                            Dynamic.beqs   (year) = Dynamic.beqs   (year) + DIST(ik,iz,ib,age,ipop)*K  (ik,iz,ib,age,year)*(1-surv(age));
                            Dynamic.labeffs(year) = Dynamic.labeffs(year) + DIST(ik,iz,ib,age,ipop)*LAB(ik,iz,ib,age,year)*zs(iz,age,idem);
                            Dynamic.labs   (year) = Dynamic.labs   (year) + DIST(ik,iz,ib,age,ipop)*LAB(ik,iz,ib,age,year);
                            Dynamic.lfprs  (year) = Dynamic.lfprs  (year) + DIST(ik,iz,ib,age,ipop)*(LAB(ik,iz,ib,age,year) > 0);
                            Dynamic.pits   (year) = Dynamic.pits   (year) + DIST(ik,iz,ib,age,ipop)*PIT(ik,iz,ib,age,year);
                            Dynamic.ssts   (year) = Dynamic.ssts   (year) + DIST(ik,iz,ib,age,ipop)*SST(ik,iz,ib,age,year);
                            Dynamic.bens   (year) = Dynamic.bens   (year) + DIST(ik,iz,ib,age,ipop)*BEN(ik,iz,ib,age,year);
                        end
                    end
                end
            end
        end
        
        
        if (year == lastyear), break, else, year = year + 1; end
        
        
        DIST_next = zeros(nk,nz,nb,T_life,3);
        
        for age = 1:T_life
            
            im_flow = sum(DIST(:))*[ (age == 1)*pgr               ;
                                     imm_age(age)*legal_rate(1)   ;
                                     imm_age(age)*illegal_rate(1) ];
            
            for ipop = 1:3
                DIST_next(1,:,1,age,ipop) = reshape(DISTz_age(:,age,ipop)*im_flow(ipop), [1,nz,1,1,1]);
            end
            
            if (age > 1)
                
                % Extract optimal k and b decision values
                k_t = max(K(:,:,:,age-1,year-1), ks(1));
                b_t = max(B(:,:,:,age-1,year-1), bs(1));
                
                % Find indices of nearest values in ks and bs series
                jk_lt = ones(size(k_t));
                for elem = 1:length(k_t(:))
                    jk_lt(elem) = find(ks(1:end-1) <= k_t(elem), 1, 'last');
                end
                jk_gt = jk_lt + 1;
                
                jb_lt = ones(size(b_t));
                for elem = 1:length(b_t(:))
                    jb_lt(elem) = find(bs(1:end-1) <= b_t(elem), 1, 'last');
                end
                jb_gt = jb_lt + 1;
                
                % Calculate linear weights for nearest values
                wk_lt = max((ks(jk_gt) - k_t) ./ (ks(jk_gt) - ks(jk_lt)), 0);
                wk_gt = 1 - wk_lt;
                
                wb_lt = max((bs(jb_gt) - b_t) ./ (bs(jb_gt) - bs(jb_lt)), 0);
                wb_gt = 1 - wb_lt;
                
                for jz = 1:nz
                    for ipop = 1:3
                        
                        % Apply survival and productivity transformations to cohort distribution from previous year
                        DIST_transz = DIST(:,:,:,age-1,ipop) * surv(age-1) .* repmat(reshape(transz(:,jz), [1,nz,1]), [nk,1,nb]);
                        
                        % Redistribute cohort distribution according to target indices and weights
                        for elem = 1:numel(DIST_transz)
                            DIST_next(jk_lt(elem), jz, jb_lt(elem), age, ipop) = DIST_next(jk_lt(elem), jz, jb_lt(elem), age, ipop) + wk_lt(elem)*wb_lt(elem)*DIST_transz(elem);
                            DIST_next(jk_gt(elem), jz, jb_lt(elem), age, ipop) = DIST_next(jk_gt(elem), jz, jb_lt(elem), age, ipop) + wk_gt(elem)*wb_lt(elem)*DIST_transz(elem);
                            DIST_next(jk_lt(elem), jz, jb_gt(elem), age, ipop) = DIST_next(jk_lt(elem), jz, jb_gt(elem), age, ipop) + wk_lt(elem)*wb_gt(elem)*DIST_transz(elem);
                            DIST_next(jk_gt(elem), jz, jb_gt(elem), age, ipop) = DIST_next(jk_gt(elem), jz, jb_gt(elem), age, ipop) + wk_gt(elem)*wb_gt(elem)*DIST_transz(elem);
                        end
                        
                    end
                end
                
            end
            
            % Increase legal immigrant population for amnesty, maintaining overall distribution over productivity levels
            DISTz_legal = sum(sum(DIST_next(:,:,:,age,2), 1), 3);
            DIST_next(:,:,:,age,2) = DIST_next(:,:,:,age,2) + repmat(sum(amnesty*DIST_next(:,:,:,age,3), 2), [1,nz,1]).*repmat(DISTz_legal, [nk,1,nb])/sum(DISTz_legal);
            
            % Reduce illegal immigrant population for amnesty and deportation
            DIST_next(:,:,:,age,3) = (1-amnesty-deportation)*DIST_next(:,:,:,age,3);
            
        end
        
        DIST_age_next = sum(sum(reshape(DIST_next, [], T_life, 3), 1), 3);
        DISTeps = max(abs(DIST_age_next(2:end)/DIST_age_next(1) - DIST_age(2:end)/DIST_age(1)));
        
        DIST     = DIST_next;
        DIST_age = DIST_age_next;
        
    end
    
end


save(fullfile(jobdir, 'trans_aggregates.mat'), '-struct', 'Dynamic');




%% Testing

trans_aggregates        = load(fullfile(jobdir  , 'trans_aggregates.mat'));
trans_aggregates_freeze = load(fullfile('Freeze', 'trans_aggregates.mat'));

fprintf('trans_aggregates\n');
valuenames = fields(trans_aggregates);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = trans_aggregates.(valuename)(:) - trans_aggregates_freeze.(valuename)(:);
    if any(isnan(delta))
        fprintf('\t%-14sNaN found\n', valuename);
    elseif any(delta)
        pdev = abs(nanmean(delta*2 ./ (trans_aggregates.(valuename)(:) + trans_aggregates_freeze.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');



end