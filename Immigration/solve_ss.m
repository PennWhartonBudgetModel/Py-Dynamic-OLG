function [] = solve_ss()


jobdir = 'Testing';



s = load('params.mat');

T_life  = s.T;
T_model = 1;

nk   = s.nk;
nz   = s.nz;
nb   = s.nb;
ndem = s.ndem;
ks   = s.kgrid;
bs   = s.bgrid;

zs     = s.z;
transz = s.tr_z;

surv = s.surv(1:T_life);
surv(T_life) = 0;

pgr = s.pgr;

DISTz_age = s.proddist_age;


s = load('Imm_Data.mat');

legal_rate   = s.legal_rate(1);
illegal_rate = s.illegal_rate(1);
imm_age      = s.imm_age;

amnesty     = 0.00;
deportation = 0.00;



% Initialize aggregates
series = {'assets', 'beqs', 'labeffs', 'labs', 'lfprs', 'pits', 'ssts', 'bens'};
for a = series, Dynamic.(a{1}) = []; end


for idem = 1:ndem
    
    
    % Solve dynamic optimization for representative cohort
    filename = sprintf('sspol%u.mat', idem);
    copyfile(fullfile('Freeze', filename), fullfile(jobdir, filename));
    
    s = load(fullfile(jobdir, sprintf('sspol%u.mat', idem)));
    
    K   = s.kopt  ;
    LAB = s.labopt;
    B   = s.bopt  ;
    PIT = zeros(nk,nz,nb,T_life,T_model);
    SST = zeros(nk,nz,nb,T_life,T_model);
    BEN = zeros(nk,nz,nb,T_life,T_model);
    
    
    
    % Initialize distribution
    DIST = ones(nk,nz,nb,T_life,3);
    DIST = DIST / numel(DIST);
    DIST_age = ones(1,T_life);
    
    DISTeps = Inf;
    DISTtol = 1e-4;
    
    year = 1;
    lastyear = Inf;
    
    
    while (DISTeps > DISTtol)
        
        fprintf('Year %3u\n', year);
        
        
        % Add values to aggregates for current year
        for a = series, if (length(Dynamic.(a{1})) < year), Dynamic.(a{1})(year) = 0; end, end
        for ipop = 1:3
            A.assets  = DIST(:,:,:,:,ipop).*K  (:,:,:,:,min(year, T_model)).*repmat(reshape(2-surv, [1,1,1,T_life]), [nk,nz,nb,1]);
            A.beqs    = DIST(:,:,:,:,ipop).*K  (:,:,:,:,min(year, T_model)).*repmat(reshape(1-surv, [1,1,1,T_life]), [nk,nz,nb,1]);
            A.labeffs = DIST(:,:,:,:,ipop).*LAB(:,:,:,:,min(year, T_model)).*repmat(reshape(zs(:,:,idem), [1,nz,1,T_life]), [nk,1,nb,1]);
            A.labs    = DIST(:,:,:,:,ipop).*LAB(:,:,:,:,min(year, T_model));
            A.lfprs   = DIST(:,:,:,:,ipop).*(LAB(:,:,:,:,min(year, T_model)) > 0);
            A.pits    = DIST(:,:,:,:,ipop).*PIT(:,:,:,:,min(year, T_model));
            A.ssts    = DIST(:,:,:,:,ipop).*SST(:,:,:,:,min(year, T_model));
            A.bens    = DIST(:,:,:,:,ipop).*BEN(:,:,:,:,min(year, T_model));
            for a = series, Dynamic.(a{1})(year) = Dynamic.(a{1})(year) + sum(A.(a{1})(:)); end
        end
        
        
        if (year < lastyear), year = year + 1; else, break, end
        
        
        % Initialize distribution for next year
        DIST_next = zeros(nk,nz,nb,T_life,3);
        
        % Grow populations
        DIST_total = sum(DIST(:));
        DIST_next(1,:,1,1,1) = DIST_total * reshape(DISTz_age(:,1,1), [1,nz,1,1]) * pgr;
        DIST_next(1,:,1,:,2) = DIST_total * reshape(DISTz_age(:,:,2), [1,nz,1,T_life]) .* repmat(reshape(imm_age, [1,1,1,T_life]), [1,nz,1,1]) * legal_rate;
        DIST_next(1,:,1,:,3) = DIST_total * reshape(DISTz_age(:,:,3), [1,nz,1,T_life]) .* repmat(reshape(imm_age, [1,1,1,T_life]), [1,nz,1,1]) * illegal_rate;
        
        for age = 2:T_life
            
            % Extract optimal k and b decision values
            k_t = max(K(:,:,:,age-1,min(year-1, T_model)), ks(1));
            b_t = max(B(:,:,:,age-1,min(year-1, T_model)), bs(1));
            
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
                    
                    % Apply survival and productivity transformations to cohort distribution from current year
                    DIST_transz = DIST(:,:,:,age-1,ipop) * surv(age-1) .* repmat(reshape(transz(:,jz), [1,nz,1]), [nk,1,nb]);
                    
                    % Redistribute cohort for next year according to target indices and weights
                    for elem = 1:numel(DIST_transz)
                        DIST_next(jk_lt(elem), jz, jb_lt(elem), age, ipop) = DIST_next(jk_lt(elem), jz, jb_lt(elem), age, ipop) + wk_lt(elem)*wb_lt(elem)*DIST_transz(elem);
                        DIST_next(jk_gt(elem), jz, jb_lt(elem), age, ipop) = DIST_next(jk_gt(elem), jz, jb_lt(elem), age, ipop) + wk_gt(elem)*wb_lt(elem)*DIST_transz(elem);
                        DIST_next(jk_lt(elem), jz, jb_gt(elem), age, ipop) = DIST_next(jk_lt(elem), jz, jb_gt(elem), age, ipop) + wk_lt(elem)*wb_gt(elem)*DIST_transz(elem);
                        DIST_next(jk_gt(elem), jz, jb_gt(elem), age, ipop) = DIST_next(jk_gt(elem), jz, jb_gt(elem), age, ipop) + wk_gt(elem)*wb_gt(elem)*DIST_transz(elem);
                    end
                    
                end
            end
            
        end
        
        % Increase legal immigrant population for amnesty, maintaining productivity distributions
        DISTz_legal = DIST_next(:,:,:,:,2) ./ repmat(sum(DIST_next(:,:,:,:,2), 2), [1,nz,1,1]);
        DISTz_legal(isnan(DISTz_legal)) = 1/nz;
        
        DIST_next(:,:,:,:,2) = DIST_next(:,:,:,:,2) + repmat(sum(amnesty*DIST_next(:,:,:,:,3), 2), [1,nz,1,1]).*DISTz_legal;
        
        % Reduce illegal immigrant population for amnesty and deportation
        DIST_next(:,:,:,:,3) = (1-amnesty-deportation)*DIST_next(:,:,:,:,3);
        
        % Calculate distribution convergence error
        DIST_age_next = sum(sum(reshape(DIST_next, [], T_life, 3), 1), 3);
        DIST_age_next = DIST_age_next / DIST_age_next(1);
        DISTeps = max(abs(DIST_age_next - DIST_age));
        
        % Update distributions
        DIST = DIST_next;
        DIST_age = DIST_age_next;
        
    end
    
    % Save distribution
    save(fullfile(jobdir, sprintf('distvars_%u.mat', idem)), 'DIST', 'DIST_age');
    
end


save(fullfile(jobdir, 'ss_aggregates.mat'), '-struct', 'Dynamic');




%% Testing

distvars_1        = load(fullfile(jobdir  , 'distvars_1.mat'));
distvars_1_freeze = load(fullfile('Freeze', 'distvars_1.mat'));

fprintf('distvars_1\n');
valuenames = fields(distvars_1);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = distvars_1.(valuename)(:) - distvars_1_freeze.(valuename)(:);
    if any(isnan(delta))
        fprintf('\t%-14sNaN found\n', valuename);
    elseif any(delta)
        pdev = abs(nanmean(delta*2 ./ (distvars_1.(valuename)(:) + distvars_1_freeze.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');


distvars_2        = load(fullfile(jobdir  , 'distvars_2.mat'));
distvars_2_freeze = load(fullfile('Freeze', 'distvars_2.mat'));

fprintf('distvars_2\n');
valuenames = fields(distvars_2);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = distvars_2.(valuename)(:) - distvars_2_freeze.(valuename)(:);
    if any(isnan(delta))
        fprintf('\t%-14sNaN found\n', valuename);
    elseif any(delta)
        pdev = abs(nanmean(delta*2 ./ (distvars_2.(valuename)(:) + distvars_2_freeze.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');


ss_aggregates        = load(fullfile(jobdir  , 'ss_aggregates.mat'));
ss_aggregates_freeze = load(fullfile('Freeze', 'ss_aggregates.mat'));

fprintf('ss_aggregates\n');
valuenames = fields(ss_aggregates);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = ss_aggregates.(valuename)(:) - ss_aggregates_freeze.(valuename)(:);
    if any(isnan(delta))
        fprintf('\t%-14sNaN found\n', valuename);
    elseif any(delta)
        pdev = abs(nanmean(delta*2 ./ (ss_aggregates.(valuename)(:) + ss_aggregates_freeze.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');


end