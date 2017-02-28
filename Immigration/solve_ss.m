function [] = solve_ss()



% --- Initialization ---

jobdir = 'Testing';



% --- Parameter loading ---

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
for o = series, Dynamic.(o{1}) = zeros(1,T_model); end


for idem = 1:ndem
    
    
    % Solve dynamic optimization for representative cohort
    filename = sprintf('sspol%u.mat', idem);
    copyfile(fullfile('Freeze', filename), fullfile(jobdir, filename));
    
    s = load(fullfile(jobdir, sprintf('sspol%u.mat', idem)));
    
    K   = s.kopt  ;
    LAB = s.labopt;
    B   = s.bopt  ;
    
    
    
    % Initialize distribution
    DIST     = zeros(nk,nz,nb,T_life,3);
    DIST_age = ones(1,T_life);
    
    DISTeps  = Inf;
    DISTtol  = 1e-4;
    
    pops(1) = 1;
    
    year = 1;
    lastyear = Inf;
    
    
    while (DISTeps > DISTtol && year < lastyear)
        
        fprintf('Year %3u\n', year);
        
        
        
        % --- Initialize distributions ---
        
        DIST_next = zeros(nk,nz,nb,T_life,3);
        
        im_flow = [ pops(year) * pgr                       ;
                    pops(year) * imm_age(1) * legal_rate   ;
                    pops(year) * imm_age(1) * illegal_rate ];
        
        for ipop = 1:3
            DIST_next(1,:,1,1,ipop) = reshape(DISTz_age(:,1,ipop)*im_flow(ipop), [1,nz,1,1,1]);
        end
        
        
        
        % --- Generate distributions through forward propagation ---
        
        for age = 2:T_life
            
            im_flow = [ 0                                          ;
                        pops(year) * imm_age(age-1) * legal_rate   ;
                        pops(year) * imm_age(age-1) * illegal_rate ];
            
            for ipop = 1:3
                DIST_next(1,:,1,age,ipop) = reshape(DISTz_age(:,age,ipop)*im_flow(ipop), [1,nz,1,1,1]);
            end
            
            
            
            % Extract optimal k and b values
            k_t = max(K(:,:,:,age-1,min(year, T_model)), ks(1));
            b_t = max(B(:,:,:,age-1,min(year, T_model)), bs(1));
            
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
                    
                    DIST_transz = DIST(:,:,:,age-1,ipop) * surv(age-1) .* repmat(reshape(transz(:,jz), [1,nz,1]), [nk,1,nb]);
                    
                    for elem = 1:numel(DIST_transz)
                        DIST_next(jk_lt(elem), jz, jb_lt(elem), age, ipop) = DIST_next(jk_lt(elem), jz, jb_lt(elem), age, ipop) + wk_lt(elem)*wb_lt(elem)*DIST_transz(elem);
                        DIST_next(jk_gt(elem), jz, jb_lt(elem), age, ipop) = DIST_next(jk_gt(elem), jz, jb_lt(elem), age, ipop) + wk_gt(elem)*wb_lt(elem)*DIST_transz(elem);
                        DIST_next(jk_lt(elem), jz, jb_gt(elem), age, ipop) = DIST_next(jk_lt(elem), jz, jb_gt(elem), age, ipop) + wk_lt(elem)*wb_gt(elem)*DIST_transz(elem);
                        DIST_next(jk_gt(elem), jz, jb_gt(elem), age, ipop) = DIST_next(jk_gt(elem), jz, jb_gt(elem), age, ipop) + wk_gt(elem)*wb_gt(elem)*DIST_transz(elem);
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
        DIST_age = DIST_age_next;
        
        pops(year+1) = sum(DIST_next(:)); %#ok<AGROW>
        DIST = DIST_next;
        year = year + 1;
        
        
        
        Dynamic.assets (min(year, T_model)) = 0;
        Dynamic.labeffs(min(year, T_model)) = 0;
        
        for ipop = 1:3
            for age = 1:T_life
                for iz = 1:nz
                    for ik = 1:nk
                        for ib = 1:nb
                            Dynamic.assets (min(year, T_model)) = Dynamic.assets (min(year, T_model)) + DIST_next(ik,iz,ib,age,ipop)*K  (ik,iz,ib,age)*(2-surv(age));
                            Dynamic.labeffs(min(year, T_model)) = Dynamic.labeffs(min(year, T_model)) + DIST_next(ik,iz,ib,age,ipop)*LAB(ik,iz,ib,age)*zs(iz,age,idem);
                        end
                    end
                end
            end
        end
        
    end
    
    
    KPR  = Dynamic.assets ;
    ELAB = Dynamic.labeffs;
    
    save(fullfile(jobdir, sprintf('distvars_%u.mat', idem)), ...
         'DIST', 'KPR', 'ELAB');
    
end




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



end