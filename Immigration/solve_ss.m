function [] = solve_ss()



% --- Initialization ---

jobdir = 'Testing';
if exist(jobdir, 'dir'), rmdir(jobdir, 's'), end
mkdir(jobdir);



% --- Parameter loading ---

load('params.mat')
surv(T) = 0;

T_life = T;

load('Imm_Data.mat')

amnesty     = 0.00;
deportation = 0.00;



% --- Dynamic aggregate generation ---

rhosseps = Inf;
rhosstol = 1e-4;



load('SSVALS.mat')   % loads the last set of outputs as the first guess



while (rhosseps > rhosstol)
    
    
    
    % --- Define prices ---
    
    DD = DEBTss;  % debt used in the calculation of prices (zero for open economy)
    
    
    
    % --- generate_aggregates ---
    
    for idem = 1:ndem
        
        
        
        % --- solve_cohort dynamic optimization ---
        
        filename = sprintf('sspol%u.mat', idem);
        copyfile(fullfile('Freeze', filename), fullfile(jobdir, filename));
        
        
        
        % --- solve_cohort distribution generation ---
        
        s = load(fullfile(jobdir, sprintf('sspol%u.mat', idem)));
        K   = s.kopt  ;
        B   = s.bopt  ;
        LAB = s.labopt;
        
        dist_previous     = zeros(nk,nz,nb,T_life,3);
        dist_age_previous = ones(1,T_life);
        
        disteps  = Inf;
        disttol  = 1e-4;
        
        pops(1) = 1;
        
        year = 1;
        lastyear = Inf;
        
        while (disteps > disttol && year < lastyear)
            
            fprintf('Year %3u\n', year);
            
            
            
            % --- Initialize distributions ---
            
            dist = zeros(nk,nz,nb,T_life,3);
            
            im_flow = [ pops(year) * pgr                          ;
                        pops(year) * imm_age(1) * legal_rate(1)   ;
                        pops(year) * imm_age(1) * illegal_rate(1) ];
            
            for ipop = 1:3
                dist(1,:,1,1,ipop) = reshape(proddist_age(:,1,ipop)*im_flow(ipop), [1,nz,1,1,1]);
            end
            
            
            
            % --- Generate distributions through forward propagation ---
            
            for t = 1:T_life-1
                
                age = t;
                
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
                        
                        DIST_step = (dist_previous(:,:,:,t,ipop) + DIST_inflow) .* repmat(reshape(tr_z(:,jz), [1,nz,1]), [nk,1,nb]) * surv(age);
                        
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
                
            end
            
            
            dist_age = sum(sum(reshape(dist, [], T_life, 3), 1), 3);
            disteps = max(abs(dist_age(2:end)/dist_age(1) - dist_age_previous(2:end)/dist_age_previous(1)));
            dist_age_previous = dist_age;
            
            pops(year+1) = sum(dist(:)); %#ok<AGROW>
            dist_previous = dist;
            year = year + 1;
            
        end
        
        
        
        % --- solve_cohort aggregate generation ---
        
        Kalive  = zeros(1,T_life);
        Kdead   = zeros(1,T_life);
        Lab     = zeros(1,T_life);
        ELab    = zeros(1,T_life);
        
        for ipop = 1:3
            for t = 1:T_life
                for iz = 1:nz
                    for ik = 1:nk
                        for ib = 1:nb
                            
                            age = t;
                            
                            Kalive(t) = Kalive(t) + K  (ik,iz,ib,t)               *dist(ik,iz,ib,t,ipop);
                            Kdead (t) = Kdead (t) + K  (ik,iz,ib,t)*(1-surv(age)) *dist(ik,iz,ib,t,ipop);
                            Lab   (t) = Lab   (t) + LAB(ik,iz,ib,t)               *dist(ik,iz,ib,t,ipop);
                            ELab  (t) = ELab  (t) + LAB(ik,iz,ib,t)*z(iz,age,idem)*dist(ik,iz,ib,t,ipop);
                            
                        end
                    end
                end
            end
        end
        
        KPR  = sum(Kalive + Kdead);
        ELAB = sum(ELab);
        
        save(fullfile(jobdir, sprintf('distvars_%u.mat', idem)), ...
             'dist', 'Kalive', 'Kdead', 'KPR', 'ELab', 'Lab', 'ELAB');
        
    end
    
    
    
    % --- Add cohort aggregates to total aggregates ---
    
    kpr     = 0;
    elab    = 0;
    lab     = 0;
    kdeadpr = 0;
    
    for idem = 1:ndem
        
        load(fullfile(jobdir, sprintf('distvars_%u.mat', idem)), ...
             'KPR', 'ELAB', 'Lab', 'Kdead');
        
        lab     = lab     + sum(Lab);
        kpr     = kpr     + KPR;
        elab    = elab    + ELAB;
        kdeadpr = kdeadpr + sum(Kdead);
        
    end
    
    
    
    % --- Calculate additional dynamic aggregates ---
    
    rhopr = (kpr-DD)/elab;
    
    beqeps = beq - kdeadpr;
    beq    = kdeadpr/pop;
    
    rhopr1   = max(0.5, kpr-DD)/elab;
    rhosseps = abs(rhopr1 - rho);
    rho      = 0.25*rhopr1 + 0.75*rho;
    
    Y  = A*((kpr-DD)^alp)*(elab^(1-alp));
    wg = A.*(1-alp).*(rhopr.^alp);
    DEBTss = 0.7*Y;
    
    K_Y = (kpr-DD)  /Y;
    I_Y = (kpr-DD)*d/Y;
    
end



% -- Save optimal decision values and distributions for baseline ---

DIST = zeros(nk,nz,nb,T_life,3,ndem);

for idem = 1:ndem
    load(fullfile(jobdir, sprintf('distvars_%u.mat', idem)), 'dist');
    DIST(:,:,:,:,:,idem) = dist;
end

save(fullfile(jobdir, 'DIST.mat'), 'DIST');




%% Testing

distvars_1        = load(fullfile(jobdir  , 'distvars_1.mat'));
distvars_1_freeze = load(fullfile('Freeze', 'distvars_1.mat'));

fprintf('distvars_1\n');
valuenames = fields(distvars_1);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = distvars_1.(valuename)(:) - distvars_1_freeze.(valuename)(:);
    if any(delta)
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