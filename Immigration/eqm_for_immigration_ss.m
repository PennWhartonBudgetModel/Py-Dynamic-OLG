clear; clc;

jobdir = 'Steady State Values';
if exist(jobdir, 'dir'), rmdir(jobdir, 's'), end
mkdir(jobdir);

load('params.mat')
surv(T) = 0;
Tr = NRA_baseline;

load('Imm_Data.mat')


rhosseps = Inf;
rhosstol = 1e-4;

im_scale = 1;    % scaling the immigration flow to target

load('SSVALS.mat')   % loads the last set of outputs as the first guess

% starting steady state equilbrium calculation
while (rhosseps > rhosstol)
    
    DD = DEBTss;  % debt used in the calculation of prices (zero for open economy)
    
    % Copy precomputed dynamic optimization values
    for idem = 1:ndem
        filename = sprintf('sspol%u.mat', idem);
        copyfile(fullfile('Freeze', filename), fullfile(jobdir, filename));
    end
    
    
    for idem = 1:ndem
        
        load(fullfile(jobdir, sprintf('sspol%u.mat', idem)));
        
        Kalive  = zeros(1,T);
        Kdead   = zeros(1,T);
        Lab     = zeros(1,T);
        ELab    = zeros(1,T);
        
        dist1ss_previous   = zeros(nk,nz,nb,T,3);
        dist1ss_r_previous = zeros(nk,   nb,T,3);
        
        disteps  = Inf;
        disttol  = 1e-4;
        pop_prev = 0;
        pop      = 1;   % initial assertion of population
        
        dist_age_previous = ones(T,1);
        
        iter = 0;
        
        while (disteps > disttol)
            
            iter = iter + 1;
            
            dist1ss   = zeros(nk,nz,nb,T,3);   % last dimension is [native, legal, illegal]
            dist1ss_r = zeros(nk,   nb,T,3);
            
            % pgr is population growth rate of existing population.  only grows youngest cohort.
            % using period 1 imm rate values for steady state
            im_flow = [ pop * pgr                                     ;
                        im_scale * pop * imm_age(1) * legal_rate(1)   ;
                        im_scale * pop * imm_age(1) * illegal_rate(1) ];
            im_flow = reshape(im_flow, [1,1,1,1,3]);
            
            for iz = 1:nz
                dist1ss(1,iz,1,1,:) = reshape(proddist_age(iz,1,:), [1,1,1,1,3]) .* im_flow;
            end
            
            
            for t = 1:Tr
                
                % using period 1 imm rate values for steady state
                im_flow = [ 0                                             ;
                            im_scale * pop * imm_age(t) * legal_rate(1)   ;
                            im_scale * pop * imm_age(t) * illegal_rate(1) ];
                im_flow = reshape(im_flow, [1,1,1,1,3]);
                
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
                            
                            dist_hold = dist1ss_previous(ik,iz,ib,t,:) + (ik == 1)*(ib == 1)*reshape(proddist_age(iz,t,:), [1,1,1,1,3]).*im_flow;
                            
                            for jz = 1:nz
                                if (t < Tr)
                                    dist1ss  (loc1  , jz, loc2  , t+1, :) = dist1ss  (loc1  , jz, loc2  , t+1, :) + surv(t)*(1-w2)*(1-w1)*tr_z(iz,jz)*dist_hold;
                                    dist1ss  (loc1+1, jz, loc2  , t+1, :) = dist1ss  (loc1+1, jz, loc2  , t+1, :) + surv(t)*(1-w2)*(w1  )*tr_z(iz,jz)*dist_hold;
                                    dist1ss  (loc1  , jz, loc2+1, t+1, :) = dist1ss  (loc1  , jz, loc2+1, t+1, :) + surv(t)*(w2  )*(1-w1)*tr_z(iz,jz)*dist_hold;
                                    dist1ss  (loc1+1, jz, loc2+1, t+1, :) = dist1ss  (loc1+1, jz, loc2+1, t+1, :) + surv(t)*(w2  )*(w1  )*tr_z(iz,jz)*dist_hold;
                                else
                                    dist1ss_r(loc1  ,     loc2  , t+1, :) = dist1ss_r(loc1  ,     loc2  , t+1, :) + surv(t)*(1-w2)*(1-w1)*tr_z(iz,jz)*reshape(dist_hold, [1,1,1,3]);
                                    dist1ss_r(loc1+1,     loc2  , t+1, :) = dist1ss_r(loc1+1,     loc2  , t+1, :) + surv(t)*(1-w2)*(w1  )*tr_z(iz,jz)*reshape(dist_hold, [1,1,1,3]);
                                    dist1ss_r(loc1  ,     loc2+1, t+1, :) = dist1ss_r(loc1  ,     loc2+1, t+1, :) + surv(t)*(w2  )*(1-w1)*tr_z(iz,jz)*reshape(dist_hold, [1,1,1,3]);
                                    dist1ss_r(loc1+1,     loc2+1, t+1, :) = dist1ss_r(loc1+1,     loc2+1, t+1, :) + surv(t)*(w2  )*(w1  )*tr_z(iz,jz)*reshape(dist_hold, [1,1,1,3]);
                                end
                            end
                            
                        end
                    end
                end
                
            end
            
            
            for t = Tr+1:T-1
                
                % using period 1 imm rate values for steady state
                im_flow = [ 0                                             ;
                            im_scale * pop * imm_age(t) * legal_rate(1)   ;
                            im_scale * pop * imm_age(t) * illegal_rate(1) ];
                im_flow = reshape(im_flow, [1,1,1,3]);
                
                for ik = 1:nk
                    for ib = 1:nb
                        
                        point_k = max(koptss(ik,ib,t-Tr), kgrid(1));
                        loc1 = find(kgrid(1:nk-1) <= point_k, 1, 'last');
                        w1 = (point_k - kgrid(loc1)) / (kgrid(loc1+1) - kgrid(loc1));
                        w1 = min(w1, 1);
                        
                        dist_hold = dist1ss_r_previous(ik,ib,t,:) + (ik == 1)*(ib == 1)*im_flow;
                        
                        dist1ss_r(loc1  , ib, t+1, :) = dist1ss_r(loc1  , ib, t+1, :) + surv(t)*(1-w1)*dist_hold;
                        dist1ss_r(loc1+1, ib, t+1, :) = dist1ss_r(loc1+1, ib, t+1, :) + surv(t)*(w1  )*dist_hold;
                        
                    end
                end
                
            end
            
            pop_prev = pop;
            pop = sum(dist1ss(:)) + sum(dist1ss_r(:));
            
            dist_age = sum(sum([ reshape(dist1ss, [], T, 3) ; reshape(dist1ss_r, [], T, 3) ], 1), 3)';
            
            disteps = max(abs(dist_age(2:end)/dist_age(1) - dist_age_previous(2:end)/dist_age_previous(1)));
            fprintf('Iteration %3u: %8.4f\n', iter, disteps)
            
            dist_age_previous  = dist_age ;
            dist1ss_previous   = dist1ss  ;   % updating working-age distribution array
            dist1ss_r_previous = dist1ss_r;   % updating retiree distribution array
            
        end
        
        % solving for aggregates by age
        for ipop = 1:3
            for t = 1:Tr
                for iz = 1:nz
                    for ik = 1:nk
                        for ib = 1:nb
                            Kalive(t) = Kalive(t) + kopt  (ik,iz,ib,t)             *dist1ss(ik,iz,ib,t,ipop);
                            Kdead (t) = Kdead (t) + kopt  (ik,iz,ib,t)*(1-surv(t)) *dist1ss(ik,iz,ib,t,ipop);
                            Lab   (t) = Lab   (t) + labopt(ik,iz,ib,t)             *dist1ss(ik,iz,ib,t,ipop);
                            ELab  (t) = ELab  (t) + labopt(ik,iz,ib,t)*z(iz,t,idem)*dist1ss(ik,iz,ib,t,ipop);
                        end
                    end
                end
            end
            for t = Tr+1:T
                for ik = 1:nk
                    for ib = 1:nb
                        Kalive(t) = Kalive(t) + koptss(ik,ib,t-Tr)            *dist1ss_r(ik,ib,t,ipop);
                        Kdead (t) = Kdead (t) + koptss(ik,ib,t-Tr)*(1-surv(t))*dist1ss_r(ik,ib,t,ipop);
                    end
                end
            end
        end
        
        KPR  = sum(Kalive + Kdead);
        ELAB = sum(ELab);
        
        save(fullfile(jobdir, sprintf('distvars_%u.mat', idem)), ...
             'dist1ss', 'dist1ss_r', 'Kalive', 'Kdead', 'KPR', 'ELab', 'Lab', 'ELAB');
        
    end
    
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



dist1 = zeros(nk,nz,nb,T,3,ndem);
distr = zeros(nk,   nb,T,3,ndem);

for idem = 1:ndem
    load(fullfile(jobdir, sprintf('distvars_%u.mat', idem)), 'dist1ss', 'dist1ss_r');
    dist1(:,:,:,:,:,idem) = dist1ss  ;
    distr(:,  :,:,:,idem) = dist1ss_r;
end

save(fullfile(jobdir, 'eqmdist.mat'), 'dist1', 'distr');




%% Testing

distvars_1   = load(fullfile(jobdir  , 'distvars_1.mat'));
distvars_1_0 = load(fullfile('Freeze', 'distvars_1.mat'));

fprintf('distvars_1\n');
valuenames = fields(distvars_1);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = distvars_1.(valuename)(:) - distvars_1_0.(valuename)(:);
    if any(delta)
        pdev = abs(nanmean(delta*2 ./ (distvars_1.(valuename)(:) + distvars_1_0.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');


distvars_2   = load(fullfile(jobdir  , 'distvars_2.mat'));
distvars_2_0 = load(fullfile('Freeze', 'distvars_2.mat'));

fprintf('distvars_2\n');
valuenames = fields(distvars_2);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = distvars_2.(valuename)(:) - distvars_2_0.(valuename)(:);
    if any(delta)
        pdev = abs(nanmean(delta*2 ./ (distvars_2.(valuename)(:) + distvars_2_0.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');



