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
    for demtype = 1:ndem
        filename = sprintf('sspol%u.mat', demtype);
        copyfile(fullfile('Freeze', filename), fullfile(jobdir, filename));
    end
    
    
    for demtype = 1:ndem
        
        load(fullfile(jobdir, sprintf('sspol%u.mat', demtype)));
        
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
            im_flow = [ pop*(pgr)                                     ;
                        im_scale * pop * imm_age(1) * legal_rate(1)   ;
                        im_scale * pop * imm_age(1) * illegal_rate(1) ];   % using period 1 imm rate values for steady state
            
            loc1 = 1;
            loc2 = 1;
            for i1=1:nz
                dist1ss(loc1,i1,loc2,1,:) = squeeze(proddist_age(i1,1,:)).*(im_flow);
            end
            
            
            for t1 = 1:Tr
                
                im_flow = [ 0                                              ;
                            im_scale * pop * imm_age(t1) * legal_rate(1)   ;
                            im_scale * pop * imm_age(t1) * illegal_rate(1) ];   % using period 1 imm rate values for steady state
                
                for j2 = 1:nz
                    for i1 = 1:nk
                        for i2 = 1:nb
                            
                            point_k = max(kopt(i1,j2,i2,t1),kgrid(1));    % placing floor at lowest gridpoint.
                            loc1 = find(kgrid(1:nk-1)<=point_k,1,'last');
                            w1 = (point_k - kgrid(loc1))/(kgrid(loc1+1)-kgrid(loc1));  % amount allocated to higher gridpoint
                            w1 = min(w1,1);

                            point_b = max(bopt(i1,j2,i2,t1),bgrid(1));    % placing floor at lowest gridpoint.
                            loc2 = find(bgrid(1:nb-1)<=point_b,1,'last');    % lower gridpoint
                            w2 = (point_b - bgrid(loc2))/(bgrid(loc2+1)-bgrid(loc2));  % amount allocated to higher gridpoint
                            w2 = min(w2,1);

                            dist_hold = squeeze(dist1ss_previous(i1,j2,i2,t1,:)) + (eq(i1,1)*eq(i2,1))*squeeze(proddist_age(j2,t1,:)).*im_flow;

                            for j4 = 1:nz
                                if (t1 < Tr)
                                    dist1ss  (loc1  , j4, loc2  , t1+1, :) = squeeze(dist1ss  (loc1  , j4, loc2  , t1+1, :)) + (surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                    dist1ss  (loc1+1, j4, loc2  , t1+1, :) = squeeze(dist1ss  (loc1+1, j4, loc2  , t1+1, :)) + (surv(t1)*(1-w2)*(w1  )*tr_z(j2,j4)).*dist_hold;
                                    dist1ss  (loc1  , j4, loc2+1, t1+1, :) = squeeze(dist1ss  (loc1  , j4, loc2+1, t1+1, :)) + (surv(t1)*(w2  )*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                    dist1ss  (loc1+1, j4, loc2+1, t1+1, :) = squeeze(dist1ss  (loc1+1, j4, loc2+1, t1+1, :)) + (surv(t1)*(w2  )*(w1  )*tr_z(j2,j4)).*dist_hold;
                                else
                                    dist1ss_r(loc1  ,     loc2  , t1+1, :) = squeeze(dist1ss_r(loc1  ,     loc2  , t1+1, :)) + (surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                    dist1ss_r(loc1+1,     loc2  , t1+1, :) = squeeze(dist1ss_r(loc1+1,     loc2  , t1+1, :)) + (surv(t1)*(1-w2)*(w1  )*tr_z(j2,j4)).*dist_hold;
                                    dist1ss_r(loc1  ,     loc2+1, t1+1, :) = squeeze(dist1ss_r(loc1  ,     loc2+1, t1+1, :)) + (surv(t1)*(w2  )*(1-w1)*tr_z(j2,j4)).*dist_hold;
                                    dist1ss_r(loc1+1,     loc2+1, t1+1, :) = squeeze(dist1ss_r(loc1+1,     loc2+1, t1+1, :)) + (surv(t1)*(w2  )*(w1  )*tr_z(j2,j4)).*dist_hold;
                                end
                            end
                            
                        end
                    end
                end
                
            end
            
            
            for t1 = Tr+1:T-1
                
                im_flow = [ 0                                              ;
                            im_scale * pop * imm_age(t1) * legal_rate(1)   ;
                            im_scale * pop * imm_age(t1) * illegal_rate(1) ];   % using period 1 imm rate values for steady state
                
                for i1 = 1:nk
                    for i2 = 1:nb
                        
                        point_k = max(koptss(i1,i2,t1-Tr),kgrid(1));    % placing floor at lowest gridpoint.
                        loc1 = find(kgrid(1:nk-1)<=point_k,1,'last');
                        w1 = (point_k - kgrid(loc1))/(kgrid(loc1+1)-kgrid(loc1));  % amount allocated to higher gridpoint
                        w1 = min(w1,1);
                        
                        dist_hold = squeeze(dist1ss_r_previous(i1,i2,t1,:)) + (eq(i1,1)*eq(i2,1)).*im_flow;
                        
                        dist1ss_r(loc1  , i2, t1+1, :) = squeeze(dist1ss_r(loc1  , i2, t1+1, :)) + (surv(t1).*(1-w1)).*dist_hold;
                        dist1ss_r(loc1+1, i2, t1+1, :) = squeeze(dist1ss_r(loc1+1, i2, t1+1, :)) + (surv(t1).*(w1  )).*dist_hold;
                        
                    end
                end
                
            end
            
            pop_prev = pop;
            pop = sum(dist1ss(:)) + sum(dist1ss_r(:));
            
            dist_age = sum(sum(reshape(dist1ss, [], T, 3), 1) + sum(reshape(dist1ss_r, [], T, 3), 1), 3)';
            
            disteps = max(abs((1/dist_age(1)).*dist_age(2:end) - (1/dist_age_previous(1)).*dist_age_previous(2:end)));
            fprintf('Iteration %3u: %8.4f\n', iter, disteps)
            
            dist_age_previous  = dist_age ;
            dist1ss_previous   = dist1ss  ;   % updating working-age distribution array
            dist1ss_r_previous = dist1ss_r;   % updating retiree distribution array
            
        end
        
        % solving for aggregates by age
        for i3 = 1:3
            for k1 = 1:Tr
                for k2 = 1:nz
                    for k4 = 1:nk
                        for k5 = 1:nb
                            Kalive(k1) = Kalive(k1) + kopt(k4,k2,k5,k1)                   *dist1ss(k4,k2,k5,k1,i3);
                            Kdead (k1) = Kdead (k1) + kopt(k4,k2,k5,k1)  *(1-surv(k1))    *dist1ss(k4,k2,k5,k1,i3);
                            Lab   (k1) = Lab   (k1) + labopt(k4,k2,k5,k1)                 *dist1ss(k4,k2,k5,k1,i3);
                            ELab  (k1) = ELab  (k1) + labopt(k4,k2,k5,k1)*z(k2,k1,demtype)*dist1ss(k4,k2,k5,k1,i3);
                        end
                    end
                end
            end
            for k1 = Tr+1:T
                for k4 = 1:nk
                    for k5 = 1:nb
                        Kalive(k1) = Kalive(k1) + koptss(k4,k5,k1-Tr)             *dist1ss_r(k4,k5,k1,i3);
                        Kdead (k1) = Kdead (k1) + koptss(k4,k5,k1-Tr)*(1-surv(k1))*dist1ss_r(k4,k5,k1,i3);
                    end
                end
            end
        end
        
        KPR  = sum(Kalive + Kdead);
        ELAB = sum(ELab);
        
        save(fullfile(jobdir, sprintf('distvars_%u.mat', demtype)), ...
             'dist1ss', 'dist1ss_r', 'Kalive', 'Kdead', 'KPR', 'ELab', 'Lab', 'ELAB');
        
    end
    
    kpr     = 0;
    elab    = 0;
    lab     = 0;
    kdeadpr = 0;
    
    for demtype = 1:ndem
        
        load(fullfile(jobdir, sprintf('distvars_%u.mat', demtype)), ...
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
    rhosseps = abs(rhopr1-rho);
    rho      = 0.25*rhopr1 + 0.75*rho;
    
    Y  = A*((kpr-DD)^alp)*(elab^(1-alp));
    wg = A.*(1-alp).*(rhopr.^alp);
    DEBTss = 0.7*Y;
    
    K_Y = (kpr-DD)  /Y;
    I_Y = (kpr-DD)*d/Y;
    
end



dist1 = zeros(nk,nz,nb,T,3,ndem);
distr = zeros(nk,   nb,T,3,ndem);

for demtype = 1:ndem
    load(fullfile(jobdir, sprintf('distvars_%u.mat', demtype)), 'dist1ss', 'dist1ss_r');
    dist1(:,:,:,:,:,demtype) = dist1ss  ;
    distr(:,  :,:,:,demtype) = dist1ss_r;
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



