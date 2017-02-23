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
        
        load(fullfile(jobdir, sprintf('sspol%u.mat', idem)));
        
        dist_previous     = zeros(nk,nz,nb,T_life,3);
        dist_age_previous = ones(1,T_life);
        
        disteps  = Inf;
        disttol  = 1e-4;
        
        pops(1) = 1;
        
        year = 1;
        lastyear = Inf;
        
        while (disteps > disttol && year < lastyear)
            
            fprintf('Year %2u\n', year);
            
            
            
            % --- Initialize distributions ---
            
            dist = zeros(nk,nz,nb,T_life,3);
            
            % pgr is population growth rate of existing population.  only grows youngest cohort.
            % using period 1 imm rate values for steady state
            im_flow = [ pops(year) * pgr                          ;
                        pops(year) * imm_age(1) * legal_rate(1)   ;
                        pops(year) * imm_age(1) * illegal_rate(1) ];
            
            for iz = 1:nz
                for ipop = 1:3
                    dist(1,iz,1,1,ipop) = proddist_age(iz,1,ipop) * im_flow(ipop);
                end
            end
            
            
            
            % --- Generate distributions through forward propagation ---
            
            for t = 1:T_life-1
                
                age = t;
                
                im_flow = [ 0                                           ;
                            pops(year) * imm_age(age) * legal_rate(1)   ;
                            pops(year) * imm_age(age) * illegal_rate(1) ];
                
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
                                    
                                    dist_hold = dist_previous(ik,iz,ib,t,ipop) + (ik == 1)*(ib == 1)*proddist_age(iz,age,ipop)*im_flow(ipop);
                                    
                                    dist(loc1  ,jz,loc2  ,t+1,ipop) = dist(loc1  ,jz,loc2  ,t+1,ipop) + surv(age)*(1-w2)*(1-w1)*tr_z(iz,jz)*dist_hold;
                                    dist(loc1+1,jz,loc2  ,t+1,ipop) = dist(loc1+1,jz,loc2  ,t+1,ipop) + surv(age)*(1-w2)*(w1  )*tr_z(iz,jz)*dist_hold;
                                    dist(loc1  ,jz,loc2+1,t+1,ipop) = dist(loc1  ,jz,loc2+1,t+1,ipop) + surv(age)*(w2  )*(1-w1)*tr_z(iz,jz)*dist_hold;
                                    dist(loc1+1,jz,loc2+1,t+1,ipop) = dist(loc1+1,jz,loc2+1,t+1,ipop) + surv(age)*(w2  )*(w1  )*tr_z(iz,jz)*dist_hold;
                                    
                                end
                            end
                            
                        end
                    end
                end
                
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
                            
                            Kalive(t) = Kalive(t) + kopt  (ik,iz,ib,t)               *dist(ik,iz,ib,t,ipop);
                            Kdead (t) = Kdead (t) + kopt  (ik,iz,ib,t)*(1-surv(age)) *dist(ik,iz,ib,t,ipop);
                            Lab   (t) = Lab   (t) + labopt(ik,iz,ib,t)               *dist(ik,iz,ib,t,ipop);
                            ELab  (t) = ELab  (t) + labopt(ik,iz,ib,t)*z(iz,age,idem)*dist(ik,iz,ib,t,ipop);
                            
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