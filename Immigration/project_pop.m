function [] = project_pop()



% --- Initialization ---

jobdir = 'Testing';



% --- Parameter loading ---

load('params.mat')
load('Imm_Data.mat')

T_life  = T;
T_model = Tss;


% policy #1 increases the flow of legal immigrants
legal_rate_scale = 1.5;
legal_rate = legal_rate * legal_rate_scale; %#ok<NODEF>


% policy #2 increases the relative flow of skilled immigrants.
% the idea here is to get the increase in avg productivity from someone and
% use that value to scale up
prem_legal = 1.117456794;


% This loop modifies proddist_age to hit productivity targets specified above.  We do this by moving 
% masses from low -> high or high -> low.  The problem is bounded by a max productivity,
% so not all productivity targets are feasible.

for i1 = 2  % *** idem = 1:2?
    for t = 1:T_life
        % We take the average over productivity shocks because the probability weights are uniform
        % (by 2-bin discretization of normal distributions)
        E_bar = sum(reshape(z(:,t,:), [], 1))/8; %#ok<NODEF>
        pub = 1;
        plb = 0;
        target_e = E_bar * prem_legal;
        e_err = 100;
        e_tol = 1e-5;
        while (e_err > e_tol)
            p_upper = (pub + plb)/2;
            p_lower = (1 - p_upper)/3;
            val1 = p_lower*z(1,t,1) + p_lower*z(2,t,1) + p_lower*z(3,t,1) + p_upper*z(4,t,1);
            val2 = p_lower*z(1,t,2) + p_lower*z(2,t,2) + p_lower*z(3,t,2) + p_upper*z(4,t,2);
            test_e = (val1 + val2)/2;
            if (test_e > target_e)
                pub = p_upper;
            else
                plb = p_upper;
            end
            e_err = abs(test_e - target_e);
        end
        proddist_age(1,t,i1) = p_lower; %#ok<AGROW>
        proddist_age(2,t,i1) = p_lower; %#ok<AGROW>
        proddist_age(3,t,i1) = p_lower; %#ok<AGROW>     % *** p_upper?
        proddist_age(4,t,i1) = p_upper; %#ok<AGROW>
    end
end  


% policy #3 grants amnesty to a percentage of illegal immigrants annually
amnesty     = 0.05;
deportation = 0.05;




idem = 1;  % by the symmetry of the demographic type sizes, population will grow equally on each "island" (i.e., idem island).  WLOG, we choose idem=1.

load(fullfile(jobdir, sprintf('sspol%u.mat'    , idem)));
load(fullfile(jobdir, sprintf('distvars_%u.mat', idem)), 'dist');

dist_previous     = dist;   %#ok<NODEF>
dist_age_previous = ones(1,T_life);
        
disteps  = Inf;
disttol  = -Inf;

pops(1) = sum(dist(:));

year = 1;
lastyear = T_model;

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
    
    
    % *** Amnesty adjustments missing?
    dist(:,:,:,:,3) = (1-deportation)*dist(:,:,:,:,3);
    
    
    dist_age = sum(sum(reshape(dist, [], T_life, 3), 1), 3);
    disteps = max(abs(dist_age(2:end)/dist_age(1) - dist_age_previous(2:end)/dist_age_previous(1)));
    dist_age_previous = dist_age;
    
    pops(year+1) = sum(dist(:)); %#ok<AGROW>
    dist_previous = dist;
    year = year + 1;
    
end


save(fullfile(jobdir, 'imm_polparams.mat'), ...
     'legal_rate', 'prem_legal', 'proddist_age', 'amnesty', 'deportation', 'pops')




%% Testing

imm_polparams        = load(fullfile(jobdir  , 'imm_polparams.mat'));
imm_polparams_freeze = load(fullfile('Freeze', 'imm_polparams.mat'));

fprintf('imm_polparams\n');
valuenames = fields(imm_polparams);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = imm_polparams.(valuename)(:) - imm_polparams_freeze.(valuename)(:);
    if any(isnan(delta))
        fprintf('\t%-14sNaN found\n', valuename);
    elseif any(delta)
        pdev = abs(nanmean(delta*2 ./ (imm_polparams.(valuename)(:) + imm_polparams_freeze.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');



end