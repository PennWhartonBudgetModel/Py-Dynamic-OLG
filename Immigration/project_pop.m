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

s = load(fullfile(jobdir, sprintf('sspol%u.mat'    , idem)));
K = s.kopt;
B = s.bopt;
load(fullfile(jobdir, sprintf('distvars_%u.mat', idem)), 'dist');

dist_previous     = dist;   %#ok<NODEF>
dist_age_previous = ones(1,T_life);

disteps  = Inf;
disttol  = -Inf;

pops(1) = sum(dist(:));

year = 1;
lastyear = T_model;


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