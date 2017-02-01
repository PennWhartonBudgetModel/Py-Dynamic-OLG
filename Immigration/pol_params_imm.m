function [] = pol_params_imm()

impolno = 1;
jobdir = 'Testing';

load('polvals_imm.mat')
load('params.mat', 'Tss', 'T', 'z', 'proddist_age')
load('Imm_Data.mat')

policyinds = policy1(impolno, 1:3);


% policy #1 increases the flow of legal immigrants
legal_rate_scales = [1.00, 1.25, 1.50, 1.75, 2.00];
legal_rate = legal_rate * legal_rate_scales(policyinds(1)); %#ok<NODEF>


% policy #2 increases the relative flow of skilled immigrants.
% the idea here is to get the increase in avg productivity from someone and
% use that value to scale up
prem_legals = [1.000000000, 1.058728397, 1.117456794, 1.176185192, 1.234913589];
prem_legal = prem_legals(policyinds(2));


%  This loop modifies proddist_age to hit productivity targets specified above.  We do this by moving 
%  masses from low -> high or high -> low.  The problem is bounded by a max productivity,
%  so not all productivity targets are feasible.

for i1 = 2  % *** idem = 1:2?
    for t = 1:T
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
amnestys     = [0.00, 0.00, 0.00, 0.05, 0.10];
deportations = [0.00, 0.10, 0.05, 0.00, 0.00];
amnesty     = amnestys    (policyinds(3)); %#ok<NASGU>
deportation = deportations(policyinds(3));



load('params.mat')
load('Imm_Data.mat')

Tr = NRA_baseline;

load('SSVALS.mat', 'pop_prev')

idem = 1;  % by the symmetry of the demographic type sizes, population will grow equally on each "island" (i.e., idem island).  WLOG, we choose idem=1.

load(fullfile(jobdir, sprintf('sspol%u.mat'    , idem)));
load(fullfile(jobdir, sprintf('distvars_%u.mat', idem)), 'dist1ss', 'dist1ss_r');

dist1ss_previous   = dist1ss  ;   %#ok<NODEF> % initiating working-age distribution array
dist1ss_r_previous = dist1ss_r;   %#ok<NODEF> % initiating retiree distribution array

pop = sum(dist1ss(:)) + sum(dist1ss_r(:));

pop_trans = zeros(Tss,1);    % fill in with transition population sizes
pop_trans(1) = pop;


for trans_year = 2:Tss
    
    fprintf('Year %2u\n', trans_year);
    
    year = trans_year - 1;
    
    dist1ss = zeros(nk,nz,nb,T,3);   % last dimension is [native, legal, illegal]
    dist1ss_r = zeros(nk,nb,T,3);
    
    % using period 1 imm rate values for steady state
    im_flow = [ pop_trans(year) * pgr                          ;
                pop_trans(year) * imm_age(1) * legal_rate(1)   ;
                pop_trans(year) * imm_age(1) * illegal_rate(1) ];
    
    for iz = 1:nz
        dist1ss(1,iz,1,1,:) = reshape(proddist_age(iz,1,:), 3, []) .* im_flow;
    end
    
    for t = 1:Tr
        
        age = t;
        
        im_flow = [ 0                                                ;
                    pop_trans(year) * imm_age(age) * legal_rate(1)   ;
                    pop_trans(year) * imm_age(age) * illegal_rate(1) ];
        
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
                            
                            if (t < Tr)
                                
                                % lack of generality in im flow indicator below
                                dist_hold = dist1ss_previous(ik,iz,ib,t,ipop) + (ik == 1)*(ib == 1)*proddist_age(iz,t,ipop)*im_flow(ipop);
                                
                                dist1ss(loc1  ,jz,loc2  ,t+1,ipop) = dist1ss(loc1  ,jz,loc2  ,t+1,ipop) + surv(t)*(1-w2)*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist1ss(loc1+1,jz,loc2  ,t+1,ipop) = dist1ss(loc1+1,jz,loc2  ,t+1,ipop) + surv(t)*(1-w2)*(w1  )*tr_z(iz,jz)*dist_hold;
                                dist1ss(loc1  ,jz,loc2+1,t+1,ipop) = dist1ss(loc1  ,jz,loc2+1,t+1,ipop) + surv(t)*(w2  )*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist1ss(loc1+1,jz,loc2+1,t+1,ipop) = dist1ss(loc1+1,jz,loc2+1,t+1,ipop) + surv(t)*(w2  )*(w1  )*tr_z(iz,jz)*dist_hold;
                                
                            else
                                
                                dist_hold = dist1ss_previous(ik,iz,ib,t,ipop) + (ik == 1)*(ib == 1)*proddist_age(iz,t,ipop)*im_flow(ipop);
                                
                                dist1ss_r(loc1  ,loc2  ,t+1,ipop) = dist1ss_r(loc1  ,loc2  ,t+1,ipop) + surv(t)*(1-w2)*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist1ss_r(loc1+1,loc2  ,t+1,ipop) = dist1ss_r(loc1+1,loc2  ,t+1,ipop) + surv(t)*(1-w2)*(w1  )*tr_z(iz,jz)*dist_hold;
                                dist1ss_r(loc1  ,loc2+1,t+1,ipop) = dist1ss_r(loc1  ,loc2+1,t+1,ipop) + surv(t)*(w2  )*(1-w1)*tr_z(iz,jz)*dist_hold;
                                dist1ss_r(loc1+1,loc2+1,t+1,ipop) = dist1ss_r(loc1+1,loc2+1,t+1,ipop) + surv(t)*(w2  )*(w1  )*tr_z(iz,jz)*dist_hold;
                            
                            end
                            
                        end
                    end
                    
                end
            end
        end
    end
    
    for t = Tr+1:T-1
        
        age = t;
        
        im_flow = [ 0                                                ;
                    pop_trans(year) * imm_age(age) * legal_rate(1)   ;
                    pop_trans(year) * imm_age(age) * illegal_rate(1) ];
        
        for ik = 1:nk
            for ib = 1:nb
                
                point_k = max(koptss(ik,ib,t-Tr), kgrid(1));
                loc1 = find(kgrid(1:nk-1) <= point_k, 1, 'last');
                w1 = (point_k - kgrid(loc1)) / (kgrid(loc1+1) - kgrid(loc1));
                w1 = min(w1, 1);
                
                for ipop = 1:3
                    
                    dist_hold = dist1ss_r_previous(ik,ib,t,ipop) + (ik == 1)*(ib == 1)*im_flow(ipop);
                    
                    dist1ss_r(loc1  ,ib,t+1,ipop) = dist1ss_r(loc1  ,ib,t+1,ipop) + surv(t)*(1-w1)*dist_hold;
                    dist1ss_r(loc1+1,ib,t+1,ipop) = dist1ss_r(loc1+1,ib,t+1,ipop) + surv(t)*(w1  )*dist_hold;
                
                end
                
            end
        end
    end
    
    dist1ss  (:,:,:,:,3) = (1 - deportation)*dist1ss  (:,:,:,:,3);  
    dist1ss_r(:,  :,:,3) = (1 - deportation)*dist1ss_r(:,  :,:,3);
    
    pop = sum(dist1ss(:)) + sum(dist1ss_r(:));
    pop_trans(trans_year) = pop;
    
    dist1ss_previous   = dist1ss  ;   % updating working-age distribution array
    dist1ss_r_previous = dist1ss_r;   % updating retiree distribution array
    
end


save(fullfile(jobdir, sprintf('imm_polparams_%u.mat', impolno)), ...
     'legal_rate', 'prem_legal', 'proddist_age', 'amnesty', 'deportation', 'pop_trans')




%% Testing

imm_polparams_1   = load(fullfile(jobdir  , 'imm_polparams_1.mat'));
imm_polparams_1_0 = load(fullfile('Freeze', 'imm_polparams_1.mat'));

fprintf('imm_polparams_1\n');
valuenames = fields(imm_polparams_1);
for i = 1:length(valuenames)
    valuename = valuenames{i};
    delta = imm_polparams_1.(valuename)(:) - imm_polparams_1_0.(valuename)(:);
    if any(delta)
        pdev = abs(nanmean(delta*2 ./ (imm_polparams_1.(valuename)(:) + imm_polparams_1_0.(valuename)(:))))*100;
        fprintf('\t%-14s%06.2f%% deviation\n', valuename, pdev);
    else
        fprintf('\t%-14sNo deviation\n', valuename);
    end
end
fprintf('\n');



end