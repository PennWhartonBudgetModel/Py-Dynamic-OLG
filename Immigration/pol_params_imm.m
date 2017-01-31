function [] = pol_params_imm()

impolno = 1;
jobdir = 'Testing';

load('polvals_imm.mat')
load('params.mat', 'Tss', 'T', 'z', 'proddist_age')
load('Imm_Data.mat')

policy1 = policy1(impolno,:);


% policy #1 increases the flow of legal immigrants
if policy1(1)==1
    legal_rate = legal_rate.*1;
elseif policy1(1)==2
    legal_rate = legal_rate.*1.25;
elseif policy1(1)==3
    legal_rate = legal_rate.*1.5;   
elseif policy1(1)==4
    legal_rate = legal_rate.*1.75;
elseif policy1(1)==5
    legal_rate = legal_rate.*2;
end


% policy #2 increases the relative flow of skilled immigrants.
% the idea here is to get the increase in avg productivity from someone and
% use that value to scale up

policy_premium = [1, 1.058728397, 1.117456794, 1.176185192, 1.234913589];

% dist_year_leg = zeros(Tss,np);
% dist_year_leg(1,:) = demdist_2015;
if policy1(2)==1
    prem_legal = policy_premium(1);
elseif policy1(2)==2
    prem_legal = policy_premium(2);
elseif policy1(2)==3
    prem_legal = policy_premium(3);
elseif policy1(2)==4
    prem_legal = policy_premium(4);
elseif policy1(2)==5
    prem_legal = policy_premium(5);
end


%  This loop modifies proddist_age to hit productivity targets specified above.  We do this by moving 
%  masses from low -> high or high -> low.  The problem is bounded by a max productivity,
%  so not all productivity targets are feasible.

for i1 = 2
    for t1 = 1:T
        % We take the average over productivity shocks because the probability weights are uniform
        % (by 2-bin discretization of normal distributions)
        E_bar = sum(sum(squeeze(z(:,t1,:))))/8;
        pub = 1;
        plb = 0;
        target_e = E_bar*(prem_legal);
        e_error = 100;
        e_tol = .00001;
        while e_error>e_tol
            p_upper = (pub+plb)/2;
            p_lower = (1-p_upper)/3;
            val1 = p_lower*z(1,t1,1) + p_lower*z(2,t1,1) + p_lower*z(3,t1,1) + p_upper*z(4,t1,1);
            val2 = p_lower*z(1,t1,2) + p_lower*z(2,t1,2) + p_lower*z(3,t1,2) + p_upper*z(4,t1,2);
            test_e = .5*(val1+val2);
            if test_e>target_e
                pub = p_upper;
            elseif test_e<target_e
                plb = p_upper;
            end
            e_error = abs(test_e-target_e);
        end
        proddist_age(1,t1,i1) = p_lower;
        proddist_age(2,t1,i1) = p_lower;
        proddist_age(3,t1,i1) = p_lower;
        proddist_age(4,t1,i1) = p_upper;
    end
end  


% policy #3 grants amnesty to a percentage of illegal immigrants annually

if policy1(3) ==1
    amnesty = 0;
    deportation = 0; 
elseif policy1(3)==2
    amnesty = 0;
    deportation = .1;
elseif policy1(3)==3
    amnesty = 0;
    deportation = .05;
elseif policy1(3)==4
    amnesty = .05;
    deportation = 0;
elseif policy1(3)==5
    amnesty = .1;
    deportation = 0;
end


load('params.mat')
load('Imm_Data.mat')

Tr = NRA_baseline;

load('SSVALS.mat', 'pop_prev')

demtype = 1;  % by the symmetry of the demographic type sizes, population will grow equally on each "island" (i.e., demtype island).  WLOG, we choose demtype=1.
filename = ['sspol' num2str(demtype) '.mat'];
totfile = fullfile(jobdir,filename);
load(totfile);

filename = ['distvars_' num2str(demtype) '.mat'];
totfile = fullfile(jobdir,filename);
load(totfile, 'dist1ss', 'dist1ss_r');

dist1ss_previous   = dist1ss  ;   % initiating working-age distribution array
dist1ss_r_previous = dist1ss_r;   % initiating retiree distribution array

pop = sum(dist1ss(:)) + sum(dist1ss_r(:));

pop_trans = zeros(Tss,1);    % fill in with transition population sizes
pop_trans(1) = pop;


for trans_year = 2:Tss
    
    fprintf('%u\n', trans_year);
    
    year = trans_year - 1;
    
    dist1ss = zeros(nk,nz,nb,T,3);   % last dimension is [native, legal, illegal]
    dist1ss_r = zeros(nk,nb,T,3);
    
    % using period 1 imm rate values for steady state
    im_flow = [ pop_trans(year) * pgr                          ;
                pop_trans(year) * imm_age(1) * legal_rate(1)   ;
                pop_trans(year) * imm_age(1) * illegal_rate(1) ];
    
    for i1 = 1:nz
        dist1ss(1,i1,1,1,:) = squeeze(proddist_age(i1,1,:)).*(im_flow);
    end
    
    for t1 = 1:Tr
        
        age = t1;
        
        im_flow = [ 0                                                ;
                    pop_trans(year) * imm_age(age) * legal_rate(1)   ;
                    pop_trans(year) * imm_age(age) * illegal_rate(1) ];
        
        for j2 = 1:nz
            for i1 = 1:nk
                for i2 = 1:nb
                    
                    point_k = max(kopt(i1,j2,i2,t1), kgrid(1));
                    loc1 = find(kgrid(1:nk-1) <= point_k, 1, 'last');
                    w1 = (point_k - kgrid(loc1)) / (kgrid(loc1+1) - kgrid(loc1));
                    w1 = min(w1, 1);
                    
                    point_b = max(bopt(i1,j2,i2,t1), bgrid(1));
                    loc2 = find(bgrid(1:nb-1) <= point_b, 1, 'last');
                    w2 = (point_b - bgrid(loc2)) / (bgrid(loc2+1) - bgrid(loc2));
                    w2 = min(w2, 1);
                    
                    for j4 = 1:nz
                        for ipop = 1:3
                            
                            if (t1 < Tr)
                                
                                % lack of generality in im flow indicator below
                                dist_hold = dist1ss_previous(i1,j2,i2,t1,ipop) + (i1 == 1)*(i2 == 1)*proddist_age(j2,t1,ipop)*im_flow(ipop);
                                
                                dist1ss(loc1  ,j4,loc2  ,t1+1,ipop) = dist1ss(loc1  ,j4,loc2  ,t1+1,ipop) + surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)*dist_hold;
                                dist1ss(loc1+1,j4,loc2  ,t1+1,ipop) = dist1ss(loc1+1,j4,loc2  ,t1+1,ipop) + surv(t1)*(1-w2)*(w1  )*tr_z(j2,j4)*dist_hold;
                                dist1ss(loc1  ,j4,loc2+1,t1+1,ipop) = dist1ss(loc1  ,j4,loc2+1,t1+1,ipop) + surv(t1)*(w2  )*(1-w1)*tr_z(j2,j4)*dist_hold;
                                dist1ss(loc1+1,j4,loc2+1,t1+1,ipop) = dist1ss(loc1+1,j4,loc2+1,t1+1,ipop) + surv(t1)*(w2  )*(w1  )*tr_z(j2,j4)*dist_hold;
                                
                            else
                                
                                dist_hold = dist1ss_previous(i1,j2,i2,t1,ipop) + (i1 == 1)*(i2 == 1)*proddist_age(j2,t1,ipop)*im_flow(ipop);
                                
                                dist1ss_r(loc1  ,loc2  ,t1+1,ipop) = dist1ss_r(loc1  ,loc2  ,t1+1,ipop) + surv(t1)*(1-w2)*(1-w1)*tr_z(j2,j4)*dist_hold;
                                dist1ss_r(loc1+1,loc2  ,t1+1,ipop) = dist1ss_r(loc1+1,loc2  ,t1+1,ipop) + surv(t1)*(1-w2)*(w1  )*tr_z(j2,j4)*dist_hold;
                                dist1ss_r(loc1  ,loc2+1,t1+1,ipop) = dist1ss_r(loc1  ,loc2+1,t1+1,ipop) + surv(t1)*(w2  )*(1-w1)*tr_z(j2,j4)*dist_hold;
                                dist1ss_r(loc1+1,loc2+1,t1+1,ipop) = dist1ss_r(loc1+1,loc2+1,t1+1,ipop) + surv(t1)*(w2  )*(w1  )*tr_z(j2,j4)*dist_hold;
                            
                            end
                            
                        end
                    end
                    
                end
            end
        end
    end
    
    for t1 = Tr+1:T-1
        
        age = t1;
        
        im_flow = [0                                                ;
                   pop_trans(year) * imm_age(age) * legal_rate(1)   ;
                   pop_trans(year) * imm_age(age) * illegal_rate(1) ];
        
        for i1 = 1:nk
            for i2 = 1:nb
                
                point_k = max(koptss(i1,i2,t1-Tr), kgrid(1));
                loc1 = find(kgrid(1:nk-1) <= point_k, 1, 'last');
                w1 = (point_k - kgrid(loc1)) / (kgrid(loc1+1) - kgrid(loc1));
                w1 = min(w1, 1);
                
                for ipop = 1:3
                    
                    dist_hold = dist1ss_r_previous(i1,i2,t1,ipop) + (i1 == 1)*(i2 == 1)*im_flow(ipop);
                    
                    dist1ss_r(loc1  ,i2,t1+1,ipop) = dist1ss_r(loc1  ,i2,t1+1,ipop) + surv(t1)*(1-w1)*dist_hold;
                    dist1ss_r(loc1+1,i2,t1+1,ipop) = dist1ss_r(loc1+1,i2,t1+1,ipop) + surv(t1)*(w1  )*dist_hold;
                
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


filename = ['imm_polparams_' num2str(impolno) '.mat'];
save(fullfile(jobdir, filename), 'pop_trans')




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