function pol_params_imm(impolno)
% impolno = 125;
% 
% iter = 0;
% policy1 = zeros(5,5^3);
% for run1 = 1:5
%     for run2 = 1:5
%         for run3 = 1:5
%                     iter = iter + 1;                    
%                     policy1(:,iter) = [run1, run2, run3,iter,iter + 0];
%         end
%     end
% end
% policy1 = policy1';
% policy4 = policy1(:,4:5);
% 
% save polvals_imm.mat
load polvals_imm.mat
load params.mat Tss T z proddist_age
load Imm_Data.mat

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

for i1 = 2 % legal =2
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

filename = ['imm_polparams_' num2str(impolno) '.mat'];
save(filename)

% Need to solve steady-state before using function below
pop_trans = population_projector(impolno);


filename = ['imm_polparams_' num2str(impolno) '.mat'];
save(filename)


