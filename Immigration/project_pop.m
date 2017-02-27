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



save(fullfile(jobdir, 'imm_polparams.mat'), ...
     'legal_rate', 'prem_legal', 'proddist_age', 'amnesty', 'deportation')




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