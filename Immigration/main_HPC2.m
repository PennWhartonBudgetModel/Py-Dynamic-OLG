
N = getenv('SGE_TASK_ID');
% run eqm_for_ss


% [~] = pol_params_imm(impolno);

iteration = str2num(N);

% load('initialrun.mat', 'runvals')

% runvals goes from 2 to 256
% polno = runvals(iteration);

% runs = [2,3,4,5,9,13,1];

% for i1 = 1:length(runs)


%     polno =runs(polno);

%     [~] = policy_params3(polno);
    
    polno=1;
    impolno = iteration;
    pol_params_imm(impolno)
    % impolno = str2num(N);

    jobdir = ['job_' num2str(polno) '_' num2str(impolno)];  % gives name
    if exist(jobdir, 'dir')==7          % 7 = directory.
        rmdir(jobdir, 's')                     % remove all files and that directory
    end
    mkdir(jobdir);

    run eqm_for_trans_par
    if impolno>1
        rmdir(jobdir, 's')
    end
% end
% [~] = eqm_for_trans(polno,impolno);

