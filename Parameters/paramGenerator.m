classdef paramGenerator

methods (Static)
    
    function [param_global] = generate_global_parameters()
        
        % This method generates the global parameters that are determined internally to the code and have no outside dependencies.
        
        %% Miscelaneous deep parameters of the model
        
        mpci = 3.25;   % Model per capita income -- obsolete.  Should be updated in every different baseline.
        rpci = 119000; % Real per capita income in dollars.  Chosen to match the percentage of individuals above the "taxmax," defined below in the Social Security section. (6%)
        
        A    = 1;    % Total factor productivity.
        alp  = 0.45;    % Capital share of output.  "alp" is short for "alpha," the Greek letter commonly reserved for this parameter.
        d = .085;   % Depreciation rate.
        T_life = 80;    % Maximum possible lifetime.
        
        T_model = 24;    % Length of transition path.
        pgr = 0.02;    % Population annual growth rate
        
        surv = calculate_survival();    % Calculating survival probabilities using a function.

        
        % Normalized stationary distribution over age, accounting for the population growth rate
        Mu2 = zeros(1,T_life);
        Mu2(1) = 1;
        for i1 = 2:T_life
            Mu2(i1) = (surv(i1-1)/(1+pgr))*Mu2(i1-1);
        end
        
        
        %% Defining the shock process. (Storesletten, Telmer, Yaron (2004))
        
        % Permanent shock
        
        np = 2;  % Number of permanent shocks
        sig_p = (.2105)^(.5);    % Standard deviation
        [z_perm, ~] = generate_normdist_shocks(np,sig_p);
        demdist_2015 = (1/np).*ones(np,1);    % Determines the initial distribution of permanent types.  After trying alternative specifications, we decided to make this stationary.
        
        
        % Persistent shock
        
        nz = 2;
        sigep_1 = (.124)^.5;
        const = 0;
        pep = .973;    % coefficient on lagged productivity 
        sigep = (.018)^(.5);    % conditional stdev of productivity 
        [z1, tr_z1, epsgrid] = generate_persistent_shocks(nz,const,pep,sigep);    % returns the transition matrix
        proddist = zeros(1,nz);
        for i1 = 1:nz
            proddist(i1) = normcdf(epsgrid(i1+1),0,sigep_1) - normcdf(epsgrid(i1),0,sigep_1);
        end
        
        
        % Transitory shock 

        nt = 2;  % number of transitory shocks
        sig_t = (.0630)^(.5);
        [z_trans, ~] = generate_normdist_shocks(nt,sig_t); 
        tr_t = (1/nt).*ones(nt);  % Initial distribution over transitory shocks.
        
        tr_z = kron(tr_z1,tr_t);    % Provides final productivity Markov transition matrix.
        
        proddist = kron(proddist,(1/nt).*ones(1,nt));    % Initial distribution over all shocks.
        
        z = zeros(nz*nt, T_life, np);
        % Deterministic component of lifecycle productivity.  Estimates are from Barro and Barnes Medicaid working paper.
        for t1 = 1:T_life
            for d1 = 1:np
                shock_node = 0;
                for idio_shock = 1:nz
                    for trans_shock = 1:nt
                        shock_node = shock_node+1;
                        z(shock_node,t1,d1) = (t1+19)* 0.2891379 -0.0081467*(t1+19)^2 +  0.000105*(t1+19)^3 -5.25E-07*(t1+19)^4 -1.203521 + z1(nz - (idio_shock-1)) + z_trans(trans_shock) + z_perm(d1);
                    end
                end
            end
        end

         z = max(0,z);
         nz = nz*nt;  % reassigning nz so that we don't have to make any other adjustments in the code.
         
         ndem = np;    % Redefining "np" as "ndem" to reflect patch for the creation of demographic types.

         
         
         
        %% Creating the capital and average earnings grids using non-linear spacing
        
        % Capital grid
        nk = 10;
        kgrid = zeros(nk-1,1);
        kub = 150;
        klb = .001;
        exppwr = 2;
        for i1 = 1:nk-1
            kgrid(i1) = (kub-klb)*((i1/nk))^exppwr + klb;
        end
        kgrid = [klb;kgrid];
        
        
        % Average earnings grid (Note the dependence of the average earnings grid on the maximum of the productivity array).
        nb = 5;
        bub = 1.5*max(max(max(z)));  % max possible ss benefit
        bgrid = zeros(nb-1,1);
        blb = 0.1;
        exppwr = 3;
        for i1 = 1:nb-1
            bgrid(i1) = (bub-blb)*(((i1+1)/nb))^exppwr + blb;
        end
        bgrid = [blb;bgrid];
        
        % Generating structure to contain all global parameters.
        param_global.mpci = mpci;
        param_global.rpci = rpci;
        param_global.A = A;
        param_global.alp = alp;
        param_global.d = d;
        param_global.T_life = T_life;
        param_global.T_model = T_model;
        param_global.pgr = pgr;
        param_global.surv = surv;
        param_global.Mu2 = Mu2;
        param_global.demdist_2015 = demdist_2015;
        param_global.nz = nz;
        param_global.z = z;
        param_global.proddist = proddist;
        param_global.nt = nt;
        param_global.tr_z = tr_z;
        param_global.ndem = ndem;
        param_global.nk = nk;
        param_global.kgrid = kgrid;
        param_global.kub = kub;
        param_global.klb = klb;
        param_global.nb = nb;
        param_global.bub = bub;
        param_global.blb = blb;
        param_global.bgrid = bgrid;
        

    end
    
    function [] = generate_basedef_parameters()
        % This method generates the baseline definition parameters.  Baseline definition parameters include any parameters that change to define different baseline projections, 
        % such as any parameters used to calibrate to match savings and labor supply elasticities.
    end
    
    function [param_inctax, param_bustax] = generate_counterfactual_parameters(plan, save_results) %, param_bustax, param_socsec, param_immigration
        % This method generates the counterfactual parameters.  This includes any parameters that might change to define counterfactuals, as well as their baseline counterparts.
        param_global = paramGenerator.generate_global_parameters();
        mpci  = param_global.mpci;
        rpci  = param_global.rpci;
        nb    = param_global.nb;
        bgrid = param_global.bgrid;
        
        
        % Set baseline ("current law") to the default.
        if ~exist('plan','var')
            plan = 'base';
        end
        
        % Determine operating system to prevent xlsread in case of Linux.  Can be removed once we eliminate dependence on TPC inputs.
        operating_system = computer();
        if ~strcmp(operating_system,'PCWIN')&&~strcmp(operating_system,'PCWIN64')
            param_file = strcat('param_inctax_',plan);
            param_inctax = load(param_file);
            
            param_file = strcat('param_bustax_',plan);
            param_bustax = load(param_file);
            
        else
            display('Windows operating system detected. Reading TPC input Excel files and generating counterfactual parameters.')
        
            % Generate tax parameters

            %% Income Taxes

            % Reading taxable income ranges
            filename = 'TPC_Inputs.xlsx';
            sheet = 1;
            xlRange = 'A3:B622';
            income_range_bins = xlsread(filename,sheet,xlRange);   % Column 1 is lower bounds, Column 2 is upper bounds.

            % Reading tax rates.
            filename = 'TPC_Inputs.xlsx';
            sheet = 1;   % "Tax Function" sheet in Excel spreadsheet.
            if strcmp(plan,'base')
                xlRange = 'C3:C622';
            elseif strcmp(plan,'trump')    
                xlRange = 'E3:E622';
            elseif strcmp(plan,'ryan')
                xlRange = 'F3:F622';
            end
            effective_tax_rates = xlsread(filename,sheet,xlRange);    % Effective tax rates on taxable income

            % Estimating income tax function.
            n_testpoints = 1000;   % Determines how finely to discretize income "data" points for estimation.
            [X, limit] = estimate_inctax_function(income_range_bins, effective_tax_rates, n_testpoints); 



            %% Deductions

            filename = 'TPC_Inputs.xlsx';
            sheet = 2;  % "Deduction Function" sheet in Excel spreadsheet.
            if strcmp(plan,'base')
                xlRange = 'C3:C622';
            elseif strcmp(plan,'trump')
                xlRange = 'E3:E622';
            elseif strcmp(plan,'ryan')
                xlRange = 'F3:F622';
            end

            income_upper_bound = 2*10e5;
            income_upper_bound = min(income_range_bins(end,2), income_upper_bound);    % Enforcing upper bound provided by TPC.
            a_log = 1;
            b_log = log(income_upper_bound)/log(10);
            n_testpoints = 1000;
            discretized_income = logspace(a_log,b_log,n_testpoints);
            n = length(discretized_income);

            taxable_income_ratio = xlsread(filename,sheet,xlRange);
            deductions = zeros(1,n);
            for i1 = 1:n
                income_bin_index = find(discretized_income(i1)<=income_range_bins(:,2),1,'first');   % Returns integer corresponding to taxable income range provided by TPC.
                deductions(i1) = discretized_income(i1)*(1-taxable_income_ratio(income_bin_index));
            end


            regressors = [discretized_income.^0; discretized_income.^1; discretized_income.^(1/2)]';   % Note the functional for specification implied here: f(x) = b_1*x^0 + b_2*x^1 + b_3*x^(1/2)
            [coefficients,~,~,~,~] = regress(deductions',regressors);


            n_regs = min(size(regressors));    % Assumes that the number of variables < number of observations.  Must be true for estimation.
            % Below changes the estimates to the order of polynomial currently hard coded
            if n_regs<5
                coefs = [coefficients(2:n_regs)', zeros(1,4-(n_regs-1))];
            elseif n_regs>=5
                warning('Deduction polynomial order potential mismatch. Curtailing at 4.')
                coefs = coefficients(2:5);
            end

            avg_deduc = coefficients(1);


            param_inctax.avg_deduc = avg_deduc;
            param_inctax.coefs     = coefs;
            param_inctax.limit     = limit;
            param_inctax.X         = X;




            %% Business Taxes
            
            filename = 'TPC_Inputs.xlsx';
            sheet = 3;  % "Capital Tax" sheet in Excel spreadsheet.
            if strcmp(plan,'base')
                xlRange = 'B2:B5';
            elseif strcmp(plan,'trump')
                xlRange = 'D2:D5';
            elseif strcmp(plan,'ryan')
                xlRange = 'E2:E5';
            end

            bus_params    = xlsread(filename,sheet,xlRange);    % gives business parameters of policy proposal
            cap_tax_share = bus_params(2); 
            exp_share     = bus_params(4); 
            tau_cap       = bus_params(1); 
            tau_capgain   = 0; 
            
            param_bustax.cap_tax_share = cap_tax_share;
            param_bustax.exp_share     = exp_share;
            param_bustax.tau_cap       = tau_cap;
            param_bustax.tau_capgain   = tau_capgain;
            

            % Generating .mat files for future access by non-Windows OS's if running on Windows OS and save_results is true.
            if exist('save_results','var')
                if save_results
                    filename = ['param_inctax_' plan '.mat'];
                    save(filename,'avg_deduc','coefs','limit','X')

                    filename = ['param_bustax_' plan '.mat'];
                    save(filename,'cap_tax_share','exp_share','tau_cap','tau_capgain')
                end
            end
        end


        %% Generate Social Security parameters
        
        b_scale     = 1.6;    % Scaling social security benefits to match total outlays as a percentage of GDP.
        v_ss_max    = 118500*(mpci/rpci);  % ("Taxmax") Maximum taxable earnings as calculated for Social Security. Also, maximum income that counts towards average earnings for Social Security benefit calculation.
        ss_tax_cred = .15;  % Percentage tax credit on Social Security benefits received.
        tau_ss      = .124;    % Social Security tax (employer + employee portion).

        % Benefit formula parameters
        v1 = 856*12*(mpci/rpci); % first threshold for social security benefit calculation
        r1 = .9; % percentage of indexed ss benefit applied to the lowest bracket
        v2 = 5157*12*(mpci/rpci); % second threshold for social security benefit calculation
        r2 = .32; % percentage of indexed ss benefit applied to the second lowest bracket
        r3 = .15; % percentage of indexed ss benefit applied to the highest bracket
        
        % Calculating the benefit.  Note dependence of benefit calculation on "bgrid," which is calculated with the global parameters.
        ben = zeros(nb,1);
        for i1 = 1:nb
            if bgrid(i1)<=v1
                ben(i1) = r1*bgrid(i1);
            elseif (bgrid(i1)>v1)&&(bgrid(i1)<=v2)
                ben(i1) = r1*v1 + r2*(bgrid(i1)-v1);
            elseif bgrid(i1)>v2
                ben(i1) = r1*v1 + r2*(v2-v1) + r3*(bgrid(i1) - v2);
            end
        end
        
        NRA_baseline = 47;    % Baseline normal retirement age.
        
        
        
        
        
        
        
        %% Generate immigration parameters
    end
    
end

end





function [X, limit] = estimate_inctax_function(income_range_bins, effective_tax_rates, n_testpoints)

% Setting optimization options.
opts1 = optimset('TolFun',1e-10,'MaxIter',10e10,'MaxFunEvals',10e10);

% Generating variables for income tax estimation.
income_upper_bound = income_range_bins(end,2);
discretized_income = linspace(0,income_upper_bound,n_testpoints);

inctax_bills = zeros(1,n_testpoints);

for i1 = 1:n_testpoints
    income_bin_index = find(discretized_income(i1)<=income_range_bins(:,2),1,'first');   % returns integer corresponding to taxable income range provided by TPC
    inctax_bills(i1) = discretized_income(i1)*effective_tax_rates(income_bin_index);
end

[X,~] = fminsearch(@(phi) gouveiastrauss(phi,discretized_income(2:end),inctax_bills(2:end)),[.8 .01 .36], opts1);
limit = X(3);
X = X(1:2); 


function sum_sq = gouveiastrauss(phi,x1,y1)

    ytax = phi(3).*(x1 - (x1.^(-phi(1)) + phi(2)).^(-1/phi(1)));    % Calculate taxable income
    ybar = ytax./x1;    % Calculate effective income
    sum_sq = sum((ybar - y1./x1).^2);    % Fit function to effective income.
end

end





function [y, epsgrid] = generate_normdist_shocks(nt,sigt)

% This function gives the cutoff points for nt slices of a normal
% distribution with equal weight.

epsgrid = zeros(nt+1,1);    % creates grid for endpoints

for i1 = 1:nt+1
    epsgrid(i1) = sigt*norminv((i1-1)/nt,0,1);
end

y = zeros(nt,1);    % grid of conditional means
for i1 = 1:nt
    
    e1 = (epsgrid(i1))/sigt;
    e2 = (epsgrid(i1+1))/sigt;
    y(i1) = nt*sigt*(normpdf(e1,0,1) - normpdf(e2,0,1));
    
end

end






function [y, yprob, epsgrid] = generate_persistent_shocks(ny,const,lambda,sigep)

% This function creates a Markov matrix that approximates an AR(1) process.  Written for C-code generation.
% ny = gridsize
% const = constant on ar1 process
% lambda = coefficient on ar1 process
% sigep = standard deviation of error terms
% y(t) = mu*(1-lambda) + lambda(t-1) + eps
% recovering mu

mu = const/(1-lambda);
sigy = sigep/(sqrt(1-lambda^2));
epsgrid = zeros(ny+1,1);    % creates grid for endpoints

for i1 = 1:ny+1    
    epsgrid(i1) = sigy*norminv((i1-1)/ny,0,1) + mu;
end

y = zeros(ny,1);    % grid of conditional means
for i1 = 1:ny    
    e1 = (epsgrid(i1)-mu)/sigy;
    e2 = (epsgrid(i1+1)-mu)/sigy;
    y(i1) = ny*sigy*(normpdf(e1,0,1) - normpdf(e2,0,1)) + mu;
end

yprob = zeros(ny,ny);
for i1 = 1:ny
    for i2 = 1:ny
        ei1 = epsgrid(i1);
        ei2 = epsgrid(i1+1);
        ej1 = epsgrid(i2);
        ej2 = epsgrid(i2+1);
        integrand([],ej1,ej2,sigy,sigep,mu,lambda);
        yprob(i1,i2) = (ny/(sqrt(2*pi*(sigy^2))))*quadgk(@integrand,ei1,ei2);        
    end    
end

end

function c2 =  integrand(x,ej1_,ej2_,sigy_,sigep_,mu_,lambda_)

% Enforce function inlining for C code generation
coder.inline('always');

% Define parameters as persistent variables
persistent initialized
persistent ej1
persistent ej2
persistent sigy
persistent sigep
persistent mu
persistent lambda

% Initialize parameters
if isempty(initialized)
    ej1    = 0;
    ej2    = 0;
    sigy   = 0;
    sigep  = 0;
    mu     = 0;
    lambda = 0;
    
    initialized = true;
    
end

% Set parameters if provided
if (nargin > 1)
    ej1    = ej1_;
    ej2    = ej2_;
    sigy   = sigy_;
    sigep  = sigep_;
    mu     = mu_;
    lambda = lambda_;
    
    c2 = [];
    return
end


nx = length(x);

c2 = zeros(1,nx);

for i1 = 1:nx

    term1 = exp(-((x(i1)-mu)^2)/(2*sigy^2));
    term2 = normcdf((ej2-mu*(1-lambda)-lambda*x(i1))/sigep,0,1);
    term3 = normcdf((ej1-mu*(1-lambda)-lambda*x(i1))/sigep,0,1);
    c2(i1) = term1*(term2-term3);
    
end

end


function surv = calculate_survival()

   surv = [0.000995136826434
           0.001030142563270
           0.001055052308246
           0.001070116441161
           0.001080337982679
           0.001086289739767
           0.001095904151824
           0.001127786393148
           0.001191607772250
           0.001279286761967
           0.001378242329718
           0.001475563110369
           0.001572042255373
           0.001663319298670
           0.001752412300092
           0.001849547284972
           0.001952489455367
           0.002066679395856
           0.002187805459109
           0.002317148776168
           0.002464180462455
           0.002619648409666
           0.002776098413999
           0.002927963005853
           0.003083332812207
           0.003260329581086
           0.003466443203404
           0.003710291930733
           0.004003007728613
           0.004327847498427
           0.004708260459706
           0.005127784084076
           0.005583927273047
           0.006069927248038
           0.006601392611336
           0.007188574767233
           0.007847959184511
           0.008585680946046
           0.009413891314214
           0.010327313510670
           0.011335560808543
           0.012422410391424
           0.013571516913235
           0.014773886553643
           0.016061597313524
           0.017437456348170
           0.018998718672376
           0.020621187057936
           0.022369112545660
           0.024208204553687
           0.026249942734302
           0.028606587513248
           0.031138947851687
           0.033805452476018
           0.036698997126841
           0.039944266267798
           0.043607996645804
           0.047653358078815
           0.051920188928337
           0.056656795889449
           0.062039388800224
           0.067849890542353
           0.074271414162276
           0.081490161948741
           0.089142356885846
           0.097846896816678
           0.107220914541972
           0.117425108249143
           0.128481061773782
           0.140342065061076
           0.153351006067180
           0.167159644542377
           0.181832943589499
           0.197560459551117
           0.214203242042907
           0.231015900125199
           0.248233084037350
           0.264967838478304
           0.281384286208328
           0
           0                ];
       
       surv = 1 - surv';
                   
end

