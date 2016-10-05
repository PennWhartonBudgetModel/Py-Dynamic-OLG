% Jorge | 2016-09-29
% 
% Static aggregate tax calculation function.  Adapted from dynamic optimization solver.
% 
% 


function [fedincome, fedincomess, fitax, fitaxss, fsstax, fsstaxss, fcaptax, fcaptaxss] ...
          ...
          = calculate_static_taxes(...
              T, Tr, Tss, birthyear, ...
              kgrid, nb, nk, nz, idem, ...
              mpci, rpci, cap_tax_share, ss_tax_cred, ...
              tau_cap, tau_capgain, ss_tax, taxmax, z, ...
              wages, cap_shares, debt_shares, rate_caps, rate_govs, exp_subsidys, q_tobin, q_tobin0, ...
              avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
              labopt) %#codegen


% Define past years and effective retirement age
% (birthyear is defined relative to the present and can take both positive and negative values)
T_past = -min(0, birthyear);
Tr_eff = max(Tr, T_past);

% Find effective numbers of working and retirement years
N_w = Tr_eff - T_past;
N_r = T      - Tr_eff;

% Initialize arrays
fedincome   = zeros(nk,nz,nb,N_w);
fitax       = zeros(nk,nz,nb,N_w);
fsstax      = zeros(nk,nz,nb,N_w);
fcaptax     = zeros(nk,nz,nb,N_w);

fedincomess = zeros(nk,nb,N_r);
fitaxss     = zeros(nk,nb,N_r);
fsstaxss    = zeros(nk,nb,N_r);
fcaptaxss   = zeros(nk,nb,N_r);


for t = N_r:-1:1

    % Determine age and year, bounded by projection period
    age  = t + Tr_eff;
    year = max(1, min(age + birthyear, Tss));

    % Extract annual parameters
    ben         = ss_benefit(:,year);
    rate_cap    = rate_caps(year);
    rate_gov    = rate_govs(year);
    cap_share   = cap_shares(year);
    debt_share  = debt_shares(year);
    exp_subsidy = exp_subsidys(year); 

    for ib = 1:nb
        for ik = 1:nk

            % Calculate tax terms
            fincome_lab = (1-ss_tax_cred)*ben(ib);
            [fincome, ftax, ~, fcap] ...
                = calculate_tax(fincome_lab, kgrid(ik), ...
                                cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy, ...
                                avg_deduc, coefs, limit, X, mpci, rpci, 0, 0, clinton, year, q_tobin, q_tobin0);

            % Record values
            fedincomess(ik,ib,t) = fincome;
            fitaxss    (ik,ib,t) = ftax;
            fsstaxss   (ik,ib,t) = 0;
            fcaptaxss  (ik,ib,t) = fcap;

        end
    end

end

% (0 iterations if T_past >= Tr -- i.e. if there is no working period)
for t = N_w:-1:1

    % Determine age and year, bounded by projection period
    age  = t + T_past;
    year = max(1, min(age + birthyear, Tss));

    % Extract annual parameters
    v_ss_max    = taxmax(year);
    tau_ss      = ss_tax(year);
    wage        = wages(year);
    rate_cap    = rate_caps(year);
    rate_gov    = rate_govs(year);
    cap_share   = cap_shares(year);
    debt_share  = debt_shares(year);
    exp_subsidy = exp_subsidys(year);

    for ib = 1:nb
        for iz = 1:nz

            eff_wage = wage*max(z(iz,age,idem), 0);

            for ik = 1:nk

                % Calculate tax terms
                fincome_lab = eff_wage*labopt(ik,iz,ib,t);
                [fincome, ftax, sstax, fcap] ...
                    = calculate_tax(fincome_lab, kgrid(ik), ...
                                    cap_share, rate_cap, debt_share, rate_gov, cap_tax_share, tau_cap, tau_capgain, exp_subsidy, ...
                                    avg_deduc, coefs, limit, X, mpci, rpci, tau_ss, v_ss_max, clinton, year, q_tobin, q_tobin0);

                % Record values
                fedincome(ik,iz,ib,t) = fincome;
                fitax    (ik,iz,ib,t) = ftax;
                fsstax   (ik,iz,ib,t) = sstax;
                fcaptax  (ik,iz,ib,t) = fcap;

            end

        end
    end

end

end