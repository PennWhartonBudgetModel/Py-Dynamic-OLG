%%
% Static aggregate tax calculator.  Adapted from dynamic optimization solver.
% 
%%


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


% Define argument properties for C code generation
T_max   = 100;
Tss_max =  75;
nb_max  =  50;
nk_max  =  50;
nz_max  =  50;

assert( isa(T,              'double') && (size(T,               1) == 1         ) && (size(T,               2) == 1         ) );
assert( isa(Tr,             'double') && (size(Tr,              1) == 1         ) && (size(Tr,              2) == 1         ) );
assert( isa(Tss,            'double') && (size(Tss,             1) == 1         ) && (size(Tss,             2) == 1         ) );
assert( isa(birthyear,      'double') && (size(birthyear,       1) == 1         ) && (size(birthyear,       2) == 1         ) );

assert( isa(kgrid,          'double') && (size(kgrid,           1) <= nk_max    ) && (size(kgrid,           2) == 1         ) );
assert( isa(nb,             'double') && (size(nb,              1) == 1         ) && (size(nb,              2) == 1         ) );
assert( isa(nk,             'double') && (size(nk,              1) == 1         ) && (size(nk,              2) == 1         ) );
assert( isa(nz,             'double') && (size(nz,              1) == 1         ) && (size(nz,              2) == 1         ) );
assert( isa(idem,           'double') && (size(idem,            1) == 1         ) && (size(idem,            2) == 1         ) );

assert( isa(mpci,           'double') && (size(mpci,            1) == 1         ) && (size(mpci,            2) == 1         ) );
assert( isa(rpci,           'double') && (size(rpci,            1) == 1         ) && (size(rpci,            2) == 1         ) );
assert( isa(cap_tax_share,  'double') && (size(cap_tax_share,   1) == 1         ) && (size(cap_tax_share,   2) == 1         ) );
assert( isa(ss_tax_cred,    'double') && (size(ss_tax_cred,     1) == 1         ) && (size(ss_tax_cred,     2) == 1         ) );

assert( isa(tau_cap,        'double') && (size(tau_cap,         1) == 1         ) && (size(tau_cap,         2) == 1         ) );
assert( isa(tau_capgain,    'double') && (size(tau_capgain,     1) == 1         ) && (size(tau_capgain,     2) == 1         ) );
assert( isa(ss_tax,         'double') && (size(ss_tax,          1) == 1         ) && (size(ss_tax,          2) <= Tss_max   ) );
assert( isa(taxmax,         'double') && (size(taxmax,          1) == 1         ) && (size(taxmax,          2) <= Tss_max   ) );
assert( isa(z,              'double') && (size(z,               1) <= nz_max    ) && (size(z,               2) <= T_max     ) && (size(z, 3) == 2) );

assert( isa(wages,          'double') && (size(wages,           1) == 1         ) && (size(wages,           2) <= Tss_max   ) );
assert( isa(cap_shares,     'double') && (size(cap_shares,      1) == 1         ) && (size(cap_shares,      2) <= Tss_max   ) );
assert( isa(debt_shares,    'double') && (size(debt_shares,     1) == 1         ) && (size(debt_shares,     2) <= Tss_max   ) );
assert( isa(rate_caps,      'double') && (size(rate_caps,       1) == 1         ) && (size(rate_caps,       2) <= Tss_max   ) );
assert( isa(rate_govs,      'double') && (size(rate_govs,       1) == 1         ) && (size(rate_govs,       2) <= Tss_max   ) );
assert( isa(exp_subsidys,   'double') && (size(exp_subsidys,    1) == 1         ) && (size(exp_subsidys,    2) <= Tss_max   ) );
assert( isa(q_tobin,        'double') && (size(q_tobin,         1) == 1         ) && (size(q_tobin,         2) == 1         ) );
assert( isa(q_tobin0,       'double') && (size(q_tobin0,        1) == 1         ) && (size(q_tobin0,        2) == 1         ) );

assert( isa(avg_deduc,      'double') && (size(avg_deduc,       1) == 1         ) && (size(avg_deduc,       2) == 1         ) );
assert( isa(clinton,        'double') && (size(clinton,         1) == 1         ) && (size(clinton,         2) == 1         ) );
assert( isa(coefs,          'double') && (size(coefs,           1) == 1         ) && (size(coefs,           2) <= 10        ) );
assert( isa(limit,          'double') && (size(limit,           1) == 1         ) && (size(limit,           2) == 1         ) );
assert( isa(X,              'double') && (size(X,               1) == 1         ) && (size(X,               2) <= 10        ) );
assert( isa(ss_benefit,     'double') && (size(ss_benefit,      1) <= nb_max    ) && (size(ss_benefit,      2) <= Tss_max   ) );

assert( isa(labopt, 'double') && all(size(labopt) <= [nk_max, nz_max, nb_max, T_max]) );


%% Initialization

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