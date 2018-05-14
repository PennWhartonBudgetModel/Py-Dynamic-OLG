%%
% Report generator for dynamic model parameters and outputs.
%
%%

wipe;

%% Params for high return economy
% params.beta =  0.976000000000000;
% params.gamma = 0.779000000000000;
% params.sigma = 1.600000000000000;
% params.modelunit_dollar = 3.952983000000000e-05;
% params.IsLowReturn = false;

%% Params for low return economy
params.beta =  0.994000000000000;
params.gamma = 0.726000000000000;
params.sigma = 1.500000000000000;
params.modelunit_dollar = 3.850000000000000e-05;
params.IsLowReturn = true;

params.TransitionFirstYear = 2018;
params.TransitionLastYear  = 2018+25;

%% Scenarios and directories

% Solve for baseline steady state
scenario   = Scenario(struct('economy', 'closed', 'beta', params.beta, 'gamma', params.gamma, ...
                             'sigma', params.sigma, 'modelunit_dollar', params.modelunit_dollar, ...
                             'IsLowReturn', params.IsLowReturn, 'bequest_phi_1', 0, 'TaxCode', 'PermanentTCJA', ...
                             'TransitionFirstYear', params.TransitionFirstYear, 'TransitionLastYear', params.TransitionLastYear));

sc_steady   = scenario.currentPolicy.steady;
sc_open     = scenario.open();
sc_closed   = scenario.closed();
scb_open    = scenario.currentPolicy.open();
scb_closed  = scenario.currentPolicy.closed();

steady_dir  = PathFinder.getWorkingDir(sc_steady);
open_dir    = PathFinder.getWorkingDir(sc_open);
closed_dir  = PathFinder.getWorkingDir(sc_closed);
bopen_dir   = PathFinder.getWorkingDir(scb_open);
bclosed_dir = PathFinder.getWorkingDir(scb_closed);

%% Solve model
ModelSolver.solve(scenario);

%% Build variables

% Import variables in dynamics file
d_steady = load(fullfile(steady_dir  , 'dynamics.mat'), ... 
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits', 'pits', 'ssts', 'bens', 'investment', 'corpTaxs' );
d_steady.caps_foreign = 0;
d_open   = load(fullfile(open_dir  , 'dynamics.mat'), ... 
            'debts', 'revs', 'outs', 'caps', 'caps_foreign', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits', 'pits', 'ssts', 'bens', 'investment', 'corpTaxs' );
d_closed = load(fullfile(closed_dir, 'dynamics.mat'), ...
            'debts', 'revs', 'outs', 'caps', 'caps_foreign', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits', 'pits', 'ssts', 'bens', 'investment', 'corpTaxs' );

% Import variables in **baseline** dynamics file
b_open   = load(fullfile(bopen_dir  , 'dynamics.mat'), ... 
            'debts', 'revs', 'outs', 'caps', 'caps_foreign', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits', 'pits', 'ssts', 'bens', 'investment', 'corpTaxs' );
b_closed = load(fullfile(bclosed_dir, 'dynamics.mat'), ...
            'debts', 'revs', 'outs', 'caps', 'caps_foreign', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits', 'pits', 'ssts', 'bens', 'investment', 'corpTaxs' );

% Import variables in statics file
s_open   = load(fullfile(open_dir  , 'statics.mat' ), ...
            'debts', 'revs', 'outs', 'caps', 'caps_foreign', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits', 'pits', 'ssts', 'bens', 'investment', 'corpTaxs');
s_closed = load(fullfile(closed_dir, 'statics.mat' ), ...
            'debts', 'revs', 'outs', 'caps', 'caps_foreign', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits', 'pits', 'ssts', 'bens', 'investment', 'corpTaxs');

% Import variables in market file
m_steady  = load(fullfile(steady_dir, 'market.mat'  ), 'caprates', 'wages', 'rhos', 'debtrates', 'equityFundDividends', 'equityFundPrices', 'capshares_1');
m_open    = load(fullfile(open_dir  , 'market.mat'  ), 'caprates', 'wages', 'rhos', 'debtrates', 'equityFundDividends', 'equityFundPrices', 'capshares_1');
mb_open   = load(fullfile(bopen_dir  , 'market.mat' ), 'caprates', 'wages', 'rhos', 'debtrates', 'equityFundDividends', 'equityFundPrices', 'capshares_1');
m_closed  = load(fullfile(closed_dir, 'market.mat'  ), 'caprates', 'wages', 'rhos', 'debtrates', 'equityFundDividends', 'equityFundPrices', 'capshares_1');
mb_closed = load(fullfile(bclosed_dir, 'market.mat' ), 'caprates', 'wages', 'rhos', 'debtrates', 'equityFundDividends', 'equityFundPrices', 'capshares_1');

% Include steady state year in all series
include_ss = @(x, xss) [xss x];
os = {'debts', 'revs', 'outs', 'caps', 'caps_foreign', 'labs', 'labeffs', 'labincs', 'assets', 'cons', 'pops', 'labpits', 'caprevs', 'cits', 'pits', 'ssts', 'bens', 'investment', 'corpTaxs'};
for o = os
    d_open.(o{1})   = include_ss(d_open.(o{1})  , d_steady.(o{1}));
    s_open.(o{1})   = include_ss(s_open.(o{1})  , d_steady.(o{1}));
    b_open.(o{1})   = include_ss(b_open.(o{1})  , d_steady.(o{1}));
    d_closed.(o{1}) = include_ss(d_closed.(o{1}), d_steady.(o{1}));
    s_closed.(o{1}) = include_ss(s_closed.(o{1}), d_steady.(o{1}));
    b_closed.(o{1}) = include_ss(b_closed.(o{1}), d_steady.(o{1}));
end
os = {'caprates', 'wages', 'rhos', 'debtrates', 'equityFundDividends', 'equityFundPrices', 'capshares_1'};
for o = os
    m_open.(o{1})    = include_ss(m_open.(o{1})  , m_steady.(o{1}));
    mb_open.(o{1})   = include_ss(mb_open.(o{1})  , m_steady.(o{1}));
    m_closed.(o{1})  = include_ss(m_closed.(o{1}), m_steady.(o{1}));
    mb_closed.(o{1}) = include_ss(mb_closed.(o{1}), m_steady.(o{1}));
end

% Create new variables

m_open.totrates = m_open.capshares_1 .* m_open.equityFundDividends + (1 - m_open.capshares_1) .* m_open.debtrates;
mb_open.totrates = mb_open.capshares_1 .* mb_open.equityFundDividends + (1 - mb_open.capshares_1) .* mb_open.debtrates;
m_closed.totrates = m_closed.capshares_1 .* m_closed.equityFundDividends + (1 - m_closed.capshares_1) .* m_closed.debtrates;
mb_closed.totrates = mb_closed.capshares_1 .* mb_closed.equityFundDividends + (1 - mb_closed.capshares_1) .* mb_closed.debtrates;

for var_name = {'tax', 'totinc', 'totincwss'}
    d_open.(var_name{1})   = build_var(var_name{1}, d_open, m_open, sc_open, 0);
    s_open.(var_name{1})   = build_var(var_name{1}, s_open, mb_open, sc_open, 1);
    d_closed.(var_name{1}) = build_var(var_name{1}, d_closed, m_closed, sc_closed, 0);
    s_closed.(var_name{1}) = build_var(var_name{1}, s_closed, mb_closed, sc_closed, 1);
end

%% Tables

% Years vector
yearsv = [(params.TransitionFirstYear - 1):1:(params.TransitionLastYear-1)];

%% Raw values reports
header      = {'year', 'MPK', 'gov_interest_rate', 'equityFundDividends', 'equityFundPrices', 'wages', 'K_L_ratio', ...
               'output_d', 'capital_d', 'labor_d', 'labor_efficient_d', 'labor_efficient_by_pop_d', 'labpits_d', 'labincs_d', 'caprevs_d', 'consumption_d', 'assets_d', 'debt_d', 'revenues_d', 'investment_d', 'corpTaxs_d', ...
               'output_s', 'capital_s', 'labor_s', 'labor_efficient_s', 'labor_efficient_by_pop_s', 'labpits_s', 'labincs_s', 'caprevs_s', 'consumption_s', 'assets_s', 'debt_s', 'revenues_s', 'investment_s', 'corpTaxs_s', ...
               'MPK_netDep_b', 'wages_b', 'output_b', 'capital_b', 'labor_efficient_b', 'labincs_b', 'revenues_b', 'investment_b', 'corpTaxs_b'};

f = @( d, s, m, mb, b ) table(yearsv', m.caprates', m.debtrates', m.equityFundDividends', m.equityFundPrices', m.wages', m.rhos', ...
                  d.outs', d.caps', d.labs', d.labeffs', (d.labeffs ./ d.pops)', d.labpits', d.labincs', d.caprevs', d.cons', d.assets', d.debts', d.revs', d.investment', d.corpTaxs', ...
                  s.outs', s.caps', s.labs', s.labeffs', (s.labeffs ./ s.pops)', s.labpits', s.labincs', s.caprevs', s.cons', s.assets', s.debts', s.revs', s.investment', s.corpTaxs', ...
                  mb.caprates', mb.wages', b.outs', b.caps', b.labeffs', b.labincs', b.revs', b.investment', b.corpTaxs', ...
                  'VariableNames', header);
           
res_open   = f(d_open, s_open, m_open, mb_open, b_open);
res_closed = f(d_closed, s_closed, m_closed, mb_closed, b_closed);

open('res_closed')
open('res_open')
writetable(res_open, 'res_open.csv')
writetable(res_closed, 'res_closed.csv')
           
%% Deltas report
header = {'year', ...
          'MPK', 'gov_interest_rate', 'equityFundDividends', 'equityFundPrices', 'wages', 'KLratio', ...
          'output', 'capital', 'caps_foreign', 'labor', ...
          'labor_efficient', 'assets', 'consumption', 'investment' ...
          'debt', 'revenue', 'tax', ...
          'labpits', 'caprevs', 'ssts', 'cits', 'pits', 'corpTaxs',  ...
          'totinc', 'totincwss', ...
          'totrates', 'MPK_b', 'gov_interest_rate_b', 'equityFundDividends_b', 'equityFundPrices_b', 'wages_b', 'KLratio_b', 'totrates_b'};
      
f = @( m, mb, d, s ) table(yearsv', ...
                    m.caprates', m.debtrates', m.equityFundDividends', m.equityFundPrices', m.wages', m.rhos', ...
                   (d.outs ./ s.outs)', (d.caps ./ s.caps)', (d.caps_foreign ./ s.caps_foreign)', (d.labs ./ s.labs)', ...
                   (d.labeffs ./ s.labeffs)', (d.assets ./ s.assets)', (d.cons ./ s.cons)', (d.investment ./ s.investment)', ...
                   (d.debts ./ s.debts)', (d.revs ./ s.revs)', (d.tax ./ s.tax)', ...
                   (d.labpits ./ s.labpits)', (d.caprevs ./ s.caprevs)', ...
                   (d.ssts ./ s.ssts)', (d.cits ./ s.cits)', (d.pits ./ s.pits)', (d.corpTaxs ./ s.corpTaxs)',...
                   (d.totinc ./ s.totinc)', (d.totincwss ./ s.totincwss)', ...
                    m.totrates', mb.caprates', mb.debtrates', mb.equityFundDividends', mb.equityFundPrices', mb.wages', mb.rhos', mb.totrates', ...
                   'VariableNames', header);
    
delta_closed  = f(m_closed, mb_closed, d_closed, s_closed);
delta_open    = f(m_open  , mb_open  , d_open  , s_open  );
open('delta_closed');
open('delta_open');

if (scenario.IsLowReturn)
    filename = 'deltas_lo_return.xlsx';
else
    filename = 'deltas_hi_return.xlsx';
end
writetable(delta_open  ,filename,'Sheet','open'  )
writetable(delta_closed,filename,'Sheet','closed')


%% TPC report
header = {'year', 'MPK', 'gov_interest_rate', 'equityFundDividends', 'equityFundPrices', 'wages', 'KLratio', ...
          'output', 'debt', 'revenue', 'capital', 'labor', 'labor_efficient', 'assets', 'assets_per_pop'};
f = @( d, s, m ) table(yearsv', m.caprates', m.debtrates', m.equityFundDividends', m.equityFundPrices', m.wages', m.rhos', ...
                      (d.outs ./ s.outs)', (d.debts ./ s.debts)', (d.revs ./ s.revs)', ...
                      (d.caps ./ s.caps)', (d.labs ./ s.labs)', (d.labeffs ./ s.labeffs)', ...
                      (d.assets ./ s.assets)', ((d.assets./d.pops) ./ (s.assets./s.pops))', ...
                      'VariableNames', header);
    
TPC_closed  = f(d_closed, s_closed, m_closed);
TPC_open    = f(d_open  , s_open  , m_open  );

if (scenario.IsLowReturn)
    filename = 'TPC_report_lo_return.xlsx';
else
    filename = 'TPC_report_hi_return.xlsx';
end
writetable(TPC_open  ,filename,'Sheet','open')
writetable(TPC_closed,filename,'Sheet','closed')
                 

%% FUNCTIONS

function x = build_var(x_name, dynamic_struct, market_struct, scenario, static)

    if strcmp(x_name, 'tax')
        
        x = dynamic_struct.pits + dynamic_struct.cits + dynamic_struct.ssts;
        % Recall that cits now mean capital income taxes paid by households
        % (so cits = cits_domestic). Corporate taxes, paid by firms, are
        % denoted by corpTaxs.
        
    elseif strcmp(x_name, 'totinc')
        
        % Define time constants
        s = ParamGenerator.timing( scenario );
        T_life = s.T_life; T_model = s.T_model;
        % Define assets grid and dimensions
        s = ParamGenerator.grids( scenario );
        nz = s.nz; nk = s.nk; nb = s.nb;
        karray = repmat(reshape(s.kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model+1]);
        % Define tax
        s = ParamGenerator.tax( scenario );
        shareCapitalCorporate = s.shareCapitalCorporate;
        s = ParamGenerator.tax( scenario.currentPolicy.steady );
        shareCapitalCorporate = [s.shareCapitalCorporate shareCapitalCorporate'];
        DIST = zeros(nz,nk,nb,T_life,T_model+1);
        s = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy.steady), 'distribution.mat'), 'DIST');
        DIST(:,:,:,:,1,:) = reshape(sum(s.DIST, 5), [nz,nk,nb,T_life,1]);

        if static
            s = load(fullfile(PathFinder.getWorkingDir(scenario), 'Static_distribution.mat'), 'Static_DIST');
            DIST(:,:,:,:,2:end,:) = reshape(sum(s.Static_DIST, 5), [nz,nk,nb,T_life,T_model]);
        else
            s = load(fullfile(PathFinder.getWorkingDir(scenario), 'distribution.mat'), 'DIST');
            DIST(:,:,:,:,2:end,:) = reshape(sum(s.DIST, 5), [nz,nk,nb,T_life,T_model]);
        end
        
        f = @(F) sum(sum(reshape(DIST .* F, [], T_model+1), 1), 3);
        x = dynamic_struct.labincs + f(karray) .* (1 - shareCapitalCorporate) .* market_struct.debtrates ...
                                   + f(karray) .* shareCapitalCorporate .* market_struct.equityFundDividends;
        
    elseif strcmp(x_name, 'totincwss')
        
        x = dynamic_struct.totinc + dynamic_struct.bens;
        
    else
        
        error('Variable does not exist.')
        
    end

end
