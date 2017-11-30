%%
% Report generator for dynamic model parameters and outputs.
%
%%

wipe;
% params = ParamGenerator.invert(struct('depreciation', 0.056, 'labelas', 0.5));

% Params for K/Y=3 with d=0.056
% params.beta =  0.98600000;
% params.gamma = 0.75000000;
% params.sigma = 1.24000000;
% params.modelunit_dollar = 4.359874681178362e-05;
% params.depreciation = 0.056;

% Params for K/Y=3 with d=0.08
params.beta =  1.003341000000000;
params.gamma = 0.680000000000000;
params.sigma = 1.500000000000000;
params.modelunit_dollar = 4.135682750000000e-05;
params.depreciation = 0.08;

% Solve for baseline steady state
scenario   = Scenario(struct('economy', 'closed', 'beta', params.beta, 'gamma', params.gamma, ...
                             'sigma', params.sigma, 'modelunit_dollar', params.modelunit_dollar, ...
                             'depreciation', params.depreciation, 'bequest_phi_1', 0, 'base_brackets', 'SenCMA'));
sc_steady  = scenario.currentPolicy.steady;
steady_dir = PathFinder.getWorkingDir(sc_steady);

% Solve for counterfactuals
sc_open     = scenario.open();
sc_closed   = scenario.closed();
scb_open    = scenario.currentPolicy.open();
scb_closed  = scenario.currentPolicy.closed();
open_dir    = PathFinder.getWorkingDir(sc_open);
closed_dir  = PathFinder.getWorkingDir(sc_closed);
bopen_dir   = PathFinder.getWorkingDir(scb_open);
bclosed_dir = PathFinder.getWorkingDir(scb_closed);

ModelSolver.solve(scenario);

% Import variables in dynamics file
d_steady = load(fullfile(steady_dir  , 'dynamics.mat'), ... 
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits', 'pits', 'ssts', 'bens' );
d_steady.cits_domestic = d_steady.cits;
d_open   = load(fullfile(open_dir  , 'dynamics.mat'), ... 
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits_domestic', 'pits', 'ssts', 'bens' );

d_closed = load(fullfile(closed_dir, 'dynamics.mat'), ...
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits_domestic', 'pits', 'ssts', 'bens' );

% Import variables in **baseline** dynamics file
b_open   = load(fullfile(bopen_dir  , 'dynamics.mat'), ... 
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits_domestic', 'pits', 'ssts', 'bens' );

b_closed = load(fullfile(bclosed_dir, 'dynamics.mat'), ...
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits_domestic', 'pits', 'ssts', 'bens' );

% Import variables in statics file
s_open   = load(fullfile(open_dir  , 'statics.mat' ), ...
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits_domestic', 'pits', 'ssts', 'bens');
s_closed = load(fullfile(closed_dir, 'statics.mat' ), ...
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'labincs', 'assets', 'pops', 'cons', 'labpits', 'caprevs', 'cits_domestic', 'pits', 'ssts', 'bens');

% Import variables in market file
m_steady = load(fullfile(steady_dir, 'market.mat' ), 'caprates', 'wages', 'rhos', 'kpricescale', 'govrates', 'totrates', 'qtobin');
m_open   = load(fullfile(open_dir  , 'market.mat' ), 'caprates', 'wages', 'rhos', 'kpricescale', 'govrates', 'totrates', 'qtobin');
mb_open   = load(fullfile(bopen_dir  , 'market.mat' ), 'caprates', 'wages', 'rhos', 'kpricescale', 'govrates', 'totrates', 'qtobin');
m_closed = load(fullfile(closed_dir, 'market.mat' ), 'caprates', 'wages', 'rhos', 'kpricescale', 'govrates', 'totrates', 'qtobin');
mb_closed = load(fullfile(bclosed_dir, 'market.mat' ), 'caprates', 'wages', 'rhos', 'kpricescale', 'govrates', 'totrates', 'qtobin');

% Reshape qtobin
m_open.qtobin   = repmat(m_open.qtobin  , [size(m_open.rhos  )]);
m_closed.qtobin = repmat(m_closed.qtobin, [size(m_closed.rhos)]);

% Include steady state year in all series
include_ss = @(x, xss) [xss x];
os = {'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'labincs', 'assets', 'cons', 'pops', 'labpits', 'caprevs', 'cits_domestic', 'pits', 'ssts', 'bens'};
for o = os
    d_open.(o{1})   = include_ss(d_open.(o{1})  , d_steady.(o{1}));
    s_open.(o{1})   = include_ss(s_open.(o{1})  , d_steady.(o{1}));
    b_open.(o{1})   = include_ss(b_open.(o{1})  , d_steady.(o{1}));
    d_closed.(o{1}) = include_ss(d_closed.(o{1}), d_steady.(o{1}));
    s_closed.(o{1}) = include_ss(s_closed.(o{1}), d_steady.(o{1}));
    b_closed.(o{1}) = include_ss(b_closed.(o{1}), d_steady.(o{1}));
end
os = {'caprates', 'wages', 'rhos', 'kpricescale', 'govrates', 'totrates', 'qtobin'};
for o = os
    m_open.(o{1})   = include_ss(m_open.(o{1})  , m_steady.(o{1}));
    mb_open.(o{1})   = include_ss(mb_open.(o{1})  , m_steady.(o{1}));
    m_closed.(o{1}) = include_ss(m_closed.(o{1}), m_steady.(o{1}));
    mb_closed.(o{1}) = include_ss(mb_closed.(o{1}), m_steady.(o{1}));
end

% Create new variables
for var_name = {'tax', 'totinc', 'totincwss'}
    d_open.(var_name{1}) = build_var(var_name{1}, d_open, m_open, sc_open, 0);
    s_open.(var_name{1}) = build_var(var_name{1}, s_open, m_open, sc_open, 1);
    d_closed.(var_name{1}) = build_var(var_name{1}, d_closed, m_closed, sc_closed, 0);
    s_closed.(var_name{1}) = build_var(var_name{1}, s_closed, m_closed, sc_closed, 1);
end

% Years vector
s = ParamGenerator.timing(sc_closed);
yearsv = [(s.first_transition_year - 1):1:(s.first_transition_year -1 + s.T_model)];

% Raw values reports
header      = {'year', 'MPK_netDep', 'gov_interest_rate', 'total_interest_rate', 'qtobin', 'wages', 'K_L_ratio', ...
               'output_d', 'capital_d', 'labor_d', 'labor_efficient_d', 'labor_efficient_by_pop_d', 'labpits_d', 'labincs_d', 'caprevs_d', 'consumption_d', 'assets_d', 'debt_d', 'revenues_d', ...
               'output_s', 'capital_s', 'labor_s', 'labor_efficient_s', 'labor_efficient_by_pop_s', 'labpits_s', 'labincs_s', 'caprevs_s', 'consumption_s', 'assets_s', 'debt_s', 'revenues_s', ...
               'MPK_netDep_b', 'wages_b', 'output_b', 'capital_b', 'labor_efficient_b', 'labincs_b', 'revenues_b'};

f = @( d, s, m, mb, b ) table(yearsv', m.caprates', m.govrates', m.totrates', m.qtobin', m.wages', m.rhos', ...
                  d.outs', d.caps', d.labs', d.labeffs', (d.labeffs ./ d.pops)', d.labpits', d.labincs', d.caprevs', d.cons', d.assets', d.debts', d.revs', ...
                  s.outs', s.caps', s.labs', s.labeffs', (s.labeffs ./ s.pops)', s.labpits', s.labincs', s.caprevs', s.cons', s.assets', s.debts', s.revs', ...
                  mb.caprates', mb.wages', b.outs', b.caps', b.labeffs', b.labincs', b.revs', ...
                  'VariableNames', header);
           
res_open   = f(d_open, s_open, m_open, mb_open, b_open);
res_closed = f(d_closed, s_closed, m_closed, mb_closed, b_closed);

open('res_closed')
open('res_open')
writetable(res_open, 'res_open.csv')
writetable(res_closed, 'res_closed.csv')
           
% Deltas report
header = {'year', 'output', 'debt', 'revenue', 'capital', 'labor', 'labor_efficient', 'labpits', 'caprevs', 'assets', 'consumption', 'pits', 'tax', 'totinc', 'totincwss'};
f = @( d, s ) table(yearsv', (d.outs ./ s.outs)', (d.debts ./ s.debts)', (d.revs ./ s.revs)',   ...
         (d.caps ./ s.caps)', (d.labs ./ s.labs)', (d.labeffs ./ s.labeffs)' , (d.labpits ./ s.labpits)' ,   ...
         (d.caprevs ./ s.caprevs)', (d.assets ./ s.assets)', (d.cons ./ s.cons)', (d.pits ./ s.pits)',...
         (d.tax ./ s.tax)', (d.totinc ./ s.totinc)', (d.totincwss ./ s.totincwss)', 'VariableNames', header);
    
delta_closed  = f(d_closed, s_closed);
delta_open    = f(d_open  , s_open  );

open('delta_closed');
open('delta_open');


% TPC report
header = {'year', 'MPK_netDep', 'gov_interest_rate', 'total_interest_rate', 'qtobin', 'wages', 'K_L_ratio', ...
          'output', 'debt', 'revenue', 'capital', 'labor', 'labor_efficient', 'assets', 'assets_per_pop'};
f = @( d, s, m ) table(yearsv', m.caprates', m.govrates', m.totrates', m.qtobin', m.wages', m.rhos', ...
                      (d.outs ./ s.outs)', (d.debts ./ s.debts)', (d.revs ./ s.revs)', ...
                      (d.caps ./ s.caps)', (d.labs ./ s.labs)', (d.labeffs ./ s.labeffs)', ...
                      (d.assets ./ s.assets)', ((d.assets./d.pops) ./ (s.assets./s.pops))', ...
                      'VariableNames', header);
    
TPC_closed  = f(d_closed, s_closed, m_closed);
TPC_open    = f(d_open  , s_open  , m_open  );

if (params.depreciation == 0.056)
    filename = 'TPC_report_low_depreciation.xlsx';
else
    filename = 'TPC_report_high_depreciation.xlsx';
end
writetable(TPC_open  ,filename,'Sheet','open')
writetable(TPC_closed,filename,'Sheet','closed')
                 

%% FUNCTIONS

function x = build_var(x_name, dynamic_struct, market_struct, scenario, static)

    if strcmp(x_name, 'tax')
        
        x = dynamic_struct.pits + dynamic_struct.cits_domestic + dynamic_struct.ssts;
        
    elseif strcmp(x_name, 'totinc')
        
        % Define time constants
        s = ParamGenerator.timing( scenario );
        T_life = s.T_life; T_model = s.T_model;
        % Define assets grid and dimensions
        s = ParamGenerator.grids( scenario );
        ndem = s.ndem; nz = s.nz; nk = s.nk; nb = s.nb;
        karray = repmat(reshape(s.kv, [1,nk,1,1,1,1]),[nz,1,nb,T_life,T_model+1,ndem]);
        DIST = zeros(nz,nk,nb,T_life,T_model+1,ndem);
        s = load(fullfile(PathFinder.getWorkingDir(scenario.currentPolicy.steady), 'distribution.mat'), 'DIST');
        DIST(:,:,:,:,1,:) = reshape(sum(s.DIST, 5), [nz,nk,nb,T_life,1,ndem]);

        if static
            s = load(fullfile(PathFinder.getWorkingDir(scenario), 'Static_distribution.mat'), 'Static_DIST');
            DIST(:,:,:,:,2:end,:) = reshape(sum(s.Static_DIST, 5), [nz,nk,nb,T_life,T_model,ndem]);
        else
            s = load(fullfile(PathFinder.getWorkingDir(scenario), 'distribution.mat'), 'DIST');
            DIST(:,:,:,:,2:end,:) = reshape(sum(s.DIST, 5), [nz,nk,nb,T_life,T_model,ndem]);
        end
        
        f = @(F) sum(sum(reshape(DIST .* F, [], T_model+1, ndem), 1), 3);
%         f = @(F) squeeze(sum(sum(sum(sum(sum(F,1),2),3),4),6));
        x = dynamic_struct.labincs + market_struct.totrates .* f(karray);
        
    elseif strcmp(x_name, 'totincwss')
        
        x = dynamic_struct.totinc + dynamic_struct.bens;
        
    else
        
        error('Variable does not exist.')
        
    end

end
