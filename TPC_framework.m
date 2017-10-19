%%
% Report generator for dynamic model parameters and outputs.
%
%%

wipe;
params = ParamGenerator.invert(struct('savelas', 0.6, 'labelas', 0.6));

% Solve for baseline steady state
scenario   = Scenario(struct('economy', 'closed', 'beta', params.beta, 'gamma', params.gamma, ...
                             'sigma', params.sigma, 'modelunit_dollar', params.modelunit_dollar, ...
                             'bequest_phi_1', 0, 'base_brackets', 'Framework'));
sc_steady  = scenario.currentPolicy.steady;
steady_dir = PathFinder.getWorkingDir(sc_steady);

% Solve for counterfactuals
sc_open    = scenario.open();
sc_closed  = scenario.closed();
open_dir   = PathFinder.getWorkingDir(sc_open);
closed_dir = PathFinder.getWorkingDir(sc_closed);

ModelSolver.solve(scenario);

% Import variables in dynamics file
d_steady = load(fullfile(steady_dir  , 'dynamics.mat'), ... 
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'assets', 'pops', 'cons' );
d_open   = load(fullfile(open_dir  , 'dynamics.mat'), ... 
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'assets', 'pops', 'cons' );
d_closed = load(fullfile(closed_dir, 'dynamics.mat'), ...
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'assets', 'pops', 'cons' );

% Import variables in statics file
s_open   = load(fullfile(open_dir  , 'statics.mat' ), ...
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'assets', 'pops', 'cons');
s_closed = load(fullfile(closed_dir, 'statics.mat' ), ...
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'assets', 'pops', 'cons');

% Import variables in market file
m_steady = load(fullfile(steady_dir, 'market.mat' ), 'caprates', 'wages', 'rhos', 'kpricescale', 'govrates', 'totrates', 'qtobin');
m_open   = load(fullfile(open_dir  , 'market.mat' ), 'caprates', 'wages', 'rhos', 'kpricescale', 'govrates', 'totrates', 'qtobin');
m_closed = load(fullfile(closed_dir, 'market.mat' ), 'caprates', 'wages', 'rhos', 'kpricescale', 'govrates', 'totrates', 'qtobin');

% Include steady state year in all series
include_ss = @(x, xss) [xss x];
os = {'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'assets', 'cons', 'pops'};
for o = os
    d_open.(o{1})   = include_ss(d_open.(o{1})  , d_steady.(o{1}));
    s_open.(o{1})   = include_ss(s_open.(o{1})  , d_steady.(o{1}));
    d_closed.(o{1}) = include_ss(d_closed.(o{1}), d_steady.(o{1}));
    s_closed.(o{1}) = include_ss(s_closed.(o{1}), d_steady.(o{1}));
end
os = {'caprates', 'wages', 'rhos', 'kpricescale', 'govrates', 'totrates', 'qtobin'};
for o = os
    m_open.(o{1})   = include_ss(m_open.(o{1})  , m_steady.(o{1}));
    m_closed.(o{1}) = include_ss(m_closed.(o{1}), m_steady.(o{1}));
end

% Years vector
s = ParamGenerator.timing(sc_closed);
yearsv = [(s.first_transition_year - 1):1:(s.first_transition_year -1 + s.T_model)];

% Raw values reports
header      = {'year', 'caprates', 'govrates', 'totrates', 'wages', 'rhos', ...
               'd_outs', 'd_caps', 'd_labs', 'd_labeffs', 'd_labeffs_by_pops', 'd_cons', 'd_assets', 'd_debts', ...
               's_outs', 's_caps', 's_labs', 's_labeffs', 's_labeffs_by_pops', 's_cons', 's_assets', 's_debts'};

f = @( d, s, m ) table(yearsv', m.caprates', m.govrates', m.totrates', m.wages', m.rhos', ...
                  d.outs', d.caps', d.labs', d.labeffs', (d.labeffs ./ d.pops)' , d.cons', d.assets', d.debts', ...
                  s.outs', s.caps', s.labs', s.labeffs', (s.labeffs ./ s.pops)' , s.cons', s.assets', s.debts',...
                  'VariableNames', header);
           
res_open   = f(d_open, s_open, m_open)      ;
res_closed = f(d_closed, s_closed, m_closed);

writetable(res_open, 'res_open.csv')
           
% Deltas report
header = {'year', 'outs', 'debts', 'revs', 'caps', 'labs', 'labeffs', 'assets'};
f = @( d, s ) table(yearsv', (d.outs ./ s.outs)', (d.debts ./ s.debts)', (d.revs ./ s.revs)',   ...
         (d.caps ./ s.caps)', (d.labs ./ s.labs)', (d.labeffs ./ s.labeffs)'          ,   ...
         (d.assets ./ s.assets)', 'VariableNames', header);
    
delta_closed  = f(d_closed, s_closed);
delta_open    = f(d_open  , s_open  );

open('delta_closed');
open('delta_open');

                        