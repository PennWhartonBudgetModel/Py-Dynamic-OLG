wipe;
do_solve = false;

params = ParamGenerator.invert(struct('savelas', 0.5, 'labelas', 0.5));

tpc_framework = struct('base_brackets', 'Framework', ...
                        'beta', params.beta, 'gamma', params.gamma, 'sigma', params.sigma, ...
                            'modelunit_dollar', params.modelunit_dollar, ... 
                            'bequest_phi_1', 0, 'economy', 'closed' );
                        
scenario = Scenario(tpc_framework);

if (do_solve )
    ModelSolver.solve(scenario);
end

sc_open     = scenario.open();
sc_closed   = scenario.closed();

open_dir    = PathFinder.getWorkingDir(sc_open);
closed_dir  = PathFinder.getWorkingDir(sc_closed);

d_open   = load(fullfile(open_dir  , 'dynamics.mat'), 'debts', 'revs', 'outs');
d_closed = load(fullfile(closed_dir, 'dynamics.mat'), 'debts', 'revs', 'outs');

s_open   = load(fullfile(open_dir  , 'statics.mat' ), 'debts', 'revs', 'outs');
s_closed = load(fullfile(closed_dir, 'statics.mat' ), 'debts', 'revs', 'outs');

m_open   = load(fullfile(open_dir  , 'market.mat' ), 'caprates');
m_closed = load(fullfile(closed_dir, 'market.mat' ), 'caprates');

Res_open    = [(d_open.outs ./ s_open.outs ) ; (d_open.debts ./ s_open.debts); ...
               (d_open.revs ./ s_open.revs)      ; m_open.caprates   ];
Res_closed  = [(d_closed.outs ./ s_closed.outs ) ; (d_closed.debts ./ s_closed.debts); ...
               (d_closed.revs ./ s_closed.revs)  ; m_closed.caprates ];
Res_both    = Res_closed * 0.6 + Res_open * 0.4;

Res_closed  = Res_closed';
Res_open    = Res_open';
Res_both    = Res_both';

open('Res_closed');
open('Res_open');
open('Res_both');

                        