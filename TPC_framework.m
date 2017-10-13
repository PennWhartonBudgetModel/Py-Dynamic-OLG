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

d_open   = load(fullfile(open_dir  , 'dynamics.mat'), 'debts', 'revs', 'outs', 'caps', 'labs');
d_closed = load(fullfile(closed_dir, 'dynamics.mat'), 'debts', 'revs', 'outs', 'caps', 'labs');

s_open   = load(fullfile(open_dir  , 'statics.mat' ), 'debts', 'revs', 'outs', 'caps', 'labs');
s_closed = load(fullfile(closed_dir, 'statics.mat' ), 'debts', 'revs', 'outs', 'caps', 'labs');

m_open   = load(fullfile(open_dir  , 'market.mat' ), 'caprates', 'wages');
m_closed = load(fullfile(closed_dir, 'market.mat' ), 'caprates', 'wages');

Res_open    = [(d_open.outs ./ s_open.outs)     ; (d_open.debts ./ s_open.debts);   ...
               (d_open.revs ./ s_open.revs)     ; (d_open.caps ./ s_open.caps)  ;   ...
               (d_open.labs ./ s_open.labs)     ;                                   ...
               m_open.caprates   ;  m_open.wages;   ];
               
Res_closed  = [(d_closed.outs ./ s_closed.outs) ; (d_closed.debts ./ s_closed.debts);   ...
               (d_closed.revs ./ s_closed.revs) ; (d_closed.caps ./ s_closed.caps)  ;   ...
               (d_closed.labs ./ s_closed.labs) ;                                       ...
               m_closed.caprates   ;  m_closed.wages;   ];
           
Res_both    = Res_closed * 0.6 + Res_open * 0.4;

Res_closed  = Res_closed';
Res_open    = Res_open';
Res_both    = Res_both';

open('Res_closed');
open('Res_open');
open('Res_both');

                        