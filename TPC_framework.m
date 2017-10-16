wipe;
do_solve = true;

params = ParamGenerator.invert(struct('savelas', 0.5, 'labelas', 0.5));

tpc_framework = struct('base_brackets', 'TEST', ...
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

d_open   = load(fullfile(open_dir  , 'dynamics.mat'), ... 
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'assets', 'pops' );
d_closed = load(fullfile(closed_dir, 'dynamics.mat'), ...
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'assets', 'pops' );

s_open   = load(fullfile(open_dir  , 'statics.mat' ), ...
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'assets');
s_closed = load(fullfile(closed_dir, 'statics.mat' ), ...
            'debts', 'revs', 'outs', 'caps', 'labs', 'labeffs', 'assets');

m_open   = load(fullfile(open_dir  , 'market.mat' ), 'caprates', 'wages', 'rhos', 'kpricescale');
m_closed = load(fullfile(closed_dir, 'market.mat' ), 'caprates', 'wages', 'rhos', 'kpricescale');

f = @( d, s, m ) ...
        [(d.outs ./ s.outs)     ; (d.debts ./ s.debts); (d.revs ./ s.revs)  ;       ...
         (d.caps ./ s.caps)     ; (d.labs ./ s.labs)  ; (d.labeffs ./ s.labeffs);   ...
         d.assets ./ s.assets   ;                                                   ...
         m.caprates     ; m.wages       ; m.rhos      ;                             ...
         d.caps         ; d.labs        ; d.assets    ;                             ...
         d.debts        ; d.labeffs     ; d.labeffs ./ d.pops ;                     ...
        ];
    
Res_closed  = f(d_closed, s_closed, m_closed);
Res_open    = f(d_open, s_open, m_open);

Res_closed  = Res_closed';
Res_open    = Res_open';

open('Res_closed');
open('Res_open');

                        