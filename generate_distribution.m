%%
% Generate population distribution for a year.
% 
%%


function [DIST_year] = generate_distribution(DIST_last, K, B, ...
             nz, nk, nb, T_life, ng, transz, ks, bs, surv, g, ...
             pgr, legal_rate, illegal_rate, imm_age, DISTz_age, amnesty, deportation) %#codegen


%% Argument verification

nz_max =  50;
nk_max =  50;
nb_max =  50;
T_max  = 100;
ng_max =  10;

assert( isa(DIST_last   , 'double'  ) && (size(DIST_last    , 1) <= nz_max  ) && (size(DIST_last    , 2) <= nk_max  ) && (size(DIST_last    , 3) <= nb_max  ) && (size(DIST_last    , 4) <= T_max   ) && (size(DIST_last    , 5) <= ng_max  ) );
assert( isa(K           , 'double'  ) && (size(K            , 1) <= nz_max  ) && (size(K            , 2) <= nk_max  ) && (size(K            , 3) <= nb_max  ) && (size(K            , 4) <= T_max   ) );
assert( isa(B           , 'double'  ) && (size(B            , 1) <= nz_max  ) && (size(B            , 2) <= nk_max  ) && (size(B            , 3) <= nb_max  ) && (size(B            , 4) <= T_max   ) );

assert( isa(nz          , 'double'  ) && (size(nz           , 1) == 1       ) && (size(nz           , 2) == 1       ) );
assert( isa(nk          , 'double'  ) && (size(nk           , 1) == 1       ) && (size(nk           , 2) == 1       ) );
assert( isa(nb          , 'double'  ) && (size(nb           , 1) == 1       ) && (size(nb           , 2) == 1       ) );
assert( isa(T_life      , 'double'  ) && (size(T_life       , 1) == 1       ) && (size(T_life       , 2) == 1       ) );
assert( isa(ng          , 'double'  ) && (size(ng           , 1) == 1       ) && (size(ng           , 2) == 1       ) );
assert( isa(transz      , 'double'  ) && (size(transz       , 1) <= nz_max  ) && (size(transz       , 2) <= nz_max  ) );
assert( isa(ks          , 'double'  ) && (size(ks           , 1) <= nk_max  ) && (size(ks           , 2) == 1       ) );
assert( isa(bs          , 'double'  ) && (size(bs           , 1) <= nb_max  ) && (size(bs           , 2) == 1       ) );
assert( isa(surv        , 'double'  ) && (size(surv         , 1) == 1       ) && (size(surv         , 2) <= T_max   ) );
assert( isstruct(g)                   && (size(g            , 1) == 1       ) && (size(g            , 2) == 1       ) && ...
        isa(g.citizen   , 'double'  ) && (size(g.citizen    , 1) == 1       ) && (size(g.citizen    , 2) == 1       ) && ...
        isa(g.legal     , 'double'  ) && (size(g.legal      , 1) == 1       ) && (size(g.legal      , 2) == 1       ) && ...
        isa(g.illegal   , 'double'  ) && (size(g.illegal    , 1) == 1       ) && (size(g.illegal    , 2) == 1       ) );

assert( isa(pgr         , 'double'  ) && (size(pgr          , 1) == 1       ) && (size(pgr          , 2) == 1       ) );
assert( isa(legal_rate  , 'double'  ) && (size(legal_rate   , 1) == 1       ) && (size(legal_rate   , 2) == 1       ) );
assert( isa(illegal_rate, 'double'  ) && (size(illegal_rate , 1) == 1       ) && (size(illegal_rate , 2) == 1       ) );
assert( isa(imm_age     , 'double'  ) && (size(imm_age      , 1) == 1       ) && (size(imm_age      , 2) <= T_max   ) );
assert( isa(DISTz_age   , 'double'  ) && (size(DISTz_age    , 1) <= nz_max  ) && (size(DISTz_age    , 2) <= T_max   ) && (size(DISTz_age    , 3) <= ng_max  ) );
assert( isa(amnesty     , 'double'  ) && (size(amnesty      , 1) == 1       ) && (size(amnesty      , 2) == 1       ) );
assert( isa(deportation , 'double'  ) && (size(deportation  , 1) == 1       ) && (size(deportation  , 2) == 1       ) );



%% Distribution generation

% Initialize population distribution array
DIST_year = zeros(nz,nk,nb,T_life,ng);

% Add population growth
% (Intermediary arrays necessarily for code generation)
P = sum(DIST_last(:));
D = reshape(DISTz_age(:,1,g.citizen,1,1), [nz,1,1,1     ,1]) * P * pgr;                                                                      DIST_year(:,1,1,1,g.citizen) = D(1:nz,1,1,1       ,1);
D = reshape(DISTz_age(:,:,g.legal  ,1,1), [nz,1,1,T_life,1]) * P * legal_rate   .* repmat(reshape(imm_age, [1,1,1,T_life,1]), [nz,1,1,1,1]); DIST_year(:,1,1,:,g.legal  ) = D(1:nz,1,1,1:T_life,1);
D = reshape(DISTz_age(:,:,g.illegal,1,1), [nz,1,1,T_life,1]) * P * illegal_rate .* repmat(reshape(imm_age, [1,1,1,T_life,1]), [nz,1,1,1,1]); DIST_year(:,1,1,:,g.illegal) = D(1:nz,1,1,1:T_life,1);


for age = 2:T_life

    % Extract optimal decision values
    k_t = K(:,:,:,age-1);
    b_t = B(:,:,:,age-1);
    
    % Find indices of nearest values in ks and bs series
    jk_lt = ones(size(k_t));
    for elem = 1:length(k_t(:))
        jk_lt(elem) = find(ks(1:end-1) <= k_t(elem), 1, 'last');
    end
    jk_gt = jk_lt + 1;
    
    jb_lt = ones(size(b_t));
    for elem = 1:length(b_t(:))
        jb_lt(elem) = find(bs(1:end-1) <= b_t(elem), 1, 'last');
    end
    jb_gt = jb_lt + 1;
    
    % Calculate linear weights for nearest values
    wk_lt = (ks(jk_gt) - k_t) ./ (ks(jk_gt) - ks(jk_lt));
    wk_gt = 1 - wk_lt;
    
    wb_lt = (bs(jb_gt) - b_t) ./ (bs(jb_gt) - bs(jb_lt));
    wb_gt = 1 - wb_lt;
    
    for jz = 1:nz
        
        % Apply survival and productivity transformations to last year's distribution
        DIST_transz = DIST_last(:,:,:,age-1,:) * surv(age-1) .* repmat(reshape(transz(:,jz), [nz,1,1,1,1]), [1,nk,nb,1,ng]);
        
        % Redistribute last year's distribution according to target indices and weights
        for ib = 1:nb, for ik = 1:nk, for iz = 1:nz %#ok<ALIGN>
            DIST_year(jz, jk_lt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) = DIST_year(jz, jk_lt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) + wk_lt(iz,ik,ib)*wb_lt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
            DIST_year(jz, jk_gt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) = DIST_year(jz, jk_gt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) + wk_gt(iz,ik,ib)*wb_lt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
            DIST_year(jz, jk_lt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) = DIST_year(jz, jk_lt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) + wk_lt(iz,ik,ib)*wb_gt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
            DIST_year(jz, jk_gt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) = DIST_year(jz, jk_gt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) + wk_gt(iz,ik,ib)*wb_gt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
        end, end, end
        
    end
    
end


% Increase legal immigrant population for amnesty, maintaining productivity distributions
% (Intermediary array necessarily for code generation)
D = repmat(sum(DIST_year(:,:,:,:,g.legal), 1), [nz,1,1,1,1]);
DISTz_legal = DIST_year(:,:,:,:,g.legal) ./ D(1:nz,1:nk,1:nb,1:T_life,1);
DISTz_legal(isnan(DISTz_legal)) = 1/nz;

DIST_year(:,:,:,:,g.legal) = DIST_year(:,:,:,:,g.legal) + repmat(sum(amnesty*DIST_year(:,:,:,:,g.illegal), 1), [nz,1,1,1,1]).*DISTz_legal;

% Reduce illegal immigrant population for amnesty and deportation
DIST_year(:,:,:,:,g.illegal) = (1-amnesty-deportation)*DIST_year(:,:,:,:,g.illegal);


end