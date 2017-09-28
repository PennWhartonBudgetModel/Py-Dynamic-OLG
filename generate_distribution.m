%%
% Generate population distribution for next year from population distribution for current year.
% 
%%


function [DIST_next] = generate_distribution(DIST_year, DIST_grow, K, B, nz, nk, nb, T_life, ng, transz, kv, bv, surv) %#codegen


%% Argument verification

nz_max =  50;
nk_max =  50;
nb_max =  50;
T_max  = 100;
ng_max =  10;

assert( isa(DIST_year   , 'double'  ) && (size(DIST_year    , 1) <= nz_max  ) && (size(DIST_year    , 2) <= nk_max  ) && (size(DIST_year    , 3) <= nb_max  ) && (size(DIST_year    , 4) <= T_max   ) && (size(DIST_year    , 5) <= ng_max  ) );
assert( isa(DIST_grow   , 'double'  ) && (size(DIST_grow    , 1) <= nz_max  ) && (size(DIST_grow    , 2) <= nk_max  ) && (size(DIST_grow    , 3) <= nb_max  ) && (size(DIST_grow    , 4) <= T_max   ) && (size(DIST_grow    , 5) <= ng_max  ) );
assert( isa(K           , 'double'  ) && (size(K            , 1) <= nz_max  ) && (size(K            , 2) <= nk_max  ) && (size(K            , 3) <= nb_max  ) && (size(K            , 4) <= T_max   ) );
assert( isa(B           , 'double'  ) && (size(B            , 1) <= nz_max  ) && (size(B            , 2) <= nk_max  ) && (size(B            , 3) <= nb_max  ) && (size(B            , 4) <= T_max   ) );

assert( isa(nz          , 'double'  ) && (size(nz           , 1) == 1       ) && (size(nz           , 2) == 1       ) );
assert( isa(nk          , 'double'  ) && (size(nk           , 1) == 1       ) && (size(nk           , 2) == 1       ) );
assert( isa(nb          , 'double'  ) && (size(nb           , 1) == 1       ) && (size(nb           , 2) == 1       ) );
assert( isa(T_life      , 'double'  ) && (size(T_life       , 1) == 1       ) && (size(T_life       , 2) == 1       ) );
assert( isa(ng          , 'double'  ) && (size(ng           , 1) == 1       ) && (size(ng           , 2) == 1       ) );
assert( isa(transz      , 'double'  ) && (size(transz       , 1) <= nz_max  ) && (size(transz       , 2) <= nz_max  ) );
assert( isa(kv          , 'double'  ) && (size(kv           , 1) <= nk_max  ) && (size(kv           , 2) == 1       ) );
assert( isa(bv          , 'double'  ) && (size(bv           , 1) <= nb_max  ) && (size(bv           , 2) == 1       ) );
assert( isa(surv        , 'double'  ) && (size(surv         , 1) == 1       ) && (size(surv         , 2) <= T_max   ) );



%% Distribution generation

% Initialize population distribution for next year using population growth distribution
DIST_next = DIST_grow;

for age = 2:T_life
    
    % Extract optimal decision values
    k_t = K(:,:,:,age-1);
    b_t = B(:,:,:,age-1);
    
    % Find indices of nearest values in decision value discretization vectors
    jk_lt = ones(size(k_t));
    for elem = 1:length(k_t(:))
        jk_lt(elem) = find(kv(1:end-1) <= k_t(elem), 1, 'last');
    end
    jk_gt = jk_lt + 1;
    
    jb_lt = ones(size(b_t));
    for elem = 1:length(b_t(:))
        jb_lt(elem) = find(bv(1:end-1) <= b_t(elem), 1, 'last');
    end
    jb_gt = jb_lt + 1;
    
    % Calculate linear weights for nearest values
    wk_lt = (kv(jk_gt) - k_t) ./ (kv(jk_gt) - kv(jk_lt));
    wk_gt = 1 - wk_lt;
    
    wb_lt = (bv(jb_gt) - b_t) ./ (bv(jb_gt) - bv(jb_lt));
    wb_gt = 1 - wb_lt;
    
    % Checks -> only work in the absence of mex file!
    assert( all(wk_lt(:)>=0) && all(wk_gt(:)>=0) && all(wb_lt(:)>=0) && all(wb_gt(:)>=0), 'Negative weights to compute households distribution.')       
    
    for jz = 1:nz
        
        % Apply survival and productivity transformations to population distribution for current year
        DIST_transz = DIST_year(:,:,:,age-1,:) * surv(age-1) .* repmat(reshape(transz(:,jz), [nz,1,1,1,1]), [1,nk,nb,1,ng]);
        assert(all(DIST_transz(:)>=0), 'Negative mass of people at DIST_transz.')
        
        % Redistribute population distribution from current year to next year according to target indices and weights
        for ib = 1:nb, for ik = 1:nk, for iz = 1:nz %#ok<ALIGN>
            DIST_next(jz, jk_lt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) = DIST_next(jz, jk_lt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) + wk_lt(iz,ik,ib)*wb_lt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
            DIST_next(jz, jk_gt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) = DIST_next(jz, jk_gt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) + wk_gt(iz,ik,ib)*wb_lt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
            DIST_next(jz, jk_lt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) = DIST_next(jz, jk_lt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) + wk_lt(iz,ik,ib)*wb_gt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
            DIST_next(jz, jk_gt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) = DIST_next(jz, jk_gt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) + wk_gt(iz,ik,ib)*wb_gt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
        end, end, end


    end
    
end


end