%%
% Generate population distribution.
% 
%%


function [DIST] = generate_distribution(DIST_last, DIST_new, K, B, nz, nk, nb, T_life, transz, ks, bs) %#codegen


%% Argument verification

T_max  = 100;
nz_max =  50;
nk_max =  50;
nb_max =  50;

assert( isa(DIST_last   , 'double'  ) && (size(DIST_last    , 1) <= nz_max  ) && (size(DIST_last    , 2) <= nk_max  ) && (size(DIST_last    , 3) <= nb_max  ) && (size(DIST_last    , 4) <= T_max   ) );
assert( isa(DIST_new    , 'double'  ) && (size(DIST_new     , 1) <= nz_max  ) && (size(DIST_new     , 2) <= nk_max  ) && (size(DIST_new     , 3) <= nb_max  ) && (size(DIST_new     , 4) <= T_max   ) );
assert( isa(K           , 'double'  ) && (size(K            , 1) <= nz_max  ) && (size(K            , 2) <= nk_max  ) && (size(K            , 3) <= nb_max  ) && (size(K            , 4) <= T_max   ) );
assert( isa(B           , 'double'  ) && (size(B            , 1) <= nz_max  ) && (size(B            , 2) <= nk_max  ) && (size(B            , 3) <= nb_max  ) && (size(B            , 4) <= T_max   ) );

assert( isa(nz          , 'double'  ) && (size(nz           , 1) == 1       ) && (size(nz           , 2) == 1       ) );
assert( isa(nk          , 'double'  ) && (size(nk           , 1) == 1       ) && (size(nk           , 2) == 1       ) );
assert( isa(nb          , 'double'  ) && (size(nb           , 1) == 1       ) && (size(nb           , 2) == 1       ) );
assert( isa(T_life      , 'double'  ) && (size(T_life       , 1) == 1       ) && (size(T_life       , 2) == 1       ) );
assert( isa(transz      , 'double'  ) && (size(transz       , 1) <= nz_max  ) && (size(transz       , 2) <= nz_max  ) );
assert( isa(ks          , 'double'  ) && (size(ks           , 1) <= nk_max  ) && (size(ks           , 2) == 1       ) );
assert( isa(bs          , 'double'  ) && (size(bs           , 1) <= nb_max  ) && (size(bs           , 2) == 1       ) );



%% Distribution generation

DIST = DIST_new;

for age = 2:T_life

    % Extract optimal k and b decision values
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

        % Apply survival and productivity transformations to cohort distribution from current year
        DIST_transz = DIST_last(:,:,:,age-1) .* repmat(reshape(transz(:,jz), [nz,1,1]), [1,nk,nb]);
        % DIST_transz = DIST__(:,:,:,age-1) * surv(age-1) .* repmat(reshape(transz(:,jz), [nz,1,1]), [1,nk,nb]);

        % Redistribute cohort for next year according to target indices and weights
        for elem = 1:numel(DIST_transz)
            DIST(jz, jk_lt(elem), jb_lt(elem), age) = DIST(jz, jk_lt(elem), jb_lt(elem), age) + wk_lt(elem)*wb_lt(elem)*DIST_transz(elem);
            DIST(jz, jk_gt(elem), jb_lt(elem), age) = DIST(jz, jk_gt(elem), jb_lt(elem), age) + wk_gt(elem)*wb_lt(elem)*DIST_transz(elem);
            DIST(jz, jk_lt(elem), jb_gt(elem), age) = DIST(jz, jk_lt(elem), jb_gt(elem), age) + wk_lt(elem)*wb_gt(elem)*DIST_transz(elem);
            DIST(jz, jk_gt(elem), jb_gt(elem), age) = DIST(jz, jk_gt(elem), jb_gt(elem), age) + wk_gt(elem)*wb_gt(elem)*DIST_transz(elem);
        end

    end

end


end