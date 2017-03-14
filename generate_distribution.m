%%
% Generate population distribution for a year.
% 
%%


function [DIST_year] = generate_distribution(DIST_last, DIST_grow, K, B, nz, nk, nb, T_life, ngroups, transz, ks, bs, surv) %#codegen


%% Argument verification

T_max       = 100;
nz_max      =  50;
nk_max      =  50;
nb_max      =  50;
ngroups_max =  10;

assert( isa(DIST_last   , 'double'  ) && (size(DIST_last    , 1) <= nz_max  ) && (size(DIST_last    , 2) <= nk_max  ) && (size(DIST_last    , 3) <= nb_max  ) && (size(DIST_last    , 4) <= T_max   ) && (size(DIST_last    , 5) <= ngroups_max ) );
assert( isa(DIST_grow   , 'double'  ) && (size(DIST_grow    , 1) <= nz_max  ) && (size(DIST_grow    , 2) <= nk_max  ) && (size(DIST_grow    , 3) <= nb_max  ) && (size(DIST_grow    , 4) <= T_max   ) && (size(DIST_grow    , 5) <= ngroups_max ) );
assert( isa(K           , 'double'  ) && (size(K            , 1) <= nz_max  ) && (size(K            , 2) <= nk_max  ) && (size(K            , 3) <= nb_max  ) && (size(K            , 4) <= T_max   ) );
assert( isa(B           , 'double'  ) && (size(B            , 1) <= nz_max  ) && (size(B            , 2) <= nk_max  ) && (size(B            , 3) <= nb_max  ) && (size(B            , 4) <= T_max   ) );

assert( isa(nz          , 'double'  ) && (size(nz           , 1) == 1       ) && (size(nz           , 2) == 1       ) );
assert( isa(nk          , 'double'  ) && (size(nk           , 1) == 1       ) && (size(nk           , 2) == 1       ) );
assert( isa(nb          , 'double'  ) && (size(nb           , 1) == 1       ) && (size(nb           , 2) == 1       ) );
assert( isa(T_life      , 'double'  ) && (size(T_life       , 1) == 1       ) && (size(T_life       , 2) == 1       ) );
assert( isa(ngroups     , 'double'  ) && (size(ngroups      , 1) == 1       ) && (size(ngroups      , 2) == 1       ) );
assert( isa(transz      , 'double'  ) && (size(transz       , 1) <= nz_max  ) && (size(transz       , 2) <= nz_max  ) );
assert( isa(ks          , 'double'  ) && (size(ks           , 1) <= nk_max  ) && (size(ks           , 2) == 1       ) );
assert( isa(bs          , 'double'  ) && (size(bs           , 1) <= nb_max  ) && (size(bs           , 2) == 1       ) );
assert( isa(surv        , 'double'  ) && (size(surv         , 1) == 1       ) && (size(surv         , 2) <= T_max   ) );



%% Distribution generation

DIST_year = DIST_grow;

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
        DIST_transz = DIST_last(:,:,:,age-1,:) * surv(age-1) .* repmat(reshape(transz(:,jz), [nz,1,1,1,1]), [1,nk,nb,1,ngroups]);
        
        % Redistribute last year's distribution according to target indices and weights
        for ib = 1:nb, for ik = 1:nk, for iz = 1:nz %#ok<ALIGN>
            DIST_year(jz, jk_lt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) = DIST_year(jz, jk_lt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) + wk_lt(iz,ik,ib)*wb_lt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
            DIST_year(jz, jk_gt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) = DIST_year(jz, jk_gt(iz,ik,ib), jb_lt(iz,ik,ib), age, :) + wk_gt(iz,ik,ib)*wb_lt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
            DIST_year(jz, jk_lt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) = DIST_year(jz, jk_lt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) + wk_lt(iz,ik,ib)*wb_gt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
            DIST_year(jz, jk_gt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) = DIST_year(jz, jk_gt(iz,ik,ib), jb_gt(iz,ik,ib), age, :) + wk_gt(iz,ik,ib)*wb_gt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
        end, end, end
        
    end
    
end


end