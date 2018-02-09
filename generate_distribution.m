%%
% Generate population distribution for a given age for next year from population distribution for current year.
% 
%%


function [DIST_next_age] = generate_distribution(DIST_year_age, DIST_grow_age, K_age, B_age, nz, nk, nb, T_life, ng, transz_age, kv, bv, surv_age) %#codegen


%% Argument verification

nz_max =  50;
nk_max =  50;
nb_max =  50;
ng_max =  10;

assert( isa(DIST_year_age   , 'double'  ) && (size(DIST_year_age    , 1) <= nz_max  ) && (size(DIST_year_age    , 2) <= nk_max  ) && (size(DIST_year_age    , 3) <= nb_max  ) && (size(DIST_year_age    , 4) == 1   ) && (size(DIST_year_age    , 5) <= ng_max  ) );
assert( isa(DIST_grow_age   , 'double'  ) && (size(DIST_grow_age    , 1) <= nz_max  ) && (size(DIST_grow_age    , 2) <= nk_max  ) && (size(DIST_grow_age    , 3) <= nb_max  ) && (size(DIST_grow_age    , 4) == 1   ) && (size(DIST_grow_age    , 5) <= ng_max  ) );
assert( isa(K_age           , 'double'  ) && (size(K_age            , 1) <= nz_max  ) && (size(K_age            , 2) <= nk_max  ) && (size(K_age            , 3) <= nb_max  ) );
assert( isa(B_age           , 'double'  ) && (size(B_age            , 1) <= nz_max  ) && (size(B_age            , 2) <= nk_max  ) && (size(B_age            , 3) <= nb_max  ) );

assert( isa(nz              , 'double'  ) && (size(nz               , 1) == 1       ) && (size(nz               , 2) == 1       ) );
assert( isa(nk              , 'double'  ) && (size(nk               , 1) == 1       ) && (size(nk               , 2) == 1       ) );
assert( isa(nb              , 'double'  ) && (size(nb               , 1) == 1       ) && (size(nb               , 2) == 1       ) );
assert( isa(T_life          , 'double'  ) && (size(T_life           , 1) == 1       ) && (size(T_life           , 2) == 1       ) );
assert( isa(ng              , 'double'  ) && (size(ng               , 1) == 1       ) && (size(ng               , 2) == 1       ) );
assert( isa(transz_age      , 'double'  ) && (size(transz_age       , 1) <= nz_max  ) && (size(transz_age       , 2) <= nz_max  ) );
assert( isa(kv              , 'double'  ) && (size(kv               , 1) <= nk_max  ) && (size(kv               , 2) == 1       ) );
assert( isa(bv              , 'double'  ) && (size(bv               , 1) <= nb_max  ) && (size(bv               , 2) == 1       ) );
assert( isa(surv_age        , 'double'  ) && (size(surv_age         , 1) == 1       ) && (size(surv_age         , 2) == 1       ) );



%% Distribution generation

% Initialize population distribution for next year using population growth distribution
DIST_next_age = DIST_grow_age;

% Extract optimal decision values
k_t = K_age;
b_t = B_age;

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
    DIST_transz = DIST_year_age(:,:,:,1,:) * surv_age .* repmat(reshape(transz_age(:,jz), [nz,1,1,1,1]), [1,nk,nb,1,ng]);
    assert(all(DIST_transz(:)>=0), 'Negative mass of people at DIST_transz.')

    % Redistribute population distribution from current year to next year according to target indices and weights
    for ib = 1:nb, for ik = 1:nk, for iz = 1:nz %#ok<ALIGN>
        DIST_next_age(jz, jk_lt(iz,ik,ib), jb_lt(iz,ik,ib), 1, :) = DIST_next_age(jz, jk_lt(iz,ik,ib), jb_lt(iz,ik,ib), 1, :) + wk_lt(iz,ik,ib)*wb_lt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
        DIST_next_age(jz, jk_gt(iz,ik,ib), jb_lt(iz,ik,ib), 1, :) = DIST_next_age(jz, jk_gt(iz,ik,ib), jb_lt(iz,ik,ib), 1, :) + wk_gt(iz,ik,ib)*wb_lt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
        DIST_next_age(jz, jk_lt(iz,ik,ib), jb_gt(iz,ik,ib), 1, :) = DIST_next_age(jz, jk_lt(iz,ik,ib), jb_gt(iz,ik,ib), 1, :) + wk_lt(iz,ik,ib)*wb_gt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
        DIST_next_age(jz, jk_gt(iz,ik,ib), jb_gt(iz,ik,ib), 1, :) = DIST_next_age(jz, jk_gt(iz,ik,ib), jb_gt(iz,ik,ib), 1, :) + wk_gt(iz,ik,ib)*wb_gt(iz,ik,ib)*DIST_transz(iz,ik,ib,1,:);
    end, end, end

end

end
