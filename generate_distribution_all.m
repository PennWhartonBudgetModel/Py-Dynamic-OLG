%%
% Generate population distribution for next year from population distribution for current year.
% 
%%


function [DIST_next] = generate_distribution_all(DIST_year, DIST_grow, K, B, nz, nk, nb, T_life, ng, transz, kv, bv, surv)

DIST_next = DIST_grow;

parfor age = 2:T_life

    DIST_year_age = DIST_year(:,:,:,age-1,:);
    DIST_grow_age = DIST_grow(:,:,:,age  ,:);

    K_age = K(:,:,:,age-1);
    B_age = B(:,:,:,age-1);

    transz_age = transz(:,:,age);
    surv_age   = surv(age-1);

    DIST_next_age = generate_distribution(DIST_year_age, DIST_grow_age, K_age, B_age, nz, nk, nb, T_life, ng, transz_age, kv, bv, surv_age);

    DIST_next(:,:,:,age,:) = DIST_next_age;
    
end

end