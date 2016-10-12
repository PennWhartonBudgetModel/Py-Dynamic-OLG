%%
% Distribution generator.  Generalized to support both steady state and transition path solvers.
% 
% For steady state:
% 
%   Tss       = T
%   birthyear = 0
% 
%   dist_w0 = nk x nz x nb x 1 matrix
%   dist_r0 = []
% 
%%


function [dist_w, dist_r, N_w, N_r, Kalive, Kdead, ELab, Lab, Lfpr, Fedincome, Fedit, SSrev, Fcaptax, SSexp] ...
          ...
            = generate_distributions(...
                kopt, koptss, labopt, bopt, ...
                fedincome, fedincomess, fitax, fitaxss, fsstax, fsstaxss, fcaptax, fcaptaxss, ...
                T, Tr, Tss, birthyear, ...
                kgrid, bgrid, nk, nz, nb, idem, ...
                z, tr_z, ben, ...
                dist_w0, dist_r0, ...
                mu2_idem, mu3_idem)

% Define past years, effective lifespan, and effective retirement age
% (Projections truncated at year Tss)
T_past = -min(0, birthyear);
T_eff  = min(T, Tss - birthyear);
Tr_eff = max(min(Tr, T_eff), T_past);

% Find effective numbers of working and retirement years
N_w = Tr_eff - T_past;
N_r = T_eff  - Tr_eff;

% Initialize distributions
dist_w = zeros(nk,nz,nb,N_w);
dist_r = zeros(nk,nb,N_r);

if (N_w > 0)
    dist_w(:,:,:,1) = dist_w0(:,:,:,T_past+1);
else
    dist_r(:,:,1)   = dist_r0(:,:,T_past-Tr+1);
end

% Find distributions for working years
for t = 1:N_w

    % Extract optimal k and b values
    % (Note that k and b should already bounded by the ranges of kgrid and bgrid respectively by the dynamic optimization solver)
    k_t = kopt(:,:,:,t);
    b_t = bopt(:,:,:,t);

    % Find indices of nearest values in kgrid and bgrid series
    % (Using arrayfun here leads to nontrivial anonymous function overhead)
    jk_lt = ones(size(k_t));
    for elem = 1:length(k_t(:))
        jk_lt(elem) = find(kgrid(1:end-1) <= k_t(elem), 1, 'last');
    end
    jk_gt = jk_lt + 1;

    jb_lt = ones(size(b_t));
    for elem = 1:length(b_t(:))
        jb_lt(elem) = find(bgrid(1:end-1) <= b_t(elem), 1, 'last');
    end
    jb_gt = jb_lt + 1;

    % Calculate linear weights for lower and upper nearest values
    wk_lt = (kgrid(jk_gt) - k_t) ./ (kgrid(jk_gt) - kgrid(jk_lt));
    wk_gt = 1 - wk_lt;

    wb_lt = (bgrid(jb_gt) - b_t) ./ (bgrid(jb_gt) - bgrid(jb_lt));
    wb_gt = 1 - wb_lt;

    if t == N_w
        break
    else
        for jz = 1:nz

            % Perform productivity transformation
            dist_step = repmat(tr_z(:,jz)', [nk,1,nb]) .* dist_w(:,:,:,t);

            % Calculate distributions for next time step
            for elem = 1:numel(dist_step)
                dist_w(jk_lt(elem), jz, jb_lt(elem), t+1) = dist_w(jk_lt(elem), jz, jb_lt(elem), t+1) + wb_lt(elem) * wk_lt(elem) * dist_step(elem);
                dist_w(jk_gt(elem), jz, jb_lt(elem), t+1) = dist_w(jk_gt(elem), jz, jb_lt(elem), t+1) + wb_lt(elem) * wk_gt(elem) * dist_step(elem);
                dist_w(jk_lt(elem), jz, jb_gt(elem), t+1) = dist_w(jk_lt(elem), jz, jb_gt(elem), t+1) + wb_gt(elem) * wk_lt(elem) * dist_step(elem);
                dist_w(jk_gt(elem), jz, jb_gt(elem), t+1) = dist_w(jk_gt(elem), jz, jb_gt(elem), t+1) + wb_gt(elem) * wk_gt(elem) * dist_step(elem);
            end

        end
    end

end

% Perform transition between working years to retirement years if there is at least 1 of each
if (N_w > 0) && (N_r > 0)

    dist_step = dist_w(:,:,:,N_w);

    for elem = 1:numel(dist_step)
        dist_r(jk_lt(elem), jb_lt(elem), 1) = dist_r(jk_lt(elem), jb_lt(elem), 1) + wb_lt(elem) * wk_lt(elem) * dist_step(elem);
        dist_r(jk_gt(elem), jb_lt(elem), 1) = dist_r(jk_gt(elem), jb_lt(elem), 1) + wb_lt(elem) * wk_gt(elem) * dist_step(elem);
        dist_r(jk_lt(elem), jb_gt(elem), 1) = dist_r(jk_lt(elem), jb_gt(elem), 1) + wb_gt(elem) * wk_lt(elem) * dist_step(elem);
        dist_r(jk_gt(elem), jb_gt(elem), 1) = dist_r(jk_gt(elem), jb_gt(elem), 1) + wb_gt(elem) * wk_gt(elem) * dist_step(elem);
    end

end

% Find distributions for retirement years
% (0 iterations if N_r <= 1 -- i.e. if there is only 1 year or less of retirement)
for t = 1:N_r-1

    % Extract optimal k values
    % (Note that k should already bounded by the range of kgrid by the dynamic optimization solver)
    k_t = koptss(:,:,t);

    % Find indices of nearest values in kgrid series
    jk_lt = ones(size(k_t));
    for elem = 1:length(k_t(:))
        jk_lt(elem) = find(kgrid(1:end-1) <= k_t(elem), 1, 'last');
    end
    jk_gt = jk_lt + 1;

    % Calculate linear weights for lower and upper nearest values
    wk_lt = (kgrid(jk_gt) - k_t) ./ (kgrid(jk_gt) - kgrid(jk_lt));
    wk_gt = 1 - wk_lt;

    % Calculate distributions for next time step
    for ib = 1:nb
        for ik = 1:nk
            dist_r(jk_lt(ik,ib), ib, t+1) = dist_r(jk_lt(ik,ib), ib, t+1) + wk_lt(ik,ib) * dist_r(ik,ib,t);
            dist_r(jk_gt(ik,ib), ib, t+1) = dist_r(jk_gt(ik,ib), ib, t+1) + wk_gt(ik,ib) * dist_r(ik,ib,t);
        end
    end

end

% Adjust distributions based on demographics
dist_w = bsxfun(@times, shiftdim(mu2_idem(T_past+1:Tr_eff), -2), dist_w);
dist_r = bsxfun(@times, shiftdim(mu2_idem(Tr_eff+1:T_eff),  -1), dist_r);

% Calculate aggregates
Kalive  = [sum(reshape(kopt     (:,:,:,1:N_w) .* dist_w, [], N_w), 1), ...
           sum(reshape(koptss   (:,:,1:N_r)   .* dist_r, [], N_r), 1)];

Kdead     = Kalive .* mu3_idem(T_past+1:T_eff) ./ mu2_idem(T_past+1:T_eff);

ELab      = sum(reshape(labopt(:,:,:,1:N_w) .* repmat(reshape(z(:,T_past+1:Tr_eff,idem), [1,nz,1,N_w]), [nk,1,nb,1]) .* dist_w, [], N_w), 1);

Lab       = sum(reshape(labopt(:,:,:,1:N_w) .* dist_w, [], N_w), 1);

Lfpr      = sum(reshape( (labopt(:,:,:,1:N_w) >= 0.01) .* dist_w, [], N_w), 1);

Fedincome = [sum(reshape(fedincome  (:,:,:,1:N_w) .* dist_w, [], N_w), 1), ...
             sum(reshape(fedincomess(:,:,1:N_r)   .* dist_r, [], N_r), 1)];

Fedit     = [sum(reshape(fitax      (:,:,:,1:N_w) .* dist_w, [], N_w), 1), ...
             sum(reshape(fitaxss    (:,:,1:N_r)   .* dist_r, [], N_r), 1)];

SSrev     = [sum(reshape(fsstax     (:,:,:,1:N_w) .* dist_w, [], N_w), 1), ...
             sum(reshape(fsstaxss   (:,:,1:N_r)   .* dist_r, [], N_r), 1)];

Fcaptax   = [sum(reshape(fcaptax    (:,:,:,1:N_w) .* dist_w, [], N_w), 1), ...
             sum(reshape(fcaptaxss  (:,:,1:N_r)   .* dist_r, [], N_r), 1)];

SSexp     = sum(reshape(repmat(ben', [nk,1,N_r]) .* dist_r, [], N_r), 1);

end