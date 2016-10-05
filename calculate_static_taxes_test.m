% Jorge | 2016-09-29
% 
% Test script for static aggregate tax calculation function.  Used to expedite C code generation.
% 
% 


clear
clc

load('calculate_static_taxes_freeze.mat')

[fedincome, fedincomess, fitax, fitaxss, fsstax, fsstaxss, fcaptax, fcaptaxss] ...
  ...
  = calculate_static_taxes(...
      T, Tr, Tss, birthyear, ...
      kgrid, nb, nk, nz, idem, ...
      mpci, rpci, cap_tax_share, ss_tax_cred, ...
      tau_cap, tau_capgain, ss_tax, taxmax, z, ...
      wages, cap_shares, debt_shares, rate_caps, rate_govs, exp_subsidys, q_tobin, q_tobin0, ...
      avg_deduc, clinton, coefs, limit, X, ss_benefit, ...
      labopt);