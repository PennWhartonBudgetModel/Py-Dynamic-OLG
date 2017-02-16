clear; clc;

T_life  = 80;
T_work  = 47;
T_model = 75;

s = load('params.mat');
bgrid = s.bgrid;

for startyear = (-T_life+1):(T_model-1)
    
    T_past   = max(-startyear, 0);
    T_shift  = max(+startyear, 0);
    T_active = min(startyear+T_life, T_model) - T_shift;
    
    for idem = 1:2
        
        file = fullfile('Freeze', 'Cohorts', sprintf('cohort=%+03d_idem=%u.mat', startyear, idem));
        
        opt = load(file);
        B = opt.B;
        
        t0 = min(max(T_work + startyear, 0), T_work) + 1;
        span = size(B(:,:,:,t0:end), 4);
        B(:,:,:,t0:end) = repmat(reshape(bgrid, [1,1,5,1]), [10,4,1,span]);
        
        save(file, '-append', 'B')
        
    end
    
end