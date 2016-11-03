%%
% Calculate static aggregates for all closed economies.
% 
%%


function [] = generate_static_aggregates_closed()

parfor inddeep = 1:16
    
    deep_params = inddeep_to_params(inddeep);
    
    for plan = {'base', 'trump', 'clinton', 'ryan'}
        for gcut = [+0.10, +0.05, +0.00, -0.05]
            
            % Skip non-zero gcut baseline runs
            if (strcmp(plan, 'base') && gcut ~= +0.00), continue, end
            
            % Generate static aggregates for closed economy run
            generate_static_aggregates(deep_params, plan{1}, gcut)
            
        end
    end
    
end

end