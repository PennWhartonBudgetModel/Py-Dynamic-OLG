% Pete | 2016-09-28
% 
% Calculate static aggregates for all transition paths.
% 
% 


function [] = generate_static_aggregates_all()

parfor inddeep = 1:16
    
    deep_params = inddeep_to_params(inddeep);
    
    for plan = {'base', 'trump', 'clinton', 'ryan'}
        
        % Open economy
        generate_static_aggregates(deep_params, plan{1})
        
        for gcut = [+0.10, +0.05, +0.00, -0.05]
            
            % Closed economy
            generate_static_aggregates(deep_params, plan{1}, gcut)
            
        end
        
    end
    
end

end