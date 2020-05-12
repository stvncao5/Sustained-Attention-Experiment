function pValue = bootp(outputForAllBoots, outputBaseline, numBoots, direction)
    % 
    % Required Scripts: {computeCohen_d.m}
    % 
    % The general idea is to take some output of all bootstrap samples (e.g. the mean of all bootstrap samples)
    % and to find how many of those bootstrap samples' outputs are more/less than the output value which we are
    % comparing it to (this can be, for example, a value that specifies pure chance). In other words, we are
    % looking for the proportion of bootstrap samples whose output "deviates" (whether above or below depends on
    % the specified direction) from the baseline output that we would "expect".
    % 
    % Alternatively, if we want to figure out whether the distribution we see in the bootstrap is significantly
    % different from a distribution that represents the output of "pure chance", then Cohen's d will be implemented
    % between the two distributions (null vs. bootstrapped).
    % (Note that the null distribution will have to be provided in this case via the outputBaseline argument.)
    
    assert( ...
        strcmp(direction,'above') || ...
        strcmp(direction,'below') || ...
        strcmp(direction,'fromNull') ...
        , 'WARNING: No valid direction specified in bootp()' );
    
    if strcmp(direction,'above')
        
        pValue = sum( outputForAllBoots > outputBaseline ) / numBoots;
        if pValue == 0, pValue = (1/numBoots); end
        
    elseif strcmp(direction,'below')
        
        pValue = sum( outputForAllBoots < outputBaseline ) / numBoots;
        if pValue == 0, pValue = (1/numBoots); end
        
        
        
    elseif strcmp(direction,'fromNull')
        
        boot_dims = ndims(outputForAllBoots);
        base_dims = ndims(outputBaseline);
        vector_outputForAllBoots = reshape( outputForAllBoots, [1,numel(outputForAllBoots)] ); % do we want mean of each bootstrap sample instead?
        vector_outputBaseline    = reshape( outputBaseline   , [1,numel(outputBaseline   )] );
        
        d = computeCohen_d(vector_outputForAllBoots, vector_outputBaseline);
        n = 32; % terrible, but will do for now
        pseudo_zscore = d * sqrt(n);
        pValue = normcdf(pseudo_zscore);
        
        
    end
    
    
    
end

