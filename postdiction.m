function [vector_triggeredAt, vector_varianceTimeCourse, vector_triggeredAt2] = ...
    postdiction(  vector_responseTimes, vector_category, vector_accuracy, ...
                  categoryIDFrequent, windowSize, numTrials, numPadTrials  )
    % Data analysis designed to aid in informing how the realtime triggering setup should be designed.
    % I.e. gives useful information to help determine what the triggering constraints should be.
    
    % Used to determine trial triggering conditions. (+/- 7.5% VTC and +/- 1std RT)
    THRESHOLD_VAR = [15,85]; % percentage-based
    THRESHOLD_RTS = [-1,1]; % standard deviation-based
    
    
    % Shorten inputs.
    ws = windowSize;
    n = numTrials;
    np = numPadTrials;
    nm = (np+1):(n-np); % i.e. all the middle (hence 'nm') trials where triggering is possible
    rts = vector_responseTimes;
    ct  = vector_category;
    acc = vector_accuracy;
    FREQUENT = categoryIDFrequent;
    
    % Preallocate outputs.
    vtc = nan(length(nm),n); % this is (nm,n) because the entirety of vtc updates with each additional trial
    trg = nan(1,n);
    trg2 = nan(1,n);
    
    
    
    
    %%% processing raw rt vector %%%
    %
    rts( rts<0.1000 ) = NaN;
    
    INFREQUENT = setdiff([1,2],FREQUENT);
    rts( acc==0 | ct==INFREQUENT ) = NaN;
    
    trialsWithResponses = ~isnan(rts); % (The RT of a trial being NaN means that it wasn't responded to.)
    trialIDs = 1:n;
    linearFit = polyfit( trialIDs(trialsWithResponses) , rts(trialsWithResponses) , 1 );
    % Subtract predicted values from actual values to get detrended values.
    rts = rts - ( linearFit(1).*trialIDs + linearFit(2) );
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    % The logic here mirrors the real-time triggering design; see RealTimeSustAttnDisplay.m for more details.
    counter=1;
    for i = nm
        
        % Calculate VTC for each trial, using an expanding window.
        rt_mean = nanmean( rts(1:i) );
        rt_std  = nanstd(  rts(1:i) );
        vtc(counter,1:i) = abs(  (rts(1:i) - rt_mean) ./ rt_std  );
        
        % Determine whether a trigger (and what kind) could have potentially happened on this trial.
        % (I.e. if all trials in the preceding window are frequent, and if all of them were correctly responded to.)
        if all( ct((i-ws):(i-1))==FREQUENT ) && all( acc((i-ws):(i-1))==1 )
            
            % Calculate triggering threshold (based on THRESHOLD_VAR)
            % (Note that VTC is itself an expanding window, so the thresholds will change through the trials)
            threshold_var = prctile( rmmissing(vtc(counter,1:i)) , THRESHOLD_VAR );
                lower_vtc_threshold = threshold_var(1);
                upper_vtc_threshold = threshold_var(2);
            
            vtc_mean = mean( vtc(counter, (i-ws):(i-1) ) );
            
            % If the VTC sliding window average has broken either threshold, indicate it with -1 (lower) or 1 (upper),
            % otherwise, indicate a 0 (no trigger *due to not exceeding thresholds*, as opposed to failing the outer requirements).
            if vtc_mean < lower_vtc_threshold
                trg(i) = -1;
            elseif vtc_mean > upper_vtc_threshold
                trg(i) = 1;
            else
                trg(i) = 0;
            end
            
            
            % Now do RT-based triggering simulation.
            threshold_rts = THRESHOLD_RTS;
                lower_rts_threshold = threshold_rts(1);
                upper_rts_threshold = threshold_rts(2);
            
            rts_window_mean = nanmean( rts((i-ws):(i-1)) );
            
            standardised_rts_i = (rts_window_mean - rt_mean) / rt_std;
            
            if standardised_rts_i < lower_rts_threshold
                trg2(i) = -1;
            elseif standardised_rts_i > upper_rts_threshold
                trg2(i) = 1;
            else
                trg2(i) = 0;
            end
            
            
        end
        
        counter=counter+1;
        
    end
    
    
    vector_varianceTimeCourse = vtc(end,:); % change this output if we want to get all vtc rows
    % also note that this will output, by default, the vtc of the 500 trials; technically the first 50 should be nan (since the last 50 is also nan)
    vector_triggeredAt = trg;
    vector_triggeredAt2 = trg2;
    
end