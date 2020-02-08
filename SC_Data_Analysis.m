%% Introduction
% This is the "master script" for the data analysis of a replication study.
% Authors: Megan deBettencourt, Monica Rosenberg, Steven Cao
% Required Scripts: {inpaint_nans.m, shadedErrorBar.m, dirP.m, A_Prime.m}
% 
% Additional Notes:
% Search '% #' to jump between various, section IDs (helps for quickly navigating through parts of the code)
% Search '%@' for all hardcoded variables that require manual specification.
% Search '% @' for all "todo" comments within in the code.
% Search '%***' for all parts of the code that are reserved for debugging (they are likely commented-out).
% 
% @ Big Todo List:
% Was recently given a heavy rewrite; need to validate all outputs
% Scaling of vtc_smooth is fairly off of what it should be
% Figure 2 replication (from Rosenberg et al, 2013) concerning comparison of errors between OTZ and ITZ is not working
    % (The reason likely being that the median split is derived from vtc_smooth which is currently behaving weirdly)
% Update A' function
% Update memory accuracy calculations to require high confidence (see above)
% Make a readme.md file mentioning how some things have been refactored, and what about this script needs to change to work with the original dataset
% Make section IDs + table of contents
% Finish filling out the big todo list
% 


%% Script Setup

% Flush session history
clear; clc; close all;

% Get main directory (it is assumed this script is running from the study's main directory)
projectDirectory = mfilename('fullpath');
% (Fragment into foldernames/filenames, get last item which is the filename of this script, subtract that, then put it all back together
% to get the script's home folder.)
projectDirectory = split(projectDirectory, {'/','\'}); projectDirectory = projectDirectory(1:end-1); projectDirectory = strjoin(projectDirectory, '/');
%@ Alternatively can specify absolute path manually for debugging purposes.
%projectDirectory = '/Users/rosenberglab/Desktop/discCPT';

% Also define the project's data directories accordingly, for convenience.
folderName_data = 'data';
folderName_rtdata = 'rtdata';
projectDirectory_data = fullfile(projectDirectory, folderName_data, '/');
projectDirectory_rtdata = fullfile(projectDirectory, folderName_rtdata, '/');


%% Preliminary Variable Setup

%%% Script Parameters
% (All subject data folders should be integer IDs)
shifts = -3:1:-1; %@ hardcoded
numShifts = numel(shifts);
stdWindow = 9; %@ hardcoded; keep in mind that the actual window size will always be +1, because we also include the RT of the target trial
% These parameters are determined by the settings of the experiment; here they are set to the default.
numAttTrialsTotal = 500;
proportionFreq = 0.90;      numTrialsFreq   = numAttTrialsTotal * proportionFreq;
proportionInfreq = 0.10;    numTrialsInfreq = numAttTrialsTotal * proportionInfreq;
numMemTrialsTotal = 200;

%%% Subject ID Selection
% (Generates a *precise* list of all available subject IDs.)
assert( ~isempty(isnumeric(str2num(char(dirP(projectDirectory_data)')))) , 'WARNING: Invalid subject data folder' ); %#ok<ST2NM>
subjectID_list = sort(str2num(char(dirP(projectDirectory_data)'))'); %#ok<ST2NM,TRSRT>

excludeID_list = [4,10,12,21,25]; %@ these are the IDs to be excluded
subjectID_list = setdiff(subjectID_list, excludeID_list);
numSubjects = length(subjectID_list);

%%% Bootstrap - Resample from Subject IDs
% (Note: this returns the ordered indices associated with each subject ID, and not the subject ID themselves,
% e.g. a vector of subjectIDs may be [1, 4, 7, 8, 9, 10], but its ordered indices are [1, 2, 3, 4, 5, 6])
numSamples_bootstrap = 100000;
subjectID_list_bootstrap = nan(numSamples_bootstrap, numSubjects); % preallocate
for i=1:size(subjectID_list_bootstrap,1)
    subjectID_list_bootstrap(i,:) = randsample(1:numSubjects, numSubjects, true)'; % randsample returns a column vector for some reason
end

%%% Preinitialise Data Collectors

%%%% (Attention Accuracy)
att_hit                 = nan(1,numSubjects); % correctly identifying frequent trials
att_falseAlarm          = nan(1,numSubjects); % mistaking infrequent trials (e.g. labeling them as frequent)
att_miss                = nan(1,numSubjects); % mistaking frequent trials (e.g. labeling them as infrequent)
att_correctReject       = nan(1,numSubjects); % correctly identifying infrequent trials (e.g. successfully overriding habitual keypress)

%%%% (Attention Accuracy - Replicating Figure 2a from 2013 AP&P Paper)
att_falseAlarm_ITZ      = nan(1,numSubjects); % out of all the infrequent trials, how many were mistaken AND "in the zone"?
att_falseAlarm_OTZ      = nan(1,numSubjects); % out of all the infrequent trials, how many were mistaken AND "out of the zone"?
att_miss_ITZ            = nan(1,numSubjects); % out of all the frequent trials, how many were mistaken AND "in the zone"?
att_miss_OTZ            = nan(1,numSubjects); % out of all the frequent trials, how many were mistaken AND "out of the zone"?

%%%% (Memory Accuracy)
memFreq_hit             = nan(1,numSubjects); % correctly remembering a "frequent-category" picture
memInfreq_hit           = nan(1,numSubjects); % correctly remembering an "infrequent-category" picture
memFreq_falseAlarm      = nan(1,numSubjects); % falsely remembering a "frequent-category" picture
memInfreq_falseAlarm    = nan(1,numSubjects); % falsely remembering an "infrequent-category" picture

%%%% (Attention Metrics - obtained from those sets of trials which precede an infrequent "target" trial)
% Note: When using nanmean() to compress a matrix into a scalar (2d into 0d), the result WILL depend on the order of
% "dimension compression"!
    % (Average ACROSS shifts first; this is the calculation used in code block #10 of the original python script)
    % ACROSS shifts == we want to know the average RT (or whatever value) of each iShift; each iShift
    % has only one value, which is an average across multiple trials
rts_before_acc          = nan(numSubjects,numShifts); % typical RTs before a correctly-responded trial
rts_before_inacc        = nan(numSubjects,numShifts); % typical RTs before an incorrectly-responded trial
var_before_acc          = nan(numSubjects,numShifts); % typical RT deviances before a correctly-responded trial
var_before_inacc        = nan(numSubjects,numShifts); % typical RT deviances before an incorrectly-responded trial
    % (Average WITHIN shifts first; this is the calculation used in code block #8 of the original python script)
    % WITHIN shifts == we want to know the average RT (or whatever value) for each target trial; each target trial
    % has only one value, which is an average of all of the target trial's shifts
comp_rts_before_acc     = cell(numSubjects,1); % average RT before each correctly-responded trial
comp_rts_before_inacc   = cell(numSubjects,1); % average RT before each incorrectly-responded trial
comp_var_before_acc     = cell(numSubjects,1); % average RT deviance before each correctly-responded trial
comp_var_before_inacc   = cell(numSubjects,1); % average RT deviance before each incorrectly-responded trial
    % (Shifts are inapplicable for standard deviation calculations.)
std_before_acc          = nan(numSubjects,1); % typical StD of RTs preceding a correctly-responded trial
std_before_inacc        = nan(numSubjects,1); % typical StD of RTs preceding an incorrectly-responded trial

%%%% (Analysis to relate Attention/RT to Memory/Accuracy - used for logistic regression)
rts_before_infreq           = nan(numSubjects,numShifts); % typical RT values for an infrequent trial
var_before_infreq           = nan(numSubjects,numShifts); % typical RT deviance for an infrequent trial
std_before_infreq           = nan(numSubjects,1        ); % typical RT StD values for an infrequent trial
comp_rts_before_infreq      = cell(numSubjects,1); % average RT value before each infrequent trial
comp_var_before_infreq      = cell(numSubjects,1); % average RT deviance before each infrequent trial
std_before_infreq_all       = cell(numSubjects,1); % average RT StD values before each infrequent trial
rts_before_infreq_estimated = cell(numSubjects,1); % logModel's RT value before each infrequent trial
var_before_infreq_estimated = cell(numSubjects,1); % logModel's RT deviation value before each infrequent trial
std_before_infreq_estimated = cell(numSubjects,1); % logModel's RT StD value before each infrequent trial
logistic_coefficients_rts   = nan(numSubjects,2); % first column == constant, second column == first coefficient
logistic_coefficients_var   = nan(numSubjects,2);
logistic_coefficients_std   = nan(numSubjects,2);

%%%% (Recent additions; experimental)
logistic_c1 = nan(numSubjects,3); % constant, first coefficient (RT), second coefficient (RT Var)
logistic_c2 = nan(numSubjects,4); % constant, first coefficient (RT), second coefficient (RT Var), third coefficient (interaction)
ttest2_pValue_std_AccVsInacc = nan(numSubjects,1);
rt_rtVar_corr_collector = cell(numSubjects,1);


%%% Set Script Flags
% plotRTs:           Plot results based on raw RTs (i.e. for comparability to deBettencourt et al, 2018).
% plotRTDevs:        Plot results based on RT Deviations (i.e. relating to RT variability).
% plotOther:         Plot results based on other metrics (e.g. standard deviation of RT windows preceding infrequent trials).
% detrend:           Remove potential time-based effects on attention and/or response times (e.g. vigilance decrement).
  % (Set to true for comparability to deBettencourt et al, 2018, otherwise leave as false.)
% includeAllTrials:  Include the RTs of infrequent and/or incorrect trials when collating the RTs of trials which precede the target trial of interest.
  % (Set to true for comparability to deBettencourt et al, 2018, otherwise leave as false.)
% interpolateRTs:    Replace trials that do not have a respone time with an interpolated value. Note that this also means VTC will be interpolated.
flag_plotRTs            = true;
flag_plotRTDevs         = true;
flag_plotOther          = true;
flag_detrend            = false;
flag_includeAllTrials   = false;
flag_interpolateRTs     = false;


%% Subject Analysis
for iSubject = 1:length(subjectID_list)
    % Notes on conventions:
    % iSubject is an index for filling subject-related data variables. This is always in order, e.g. [1, 2, 3, ...]
    % subjectID is the actual subject's ID. This is not always in order, e.g. [51, 52, 53, 59, ...]
    % Accuracy is always encoded as either 0 or 1, and not some other value (e.g. NaN)
    
    subjectID = subjectID_list(iSubject);
    subjectID_dir = fullfile(projectDirectory_data, num2str(subjectID), '/');
    
    % Load subject's attention data.
        file = dir([subjectID_dir, 'attndata_*']);
        assert(numel(file)<=1,['ERROR: More than one attention datafile for subject ID: ' num2str(subjectID)])
        load([subjectID_dir, file.name]); % this will load in the variable called 'attnData'
    % Load subject's memory data.
        file = dir([subjectID_dir, 'memdata_*']);
        assert(numel(file)<=1,['ERROR: More than one memory datafile for subject ID: ' num2str(subjectID)])
        load([subjectID_dir, file.name]); % this will load in the variable called 'memData'
    
    %%% Attention Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Pull relevant data to be analysed.
    att_accs = attnData.accs;
    att_categs = attnData.categs;
    att_RTs = attnData.rts;
    att_trialIDs = attnData.trial;
    
    % If we want detrended RTs instead (to account for natural attention decrements over time/trials), then
    % convert RTs into RT-residuals (THIS MEANS 'att_RTs' WILL NOW REFER TO THE RESIDUALS AND NOT TO THE RAW RT VALUES).
    if flag_detrend==true
        % detrend() not recommended because it doesn't work with nan values.
        
        % Generated predicted values (based on a linear model).
        % (Think of the trialIDs as time points.)
        trialsWithResponses = ~isnan(att_RTs); % (The RT of a trial being NaN means that it wasn't responded to.)
        linearFit = polyfit( att_trialIDs(trialsWithResponses) , att_RTs(trialsWithResponses) , 1 );
        % Subtract predicted values from actual values to get detrended values.
        att_RTs = att_RTs - ( linearFit(1).*att_trialIDs + linearFit(2) );
    end

    % Identify frequent and infrequent category (category ID is either 1 or 2)
    FREQUENT = mode(attnData.categs);
    INFREQUENT = setdiff([1,2],FREQUENT);
    
    % Get subject's accuracy-related stats from the attention task.
    att_hit(iSubject)            = mean( att_accs(att_categs==FREQUENT  )==1 );
    att_falseAlarm(iSubject)     = mean( att_accs(att_categs==INFREQUENT)==0 );
    att_miss(iSubject)           = mean( att_accs(att_categs==FREQUENT  )==0 );
    att_correctReject(iSubject)  = mean( att_accs(att_categs==INFREQUENT)==1 );
    
    % Find the trial IDs where subject correctly identified and where subject incorrectly identified infrequent trials.
    trial_infreq_acc   = att_trialIDs( att_accs==1 & att_categs==INFREQUENT );
    trial_infreq_inacc = att_trialIDs( att_accs~=1 & att_categs==INFREQUENT );
    
    % When looking at any infrequent trial, we will look at the response times of the previous (3) trials.
    % If we don't want to include all trials, then we only consider the ones that were both correctly responded to and also in the frequent category.
    % The others will be "ignored" by replacing them with NaN (which gets ignored in e.g., 'nanmean()').
    if flag_includeAllTrials==false
        
        validShift_RTs = att_RTs;
        validShift_RTs( att_accs==0 | att_categs==INFREQUENT ) = NaN;
        
    % If we do want to include all trials, then all RTs are considered "valid": none of them will be replaced with NaN.
    elseif flag_includeAllTrials==true
        
        validShift_RTs = att_RTs;
        
    end
    
    % Calculate the variance time course (VTC) for an alternative measure (the main measure being the raw RTs).
    % (Method of calculating VTC is derived from Rosenberg et al, 2013.)
    rt_mean   = nanmean(validShift_RTs);
    rt_std    = nanstd(validShift_RTs);
    
    % (Old way of calculating VTC; deprecated)
    %rt_median = nanmedian(validShift_RTs);
    %vtc = abs(validShift_RTs - rt_median);
    
    % Interpolate RT (and consequently, VTC) values.
    % Note: this affects all data analysis. Leave this on 'false' to use interpolation only for calculating smoothed VTC.
    if flag_interpolateRTs==true
        
        validShift_RTs = inpaint_nans(validShift_RTs,4);
        
    end
    
    % Calculate VTC.
    % (Can't use zscore() because it is affected by NaNs.)
    vtc = abs(  (validShift_RTs - rt_mean) ./ rt_std  );
    
    % Calculate smoothed VTC (mainly for Figure 2 of Rosenberg et al, 2013).
    interpolated_RTs = inpaint_nans(att_RTs,4); % note that this will include RTs of infrequent and/or error trials, i.e. ignores validShift_RTs
    interpolated_vtc = abs( (interpolated_RTs - rt_mean) ./ rt_std ); % note that this is not the same as abs(zscore(interpolated_RTs))
    kernelSize  = 10;
    gaussKernel = gausswin(kernelSize);
    vtc_smooth  = zscore(filtfilt(gaussKernel,1,interpolated_vtc)); % @todo: figure out how to fix the scaling so we can remove the need to zscore again
    vtc_zoneSplit = median(vtc); % for getting the median split of the VTC in order to determine which values correspond to ITZ and OTZ
    
    
    % (Data array preallocation)
    numTrial_infreq_acc   = length(trial_infreq_acc);
    numTrial_infreq_inacc = length(trial_infreq_inacc);
    
    all_rts_before_acc   = nan(numTrial_infreq_acc,   numShifts);
    all_rts_before_inacc = nan(numTrial_infreq_inacc, numShifts);
    all_var_before_acc   = nan(numTrial_infreq_acc,   numShifts);
    all_var_before_inacc = nan(numTrial_infreq_inacc, numShifts);
    all_std_before_acc   = nan(numTrial_infreq_acc,   1        );
    all_std_before_inacc = nan(numTrial_infreq_inacc, 1        );
    
    % Now retrieve our metrics of interest.
    counter=1;
    for iShift = shifts
        
        % Grab all the "shifted trials" (i.e. the trials preceding their respective target trial by "i" shifts).
        shifted_trials_acc   = trial_infreq_acc    + iShift;
        shifted_trials_inacc = trial_infreq_inacc  + iShift;
        % Find out which of the shift-trial IDs are valid indices (i.e. won't make us index out of bounds by having a value of 0 or less).
        whereValid_shiftedTrials_acc       = shifted_trials_acc>=1;
        whereValid_shiftedTrials_inacc     = shifted_trials_inacc>=1;
        % For all of the shift-trial IDs that are valid, get its corresponding ordered index (i.e. for data entry into a new array).
        validIndicesOf_shiftedTrials_acc   = find(whereValid_shiftedTrials_acc);
        validIndicesOf_shiftedTrials_inacc = find(whereValid_shiftedTrials_inacc);
        % Truncate all invalid shift-trial IDs from the list of shift-trial IDs.
        valid_shifted_trials_acc            = shifted_trials_acc  (whereValid_shiftedTrials_acc);
        valid_shifted_trials_inacc          = shifted_trials_inacc(whereValid_shiftedTrials_inacc);
        
        % Fill in the data array's ordered index with its respective trial's corresponding RT value.
        % Anything which is ignored (those indices which correspond to shifted indices that are negative) is left as a NaN.
            % (Note: these operations involve a subtle coercion from a row vector to a column vector)
        all_rts_before_acc  (validIndicesOf_shiftedTrials_acc  ,counter) = validShift_RTs(valid_shifted_trials_acc);
        all_rts_before_inacc(validIndicesOf_shiftedTrials_inacc,counter) = validShift_RTs(valid_shifted_trials_inacc);
        % Do the same but for RT variability.
        all_var_before_acc  (validIndicesOf_shiftedTrials_acc  ,counter) = vtc(valid_shifted_trials_acc);
        all_var_before_inacc(validIndicesOf_shiftedTrials_inacc,counter) = vtc(valid_shifted_trials_inacc);
        
        counter=counter+1; % counter corresponds to iShift, except they're positive integers (for indexing)
    end
    
    % Now compute the standard deviations of the preceding trials (default window size is 9+1) for each infrequent trial.
    for i = 1:numTrial_infreq_acc
        % If we can index backwards by the full stdWindow length, then get the StD of the window; otherwise leave as NaN.
        if (trial_infreq_acc(i) - stdWindow) >= 1
            all_std_before_acc(i) = nanstd(   validShift_RTs( trial_infreq_acc(i)-stdWindow : trial_infreq_acc(i) )   );
        end
    end
    for i = 1:numTrial_infreq_inacc
        % If we can index backwards by the full stdWindow length, then get the StD of the window; otherwise leave as NaN.
        if (trial_infreq_inacc(i) - stdWindow) >= 1
            all_std_before_inacc(i) = nanstd(   validShift_RTs( trial_infreq_inacc(i)-stdWindow : trial_infreq_inacc(i) )   );
        end
    end
    
    
    % Calculate some summary statistics for RTs on the current subject ID.
    rts_before_acc(iSubject,:)      = nanmean(all_rts_before_acc      );
    rts_before_inacc(iSubject,:)    = nanmean(all_rts_before_inacc    );
    comp_rts_before_acc{iSubject}   = nanmean(all_rts_before_acc   , 2);
    comp_rts_before_inacc{iSubject} = nanmean(all_rts_before_inacc , 2);
   
    var_before_acc(iSubject,:)      = nanmean(all_var_before_acc      );
    var_before_inacc(iSubject,:)    = nanmean(all_var_before_inacc    );
    comp_var_before_acc{iSubject}   = nanmean(all_var_before_acc   , 2);
    comp_var_before_inacc{iSubject} = nanmean(all_var_before_inacc , 2);
    
    std_before_acc(iSubject)        = nanmean(all_std_before_acc      );
    std_before_inacc(iSubject)      = nanmean(all_std_before_inacc    );
    
    % @todo: write this better
    % append to an array a p-value of a ttest2 for this subject
    % the ttest2 asks, "is there a significant difference in RT std values between the accurate and the inaccurate trials?"
    % aka, "are the RT standard deviations of accurate trials significantly different from those of inaccurate trials?"
    % there will be some nan values that are removed first - these are the trials which we couldn't get an stdWindow for
    [~,ttest2_pValue_std_AccVsInacc(iSubject),~,~]=ttest2(rmmissing(all_std_before_acc),rmmissing(all_std_before_inacc));
    
    
    %%% Memory Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mem_accs = memData.accs;
    mem_categs = memData.memCategs; % 1 == frequent+old, 2 == infrequent+old, 3 == frequent+new, 4 == infrequent+new
    mem_ratings = memData.rating; % 1 == definitely new, 2 == maybe new, 3 == maybe old, 4 == definitely old
    mem_trialIDs = memData.trial;
    
    memFreq_hit         (iSubject) = mean(mem_accs(mem_categs==1)==1);
    memInfreq_hit       (iSubject) = mean(mem_accs(mem_categs==2)==1);
    memFreq_falseAlarm  (iSubject) = mean(mem_accs(mem_categs==3)==0);
    memInfreq_falseAlarm(iSubject) = mean(mem_accs(mem_categs==4)==0);
    % Trials are only considered "remembered" if the responses were also confident:
    memFreq_remembered  (iSubject) = mean(mem_accs(mem_categs==1 & mem_ratings==4)==1);
    memInfreq_remembered(iSubject) = mean(mem_accs(mem_categs==2 & mem_ratings==4)==1);
        
    % Get only the infrequent+old trial IDs from the attention task.
    % (Named 'both' because it's both acc and inacc.)
    trial_infreq_both = memData.attOrder(mem_categs==2);
    
    % (Data array preallocation)
    numTrial_infreq_both = size(trial_infreq_both,2);
    
    all_rts_before_infreq = nan(numTrial_infreq_both , numShifts);
    all_var_before_infreq = nan(numTrial_infreq_both , numShifts);
    all_std_before_infreq = nan(numTrial_infreq_both , 1        );
    
    % Now retrieve our metrics of interest.
    counter=1;
    for iShift = shifts
        
        % Grab all the "shifted trials" (i.e. the trials preceding their respective target trial by "i" shifts).
        shifted_trials_both = trial_infreq_both + iShift;
        % Find out which of the shift-trial IDs are valid indices (i.e. won't make us index out of bounds by having a value of 0 or less).
        whereValid_shiftedTrials_both = shifted_trials_both>=1;
        % For all of the shift-trial IDs that are valid, get its corresponding ordered index (i.e. for data entry into a new array).
        validIndicesOf_shiftedTrials_both = find(whereValid_shiftedTrials_both);
        % Truncate all invalid shift-trial IDs from the list of shift-trial IDs.
        valid_shifted_trials_both = shifted_trials_both(whereValid_shiftedTrials_both);
        
        % Fill in the data array's ordered index with its respective trial's corresponding RT value.
        % Anything which is ignored (those indices which correspond to shifted indices that are negative) is left as a NaN.
            % (Note: these operations involve a subtle coercion from a row vector to a column vector)
        all_rts_before_infreq(validIndicesOf_shiftedTrials_both , counter) = validShift_RTs(valid_shifted_trials_both);
        % Do the same but for RT variability.
        all_var_before_infreq(validIndicesOf_shiftedTrials_both , counter) = vtc(valid_shifted_trials_both);
        
        counter=counter+1; % counter corresponds to ishift, except with positive integers
    end
    
    % Now compute the standard deviations of the preceding trials (default window size is 9+1) for each infrequent trial.    
    for i = 1:numTrial_infreq_both
        % If we can index backwards by the full stdWindow length, then get the StD of the window; otherwise leave as NaN.
        if (trial_infreq_both(i) - stdWindow) >= 1
            all_std_before_infreq(i) = nanstd(   validShift_RTs( trial_infreq_both(i)-stdWindow : trial_infreq_both(i) )   );
        end
    end
    
    
    % Calculate some summary statistics for RTs on the current subject ID.
    rts_before_infreq(iSubject,:) = nanmean(all_rts_before_infreq);   
    comp_rts_before_infreq{iSubject} = nanmean(all_rts_before_infreq,2);
    
    var_before_infreq(iSubject,:) = nanmean(all_var_before_infreq);
    comp_var_before_infreq{iSubject} = nanmean(all_var_before_infreq,2);
    
    std_before_infreq(iSubject) = nanmean(all_std_before_infreq);
    std_before_infreq_all{iSubject} = all_std_before_infreq; % will need this for comparing predicted values against actual values wit logistic regression
    
    
    % Arrange vectors in the correct order for use in logistic regression modeling.
    % We need an inputVariable vector (e.g. RTs), an outputVariable vector (i.e. accuracy), and a vector of ones (because MATLAB is picky).
    % Additionally, all three vectors should have the same length (numInfreqTrials==50) and be column vectors.
    
    % Output: Accuracy Vector (boolean vector)
    accuracyVector = mem_ratings(mem_categs==2); % get only the trials of interest (old+infrequent)
    accuracyVector(accuracyVector~=4)=0; accuracyVector(accuracyVector==4)=1; % only confident responses (4) are considered correct
    accuracyVector=accuracyVector'; % transpose from row to column
    
    % Input: RT Vector (numeric vector)
    inputVector_rts = comp_rts_before_infreq{iSubject};
    
    % Input: VTC Vector (numeric vector)
    inputVector_var = comp_var_before_infreq{iSubject};
    
    % Input: RT StD Vector (numeric vector)
    inputVector_std = std_before_infreq_all{iSubject};
    
    % Note: All the vectors above already correspond to the correct (i.e. the same) trial IDs, so we can pair the vectors off just fine.
    
    % Auxiliary: "Number of times each value of the input variable is repeated" Vector (numeric vector)
    testedVector = ones(length(accuracyVector),1);
    
    
    % Logistic Regression, RTs
    % Resort in ascending order for the inputVector, so we can get the model's estimated values in order (for clean plotting).
    resorted = sortrows([inputVector_rts, accuracyVector, testedVector]);
    % Delete invalid trials (i.e. rows which have NaNs, i.e. rows which have no input value)
    resorted = rmmissing(resorted);
    % Generate arguments for glmfit() and glmval()
    logModel_rts_input  = resorted(:,1);
    logModel_rts_output = resorted(:,2);
    logModel_rts_aux    = resorted(:,3);
    % Generate the logistic regression model
    logModel_rts_coefficients = glmfit( logModel_rts_input, [logModel_rts_output, logModel_rts_aux], 'binomial', 'link', 'logit');
    % Generate the model's predicted values
    logModel_rts_predicted    = glmval( logModel_rts_coefficients, logModel_rts_input, 'logit');
    
    % Logistic Regression, VTC
    % Resort in ascending order for the inputVector, so we can get the model's estimated values in order (for clean plotting).
    resorted = sortrows([inputVector_var, accuracyVector, testedVector]);
    % Delete invalid trials (i.e. rows which have NaNs, i.e. rows which have no input value)
    resorted = rmmissing(resorted);
    % Generate arguments for glmfit() and glmval()
    logModel_var_input  = resorted(:,1);
    logModel_var_output = resorted(:,2);
    logModel_var_aux    = resorted(:,3);
    % Generate the logistic regression model
    logModel_var_coefficients = glmfit( logModel_var_input, [logModel_var_output, logModel_var_aux], 'binomial', 'link', 'logit');
    % Generate the model's predicted values
    logModel_var_predicted    = glmval( logModel_var_coefficients, logModel_var_input, 'logit');
    
    % Logistic Regression, StD
    % Resort in ascending order for the inputVector, so we can get the model's estimated values in order (for clean plotting).
    resorted = sortrows([inputVector_std, accuracyVector, testedVector]);
    % Delete invalid trials (i.e. rows which have NaNs, i.e. rows which have no input value)
    resorted = rmmissing(resorted);
    % Generate arguments for glmfit() and glmval()
    logModel_std_input  = resorted(:,1);
    logModel_std_output = resorted(:,2);
    logModel_std_aux    = resorted(:,3);
    % Generate the logistic regression model
    logModel_std_coefficients = glmfit( logModel_std_input, [logModel_std_output, logModel_std_aux], 'binomial', 'link', 'logit');
    % Generate the model's predicted values
    logModel_std_predicted    = glmval( logModel_std_coefficients, logModel_std_input, 'logit');
    
    % Append model results.
    logistic_coefficients_rts(iSubject,:) = logModel_rts_coefficients;
    rts_before_infreq_estimated{iSubject} = logModel_rts_predicted;
    
    logistic_coefficients_var(iSubject,:) = logModel_var_coefficients;
    var_before_infreq_estimated{iSubject} = logModel_var_predicted;
    
    logistic_coefficients_std(iSubject,:) = logModel_std_coefficients;
    std_before_infreq_estimated{iSubject} = logModel_std_predicted;
    
    
    % @todo: insert more plotting-related flags here
    
    % Plot logistic regression on raw RTs for each individual.
    if false
        figure('Position', [10 10 900 600]);
        hold on;
        plot(logModel_rts_input, logModel_rts_output, 'bs', logModel_rts_input, logModel_rts_predicted, 'r-', 'LineWidth', 2)
        set(gca,'fontsize',24);
        title("RT Attention & Memory Performance")
        xlabel("RT"); ylabel("Memory Performance");
        set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);
    end
    
    % Plot logistic regression on var for each individual.
    if false
        figure('Position', [10 10 900 600]);
        hold on;
        plot(logModel_var_input, logModel_var_output, 'bs', logModel_var_input, logModel_var_predicted, 'r-', 'LineWidth', 2)
        set(gca,'fontsize',24);
        title("RT Var Attention & Memory Performance")
        xlabel("RT Var"); ylabel("Memory Performance");
        set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);
    end
    
    % Plot logistic regression on std for each individual.
    if false
        figure('Position', [10 10 900 600]);
        hold on;
        plot(logModel_std_input, logModel_std_output, 'bs', logModel_std_input, logModel_std_predicted, 'r-', 'LineWidth', 2)
        set(gca,'fontsize',24);
        title("RT Std Attention & Memory Performance")
        xlabel("RT Std"); ylabel("Memory Performance");
        set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);
    end
    
    
    %% Memory Analysis (Extra)
    
    % Logistic regression, RTs+Var & RTs+Var+Interaction
    % Resort in ascending order for the inputVector (RTs), so we can get the model's estimated values in order (in case we want to plot for some reason).
    resorted = sortrows([inputVector_rts, inputVector_var, accuracyVector, testedVector]);
    % Delete invalid trials (i.e. rows which have NaNs, i.e. rows which have no input value)
    resorted = rmmissing(resorted);
    % Generate arguments for glmfit() and glmval()
    logModel_combined_input_rts = resorted(:,1);
    logModel_combined_input_var = resorted(:,2);
    logModel_combined_output    = resorted(:,3);
    logModel_combined_aux       = resorted(:,4);
    % Generate both logistic regression models
    logModel_c1_coefficients = glmfit( [logModel_combined_input_rts, logModel_combined_input_var], ...
        [logModel_combined_output, logModel_combined_aux], 'binomial', 'link', 'logit');
    logModel_c2_coefficients = glmfit( [logModel_combined_input_rts, logModel_combined_input_var, logModel_combined_input_rts .* logModel_combined_input_var], ...
        [logModel_combined_output, logModel_combined_aux], 'binomial', 'link', 'logit');
    % @todo: generate predicted values, in case we want to plot or something
    
    % Append model results.
    logistic_c1(iSubject,:) = logModel_c1_coefficients;
    logistic_c2(iSubject,:) = logModel_c2_coefficients;
    
    
    %% Experimental/Other Analyses
    
    % Figure 2a from Rosenberg et al 2013 APP paper
    % @todo: set plotting flags here
    if true
        figure('Position', [10 10 900 600]);
        hold on;
        % plot the VTC using a light-gray line
        plot(att_trialIDs, vtc                                    , '-', 'Color', [0.75,0.75,0.75]);
        % plot a smoothed VTC using a red, thick line
        plot(att_trialIDs, vtc_smooth                             , 'r-', 'LineWidth', 3);
        % draw a black horizontal line through the time series, showing the median split of the vtc
        plot(att_trialIDs, vtc_zoneSplit*ones(1,numAttTrialsTotal), '-', 'Color', [0,0,0], 'LineWidth', 3);
        % other plotting stuff
        set(gca,'fontsize',24,'ylim',[-1,6]);
        title("Variance Time Course");
        xlabel("Trial"); ylabel("VTC");
        legend('VTC','Smoothed VTC');
        legend boxoff
    end
    
    % Get the commission/omission (i.e. falseAlarm/miss) error rates of a given category (i.e. frequent or infrequent).
    att_falseAlarm_ITZ(iSubject) = sum( att_accs(att_categs==INFREQUENT & vtc_smooth<vtc_zoneSplit)==0 ) / numTrialsInfreq;
    att_falseAlarm_OTZ(iSubject) = sum( att_accs(att_categs==INFREQUENT & vtc_smooth>vtc_zoneSplit)==0 ) / numTrialsInfreq;
    att_miss_ITZ(iSubject)       = sum( att_accs(att_categs==FREQUENT   & vtc_smooth<vtc_zoneSplit)==0 ) / numTrialsFreq;
    att_miss_OTZ(iSubject)       = sum( att_accs(att_categs==FREQUENT   & vtc_smooth>vtc_zoneSplit)==0 ) / numTrialsFreq;
    
    % @todo: write this better
    % basically gets the correlation coefficient between RT and RT var
    f = fit( rmmissing(validShift_RTs)', rmmissing(vtc)', 'poly1' );
    rt_rtVar_corr_collector{iSubject} = f;
    % plot RT against RT var, along with the correlation coefficient line
    if false
        figure('Position', [10 10 900 600]);
        hold on;
        plot(validShift_RTs, vtc, 'o')
        xSpace = xlim;
        line_xSpace = linspace(xSpace(1), xSpace(2));
        line_ySpace = linspace(xSpace(1)*f.p1, xSpace(2)*f.p1);
        plot(line_xSpace, line_ySpace, 'r-', 'LineWidth', 3)
        set(gca,'fontsize',24);
        title("RTs & RT Variances")
        xlabel("RTs"); ylabel("RT Variances");
    end
    
    % #QBA
    % quartile based analyses
    numAttTrialsQuartile = numAttTrialsTotal/4; % will probably do something weird if this doesn't return an integer value
    for iQuartile = 1:4
        t = numAttTrialsQuartile*(iQuartile-1)+(1:numAttTrialsQuartile);
        temp_accs = attnData.accs(t);
        temp_rts  = attnData.rts(t);
        temp_rare = find(attnData.categs(t)==INFREQUENT);
        temp_freq = find(attnData.categs(t)==FREQUENT);
        
        quartile_commission_errors(iSubject,iQuartile) = mean(temp_accs(temp_rare)==0)*100;
        quartile_ommission_errors(iSubject,iQuartile)  = mean(temp_accs(temp_freq)==0)*100;
        quartile_rt(iSubject,iQuartile)                = nanmean(temp_rts)*1000;
        quartile_rtvar(iSubject,iQuartile)             = nanstd(temp_rts)*1000;
    end

    
end


%% Plotting

% Figures 2b and 2c from Rosenberg et al, 2013 (MATLAB does not support error bars for histograms)
figure('Position', [10 10 900 600]);

subplot(1,2,1);
set(gca,'fontsize',24);
hold on;
errorbar(1,nanmean(att_falseAlarm_ITZ),nanstd(att_falseAlarm_ITZ)/sqrt(numSubjects),'go','linewidth',4)
errorbar(2,nanmean(att_falseAlarm_OTZ),nanstd(att_falseAlarm_OTZ)/sqrt(numSubjects),'ro','linewidth',4)
title("False Alarms: ITZ vs OTZ")
set(gca,'xtick',[],'xlim',[0.0,3.0],'ylim',[0.0,0.5],'ytick',[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
ylabel('Error Rate (%)');
xlabel('');
legend('In the zone','Out of the zone');
legend boxoff

subplot(1,2,2);
set(gca,'fontsize',24);
hold on;
errorbar(1,nanmean(att_miss_ITZ),nanstd(att_miss_ITZ)/sqrt(numSubjects),'go','linewidth',4)
errorbar(2,nanmean(att_miss_OTZ),nanstd(att_miss_OTZ)/sqrt(numSubjects),'ro','linewidth',4)
set(gca,'xtick',[],'xlim',[0.0,3.0],'ylim',[0.0,0.05],'ytick',[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
title("Misses: ITZ vs OTZ")
ylabel('Error Rate (%)');
xlabel('');

hold off;
% @todo: use a paired ttest to compare ITZ against OTZ?


% Figure 3 from Rosenberg et al, 2013 (code was largely taken from MdB)
% (Note that the detrend flag MUST be off for this to make sense)
% #QBA
figure;

%figure 3a from rosenberg ap&p
subplot(2,2,1);
hold on;
shadedErrorBar(1:4,mean(quartile_commission_errors),std(quartile_commission_errors)/sqrt(numSubjects),'lineprops','k');
plot(1:4,mean(quartile_commission_errors),'k');
set(gca,'ylim',[25,40],'ytick',[25 30 35 40]);
ylabel('Commission errors')
xlabel('Quartile')

%figure 3b from rosenberg ap&p
subplot(2,2,2);
hold on;
shadedErrorBar(1:4,mean(quartile_rtvar),std(quartile_rtvar)/sqrt(numSubjects),'lineprops','k');
plot(1:4,mean(quartile_rtvar),'k');
ylabel('RT variability')
xlabel('Quartile')

%figure 3c from rosenberg ap&p
subplot(2,2,3);
hold on;
shadedErrorBar(1:4,mean(quartile_ommission_errors),std(quartile_ommission_errors)/sqrt(numSubjects),'lineprops','k');
plot(1:4,mean(quartile_ommission_errors),'k');
set(gca,'ylim',[0,5],'ytick',[1 2 3 4 5]);
ylabel('Omission errors')
xlabel('Quartile')

%figure 3d from rosenberg ap&p
subplot(2,2,4);
hold on;
shadedErrorBar(1:4,mean(quartile_rt),std(quartile_rt)/sqrt(numSubjects),'lineprops','k');
plot(1:4,mean(quartile_rt),'k');
set(gca,'xlim',[.5,4.5],'xtick',[1 2 3 4]);
ylabel('RT')
xlabel('Quartile')


% Plot logistic regression for all subjects (RTs)
% Note: rmmissing() needed for the actual values (the inputVector) in order to match the size of the estimated values vector
figure('Position', [10 10 900 600]);
hold on;
for iSubject = 1:numSubjects
    sorted_actualValues = rmmissing(sort(comp_rts_before_infreq{iSubject}));
    plot(sorted_actualValues, rts_before_infreq_estimated{iSubject}, '-','LineWidth',2);
end
set(gca,'fontsize',24);
title("Logistic Regression, all subject RTs"); xlabel("RTs (s)"); ylabel("Memory Performance"); set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);

% Plot logistic regression for all subjects (Var)
figure('Position', [10 10 900 600]);
hold on;
for iSubject = 1:numSubjects
    sorted_actualValues = rmmissing(sort(comp_var_before_infreq{iSubject}));
    plot(sorted_actualValues, var_before_infreq_estimated{iSubject}, '-','LineWidth',2);
end
set(gca,'fontsize',24);
title("Logistic Regression, all subject RT Vars"); xlabel("RT Vars (s)"); ylabel("Memory Performance"); set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);

% Plot logistic regression for all subjects (StD)
figure('Position', [10 10 900 600]);
hold on;
for iSubject = 1:numSubjects
    sorted_actualValues = rmmissing(sort(std_before_infreq_all{iSubject}));
    plot(sorted_actualValues, std_before_infreq_estimated{iSubject}, '-','LineWidth',2);
end
set(gca,'fontsize',24);
title("Logistic Regression, all subject RT StDs"); xlabel("RT StDs (s)"); ylabel("Memory Performance"); set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);


%plot across all subjects (note: the difference between detrended and non-detrended plots is just the y-axis scaling)
titleAverage = ["Average"];
titleSpaghetti = ["All Subjects"];

% @todo: probably should rework this code block
if flag_detrend == true

    if flag_plotRTDevs == true
        
        figure('Position', [10 10 900 600]);
        
        subplot(1,2,1);
        set(gca,'fontsize',24);
        hold on;
        errorbar(shifts,nanmean(var_before_acc),nanstd(var_before_acc)/sqrt(numSubjects),'b','linewidth',4)
        errorbar(shifts,nanmean(var_before_inacc),nanstd(var_before_inacc)/sqrt(numSubjects),'m','linewidth',4)
        title(titleAverage)
        set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',[0.000,1.200],'ytick',[0.000,0.400,0.800,1.200])
        ylabel('Detrended RT Devs, avg across subjects [s]');
        xlabel('Trials before infrequent trial');
        legend('Accurate','Inaccurate');
        legend boxoff
        
        subplot(1,2,2);
        set(gca,'fontsize',24);
        hold on;
        plot(shifts,var_before_acc,'b','linewidth',1)
        plot(shifts,var_before_inacc,'m','linewidth',1)
        plot(shifts,nanmean(var_before_acc),'b','linewidth',4)
        plot(shifts,nanmean(var_before_inacc),'m','linewidth',4)
        set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',[0.000,1.200],'ytick',[0.000,0.400,0.800,1.200])
        title(titleSpaghetti)
        xlabel('Trials before infrequent trial');
        
        hold off;
        
    end
    
    if flag_plotRTs == true
        
        figure('Position', [10 10 900 600]);
        
        subplot(1,2,1);
        set(gca,'fontsize',24);
        hold on;
        errorbar(shifts,nanmean(rts_before_acc),nanstd(rts_before_acc)/sqrt(numSubjects),'b','linewidth',4)
        errorbar(shifts,nanmean(rts_before_inacc),nanstd(rts_before_inacc)/sqrt(numSubjects),'m','linewidth',4)
        title(titleAverage)
        set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',[-0.200,0.200],'ytick',[-0.200,-0.100,0.000,0.100,0.200])
        ylabel('Detrended RTs, avg across subjects [s]');
        xlabel('Trials before infrequent trial');
        legend('Accurate','Inaccurate');
        legend boxoff
        
        subplot(1,2,2);
        set(gca,'fontsize',24);
        hold on;
        plot(shifts,rts_before_acc,'b','linewidth',1)
        plot(shifts,rts_before_inacc,'m','linewidth',1)
        plot(shifts,nanmean(rts_before_acc),'b','linewidth',4)
        plot(shifts,nanmean(rts_before_inacc),'m','linewidth',4)
        set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',[-0.200,0.200],'ytick',[-0.200,-0.100,0.000,0.100,0.200])
        title(titleSpaghetti)
        xlabel('Trials before infrequent trial');
        
        hold off;
        
    end

    
    
elseif flag_detrend == false
    
    if flag_plotRTDevs == true
        
        figure('Position', [10 10 900 600]);
        
        subplot(1,2,1);
        set(gca,'fontsize',24);
        hold on;
        errorbar(shifts,nanmean(var_before_acc),nanstd(var_before_acc)/sqrt(numSubjects),'b','linewidth',4)
        errorbar(shifts,nanmean(var_before_inacc),nanstd(var_before_inacc)/sqrt(numSubjects),'m','linewidth',4)
        title(titleAverage)
        set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',[0.000,1.200],'ytick',[0.000,0.400,0.800,1.200])
        ylabel('RT Deviations, avg across subjects [s]');
        xlabel('Trials before infrequent trial');
        legend('Accurate','Inaccurate');
        legend boxoff
        
        subplot(1,2,2);
        set(gca,'fontsize',24);
        hold on;
        plot(shifts,var_before_acc,'b','linewidth',1)
        plot(shifts,var_before_inacc,'m','linewidth',1)
        plot(shifts,nanmean(var_before_acc),'b','linewidth',4)
        plot(shifts,nanmean(var_before_inacc),'m','linewidth',4)
        set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',[0.000,1.200],'ytick',[0.000,0.400,0.800,1.200])
        title(titleSpaghetti)
        xlabel('Trials before infrequent trial');
        
        hold off;
        
    end
    
    if flag_plotRTs == true
        
        figure('Position', [10 10 900 600]);
        
        subplot(1,2,1);
        set(gca,'fontsize',24);
        hold on;
        errorbar(shifts,nanmean(rts_before_acc),nanstd(rts_before_acc)/sqrt(numSubjects),'b','linewidth',4)
        errorbar(shifts,nanmean(rts_before_inacc),nanstd(rts_before_inacc)/sqrt(numSubjects),'m','linewidth',4)
        title(titleAverage)
        set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',[0.000,0.700],'ytick',[0.000,0.100,0.200,0.300,0.400,0.500,0.600,0.700])
        ylabel('Response Times, avg across subjects [s]');
        xlabel('Trials before infrequent trial');
        legend('Accurate','Inaccurate');
        legend boxoff
        
        subplot(1,2,2);
        set(gca,'fontsize',24);
        hold on;
        plot(shifts,rts_before_acc,'b','linewidth',1)
        plot(shifts,rts_before_inacc,'m','linewidth',1)
        plot(shifts,nanmean(rts_before_acc),'b','linewidth',4)
        plot(shifts,nanmean(rts_before_inacc),'m','linewidth',4)
        set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',[0.000,0.700],'ytick',[0.000,0.100,0.200,0.300,0.400,0.500,0.600,0.700])
        title(titleSpaghetti)
        xlabel('Trials before infrequent trial');
        
        hold off;
        
    end
    
end
% @todo: write in the todo for this todo
if flag_plotOther == true
    
    figure('Position', [10 10 450 600]);
    
    %subplot(1,2,1);
    set(gca,'fontsize',24);
    hold on;
    errorbar(1,nanmean(std_before_acc),nanstd(std_before_acc)/sqrt(numSubjects),'b','linewidth',4)
    errorbar(1,nanmean(std_before_inacc),nanstd(std_before_inacc)/sqrt(numSubjects),'m','linewidth',4)
    title(titleAverage)
    set(gca,'xtick',1,'xlim',[0.5,1.5],'ylim',[0.000,0.200],'ytick',[0.000,0.100,0.200])
    ylabel('RT StDs, avg across subjects [s]');
    xlabel(['StD of Preceding Trials (n=' num2str(stdWindow) ')']);
    legend('Accurate','Inaccurate');
    legend boxoff
    
    hold off;
    
end


%% Summary Statistics

% @todo: replicate all of the numbers in deBettencourt et al, 2018

chance = 0.5;

A_Prime_attention = A_Prime( att_hit, att_falseAlarm, chance );
[h,p,ci,stats] = ttest2(std_before_acc,std_before_inacc);

% @todo: add in bootstrapping for attention stats
fprintf(['--- NON-BOOTSTRAPPED RESULTS [ATTENTION] ---', '\n']);
fprintf(['CATEGORY JUDGMENT SENSITIVITY: ', ...
    num2str(nanmean(A_Prime_attention)), '\n']);
fprintf(['ERROR RATE - INFREQUENT: ', ...
    num2str(mean(att_falseAlarm)), '\n']);
fprintf(['ERROR RATE - FREQUENT: ', ...
    num2str(mean(att_miss)), '\n']);
fprintf(['OVERALL AVERAGE RT BEFORE CORRECT RESPONSE: ', ...
    num2str(mean(nanmean(rts_before_acc))), ' seconds', '\n']);
fprintf(['OVERALL AVERAGE RT BEFORE INCORRECT RESPONSE: ', ...
    num2str(mean(nanmean(rts_before_inacc))), ' seconds', '\n']);
fprintf(['OVERALL AVERAGE RT VARIANCE BEFORE CORRECT RESPONSE: ', ...
    num2str(mean(nanmean(var_before_acc))), '\n']);
fprintf(['OVERALL AVERAGE RT VARIANCE BEFORE INCORRECT RESPONSE: ', ...
    num2str(mean(nanmean(var_before_inacc))), '\n']);
fprintf(['OVERALL AVERAGE RT STANDARD DEVIATION BEFORE CORRECT RESPONSE: ', ...
    num2str(mean(std_before_acc)), '\n']);
fprintf(['OVERALL AVERAGE RT STANDARD DEVIATION BEFORE INCORRECT RESPONSE: ', ...
    num2str(mean(std_before_inacc)), '\n']);
fprintf(['SEPARATE TRIAL ANALYSES, RT DIFFERENCE BETWEEN CORRECT AND INCORRECT, index_distance=[3,2,1]: ', ...
    num2str(nanmean(rts_before_acc) - nanmean(rts_before_inacc)), ' (seconds)', '\n']); % intuitively, should be all positive values because accurate is slower
fprintf(['SEPARATE TRIAL ANALYSES, RT VARIANCE DIFFERENCE BETWEEN CORRECT AND INCORRECT, index_distance=[3,2,1]: ', ...
    num2str(nanmean(var_before_acc) - nanmean(var_before_inacc)), '\n']); % intuitively, should be all negative values because accurate is lower variability
fprintf(['STATISTICAL ANALYSES, RT STANDARD DEVIATION DIFFERENCE BETWEEN CORRECT AND INCORRECT, subject-level: ', ...
    'p-value = ', num2str(p), '\n']);

% @todo: add in non-bootstrapping for memory stats
% @todo: add in fprintf() to display the results
%fprintf(['--- BOOTSTRAPPED RESULTS [MEMORY] ---', '\n']);
mean_coefficient_rts = mean(logistic_coefficients_rts);
mean_coefficient_var = mean(logistic_coefficients_var);
mean_coefficient_std = mean(logistic_coefficients_std);
bootci(numSamples_bootstrap , @mean , logistic_coefficients_rts);
bootci(numSamples_bootstrap , @mean , logistic_coefficients_var);
bootci(numSamples_bootstrap , @mean , logistic_coefficients_std);

mean_coefficient_rts_var = mean(logistic_c1);
mean_coefficient_rts_var_interaction = mean(logistic_c2);
bootci(numSamples_bootstrap , @mean , logistic_c1);
bootci(numSamples_bootstrap , @mean , logistic_c2);

% @todo: change first argument from mem_hit to mem_remembered
A_Prime_memoryFrequent   = A_Prime( memFreq_hit   , memFreq_falseAlarm   , chance );
A_Prime_memoryInfrequent = A_Prime( memInfreq_hit , memInfreq_falseAlarm , chance );


%% Unfinished/Extra Things

%%%% nonparametric statistics %%%%
% (BS stands for bootstrapped)

% A_Prime_attention, again
BS_att_hit = att_hit(subjectID_list_bootstrap);
BS_att_falseAlarm = att_falseAlarm(subjectID_list_bootstrap);

BS_A_Prime_attention = A_Prime(BS_att_hit, BS_att_falseAlarm, chance);
BS_A_Prime_attention = nanmean(BS_A_Prime_attention,2); % yields 10000 independent random samples

% some copy+pasted code
p_A_Prime_attention = sum( BS_A_Prime_attention( BS_A_Prime_attention<chance) ) / numSamples_bootstrap;
if p_A_Prime_attention == 0, p_A_Prime_attention=1./numSamples_bootstrap; end % not sure what this is for, just copied it from the original codebase

% see if our alternative metric for RT variability is any good for telling apart acc from inacc trials on the basis of RT
% (this gets all "significant" p-values from each subject)
ttest2_pValue_std_AccVsInacc(ttest2_pValue_std_AccVsInacc<0.05)

%%%% proof of concept that direction of means matters %%%%
prevRTsInfreq_ShiftMean_Cor = [];
prevRTsInfreq_ShiftMean_Err = [];
for i = 1:numSubjects
    prevRTsInfreq_ShiftMean_Cor = cat(2, prevRTsInfreq_ShiftMean_Cor, nanmean(comp_rts_before_acc{i}));
    prevRTsInfreq_ShiftMean_Err = cat(2, prevRTsInfreq_ShiftMean_Err, nanmean(comp_rts_before_inacc{i}));
end
% compare prevRTsInfreq_ShiftMean_Cor and prevRTsInfreq_ShiftMean_Err against rts_before_acc and rts_before_inacc

%%%% correlation between RTs and RT_Vars %%%%
corr_coefficients = [];
for i = 1:length(rt_rtVar_corr_collector)
    corr_coefficients = cat(2, corr_coefficients, rt_rtVar_corr_collector{i}.p1);
end
corr_coefficients_mean = mean(corr_coefficients)
bootci(numSamples_bootstrap, @mean, corr_coefficients)

