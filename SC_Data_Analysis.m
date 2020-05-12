%% Introduction

% This is the "master script" for the data analysis of a replication study.
% Authors: Megan deBettencourt, Monica Rosenberg, Steven Cao
% Required Scripts: {inpaint_nans.m, shadedErrorBar.m, dirP.m, A_Prime.m, bootp.m, cohenD.m, bootstrapper.m}
% Other Requirements: MATLAB R2018a (tested up to R2019b), and maybe also the full MATLAB suite

% Additional Notes:
% Search '% #' to jump between various, section IDs (helps for quickly navigating through parts of the code).
% Search '%@' for all hardcoded variables that require manual specification.
% Search '% @' for all "todo" comments within in the code.
% Search '%***' for all parts of the code that are reserved for debugging (they are likely commented-out).
% 
% Bootstrapped confidence intervals are manually calculated via prctile(), which may result in a small discrepancy from using bootci() directly.
% When applying nanmean() to a matrix, the direction matters (i.e. non-commutative), which may sometimes matter when doing exact comparisons of results.
% 
% To adapt this to the original dataset (deBettencourt et al, 2018), make the following changes in this script:
% - The variable name 'attndata' should be 'blockdata' (case sensitive)
% - (And 'attnData' should be 'blockData')
% - For regular (non-realtime) data, exclude ID=60 instead of ID=4,10,12,21,25
% - For realtime data, exclude ID=26 instead of ID=(to be determined)

% @ Big Todo List:
% - Was recently given a heavy rewrite; need to validate all outputs
% - Confirm that when calculating the "overall average" of some metric for preceding trials, that it's calculating means in the correct order
% - (Just a small graphing bug) Scaling of vtc_smooth is fairly displaced; as far as I know, the original figure required some hand-adjusting
% - Make a readme.md file
% - Finish filling out the big todo list

%% #Table of Contents:
% #WS  : Workspace Setup                                Sets up the workspace, e.g. the project directory.
% #SSI : Set Subject IDs                                Fetches all the correct subject IDs.
% #PDM : Preallocate Data Matrices                      Preallocate matrices for all data of all subjects
% #SF  : Set Flags                                      Sets all flags. Also implements meta-flag logic (cf. master flags).
% 
% #ISDa  : Initialise Subject Data (attention)          Get subject's data to later preprocess and analyse.
% #PPSDa : Preprocess Subject Data (attention)          Perform (appropriate) transformations of subject's data so we can do useful analyses on it.
% #PSDa  : Process Subject Data (attention)             Structure the data according to how we will want to analyse it (e.g. get all windows preceding target trials).
% #ISDm  : Initialise Subject Data (memory)             Get subject's data to later preprocess and analyse.
% #PPSDm : Preprocess Subject Data (memory)             Perform (appropriate) transformations of subject's data so we can do useful analyses on it.
% #PSDm  : Process Subject Data (memory)                Structure the data according to how we will want to analyse it (e.g. get all windows preceding target trials).
% 
% #QBA : Quartile-Based Analysis                        (Rosenberg et al, 2013, Figure 3a-3d) Looks at vigilance decrement effects over the session.
% #ASC : Attentional-State Comparison                   (Rosenberg et al, 2013, Figure 2b-2c) Compares error rates of "In the Zone" vs. "Out of the Zone".
% #AMC : Attention Metrics Correlation                  Correlates various attention metrics with each other to indicate how strongly they overlap as measures of attention.
% #AMPTP : Attention Metric Predicts Task Performance   (Numerous articles, ex. deBettencourt et al, 2018, Figure 2a) Looks at whether a measure of attention can reliably predict outcomes of responses to infrequent ("target") trials
% #AMPMP : Attention Metric Predicts Memory Performance (deBettencourt et al, 2018, Figure 2b-2c) Looks at whether a measure of attention can reliably predict whether some stimulus/trial will later be well-remembered.
% #PLOTS : Plot all the things                          In actuality, where most - not all - of the plotting code is located.
% #STATS : Summary statistics                           Where most of the meaningful output for the statistical analyses is located.
% #STATSA  : Attention analyses
% #STATSAB : Attention analyses, bootstrapped
% #STATSM  : Memory analyses
% #STATSMB : Memory analyses, bootstrapped
% 
% #SAVE : Save outputs
% #OTHER : Miscellaneous and/or unlabelled sections
% 


%% #WS - Script Setup

% Flush session history.
clear; clc; close all;

% Set seed for bootstrapping (the value of 999 is arbitrary).
stream = RandStream('mt19937ar','Seed',999);
RandStream.setGlobalStream( stream );

% Get main directory (it is assumed this script is running from the study's main directory).
projectDirectory = mfilename('fullpath');
% (Fragment into foldernames/filenames, get last item which is the filename of this script, subtract that, then put it all back together to get the script's home folder.)
projectDirectory = split(projectDirectory, {'/','\'}); projectDirectory = projectDirectory(1:end-1); projectDirectory = strjoin(projectDirectory, '/');
%@ Alternatively can specify absolute path manually for debugging purposes.
%projectDirectory = '/Users/rosenberglab/Desktop/discCPT';

% Also define the project's subfolders accordingly, for convenience.
% (Subject data folders to be parsed through)
folderName_data   = 'data';
folderName_rtdata = 'rtdata';
projectDirectory_data   = fullfile( projectDirectory, folderName_data  , '/' );
projectDirectory_rtdata = fullfile( projectDirectory, folderName_rtdata, '/' );
% (Folder for storing outputted figures)
folderName_plots = 'figures';
projectDirectory_plots = fullfile( projectDirectory, folderName_plots  , '/' );
    if ~exist(projectDirectory_plots, 'dir')
        mkdir(projectDirectory_plots);
        fprintf(['FOLDER FOR FIGURES NOT FOUND', '\n\n\n', 'Creating one now...', '\n\n']);
    end
% (Folder for storing outputted stats)
folderName_stats = 'mat';
projectDirectory_stats = fullfile( projectDirectory, folderName_stats  , '/' );
    if ~exist(projectDirectory_stats, 'dir')
        mkdir(projectDirectory_stats);
        fprintf(['FOLDER FOR ANALYSIS OUTPUTS NOT FOUND', '\n\n\n', 'Creating one now...', '\n\n']);
    end


%% #SSI - Preliminary Variable Setup

%%% Script Parameters
% (All subject data folders should be integer IDs)
shifts = -3:1:-1; %@ hardcoded
numShifts = numel(shifts);
stdWindow = 9; %@ hardcoded; keep in mind that the actual window size will always be +1, because we also include the RT of the target trial
%@ hardcoded; these parameters are determined by the settings of the experiment - here they are set to the default
numAttTrialsTotal = 500;
proportionFreq   = 0.90;    numTrialsFreq   = numAttTrialsTotal * proportionFreq;
proportionInfreq = 0.10;    numTrialsInfreq = numAttTrialsTotal * proportionInfreq;
numMemTrialsTotal = 200;
%@ hardcoded; see 'rememberThreshold' below under the analysis flags
rememberRateThreshold = 0.10;

%%% Subject ID Selection
% (Robustly generates a list of all available subject IDs from the appropriate data folder.)
realtime = false; %@ analyse experimental setup with no trial triggering or with trial triggering?
if realtime == true
    assert( ~isempty(isnumeric(str2num(char(dirP(projectDirectory_rtdata)')))) , 'WARNING: Invalid subject data folder' ); %#ok<ST2NM>
    subjectID_list = sort(str2num(char(dirP(projectDirectory_rtdata)'))'); %#ok<ST2NM,TRSRT>
elseif realtime == false
    assert( ~isempty(isnumeric(str2num(char(dirP(projectDirectory_data)')))) , 'WARNING: Invalid subject data folder' ); %#ok<ST2NM>
    subjectID_list = sort(str2num(char(dirP(projectDirectory_data)'))'); %#ok<ST2NM,TRSRT>
end

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

%%% #PDM - Preinitialise Data Collectors

%%%% (Attention Accuracy)
att_hit                 = nan(1,numSubjects); % correctly identifying frequent trials
att_falseAlarm          = nan(1,numSubjects); % mistaking infrequent trials (e.g. labeling them as frequent)
att_miss                = nan(1,numSubjects); % mistaking frequent trials (e.g. labeling them as infrequent)
att_correctReject       = nan(1,numSubjects); % correctly identifying infrequent trials (e.g. successfully overriding habitual keypress)

%%%% (Attention Accuracy - #ASC Initialisation)
att_falseAlarm_ITZ      = nan(1,numSubjects); % out of all the infrequent trials, how many were mistaken AND "in the zone"?
att_falseAlarm_OTZ      = nan(1,numSubjects); % out of all the infrequent trials, how many were mistaken AND "out of the zone"?
att_miss_ITZ            = nan(1,numSubjects); % out of all the frequent trials, how many were mistaken AND "in the zone"?
att_miss_OTZ            = nan(1,numSubjects); % out of all the frequent trials, how many were mistaken AND "out of the zone"?

%%%% (Attention Other)
std_accVsInacc_ttest2_pValue = nan(numSubjects,1); % for each subject, are the standard deviations of the windows of preceding trials significantly different between accurate and inaccurate target trials? (i.e. is this a reliable metric for attention?)

%%%% (Memory Accuracy)
memFreq_hit             = nan(1,numSubjects); % correctly remembering a "frequent-category" picture
memInfreq_hit           = nan(1,numSubjects); % correctly remembering an "infrequent-category" picture
memFreq_falseAlarm      = nan(1,numSubjects); % falsely remembering a "frequent-category" picture
memInfreq_falseAlarm    = nan(1,numSubjects); % falsely remembering an "infrequent-category" picture
memFreq_remembered      = nan(1,numSubjects); % correctly and strongly remembering a "frequent-category" picture
memInfreq_remembered    = nan(1,numSubjects); % correctly and strongly remembering an "infrequent-category" picture
memFreq_forgotten       = nan(1,numSubjects); % otherwise failing to meet the criteria for remembered for "frequent-category" pictures
memInfreq_forgotten     = nan(1,numSubjects); % otherwise failing to meet the criteria for remembered for "infrequent-category" pictures
memFreq_cFalseAlarm     = nan(1,numSubjects); % (to ensure consistency with particular A_Prime calculations) for cases where subject confidently remembers a completely new "frequent-category" picture
memInfreq_cFalseAlarm   = nan(1,numSubjects); % (to ensure consistency with particular A_Prime calculations) for cases where subject confidently remembers a completely new "infrequent-category" picture

%%%% (Memory Other)
mem_sessionLength       = nan(1,numSubjects); % the sum of "reaction times" to all (200) trials in the memory task, in seconds

%%%% (Attention Metrics - obtained from those sets of trials which precede an infrequent "target" trial)
% Note: When using nanmean() to compress a matrix into a scalar (2d into 0d), the result WILL depend on the order of "dimension compression"!
    % (Average ACROSS shifts first; this is the calculation used in code block #10 of the original python script)
    % ACROSS shifts == we want to know the average RT (or whatever value) of each iShift; each iShift has only one value, which is an average across multiple trials
rts_before_acc          = nan(numSubjects,numShifts); % typical RTs before a correctly-responded trial
rts_before_inacc        = nan(numSubjects,numShifts); % typical RTs before an incorrectly-responded trial
var_before_acc          = nan(numSubjects,numShifts); % typical RT deviances before a correctly-responded trial
var_before_inacc        = nan(numSubjects,numShifts); % typical RT deviances before an incorrectly-responded trial
    % (Average WITHIN shifts first; this is the calculation used in code block #8 of the original python script)
    % WITHIN shifts == we want to know the average RT (or whatever value) for each target trial; each target trial has only one value, which is an average of all of the target trial's shifts
comp_rts_before_acc     = cell(numSubjects,1); % average RT before each correctly-responded trial
comp_rts_before_inacc   = cell(numSubjects,1); % average RT before each incorrectly-responded trial
comp_var_before_acc     = cell(numSubjects,1); % average RT deviance before each correctly-responded trial
comp_var_before_inacc   = cell(numSubjects,1); % average RT deviance before each incorrectly-responded trial
    % (Shifts are inapplicable for standard deviation calculations.)
std_before_acc          = nan(numSubjects,1); % typical StD of RTs preceding a correctly-responded trial
std_before_inacc        = nan(numSubjects,1); % typical StD of RTs preceding an incorrectly-responded trial

%%%% (Attention Metrics - #QBA Initialisation)
quartile_commission_errors = nan(numSubjects,4);
quartile_ommission_errors  = nan(numSubjects,4);
quartile_rt                = nan(numSubjects,4);
quartile_rtvar             = nan(numSubjects,4);

%%%% (Attention Metrics Comparison - #AMC Initialisation)
rts_var_rho = nan(numSubjects,1);
rts_var_rho_pValue = nan(numSubjects,1);

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
logistic_coefficients_rts_var = nan(numSubjects,3); % constant, first coefficient (RT), second coefficient (RT Var)
logistic_coefficients_rts_var_interact = nan(numSubjects,4); % constant, first coefficient (RT), second coefficient (RT Var), third coefficient (interaction)

%%%% (Recent additions; experimental)
% @todo: make these not-so-experimental
% preceding and frequent trial analyses
i1_logistic_coefficients_rts   = nan(numSubjects,2);
i1_logistic_coefficients_var   = nan(numSubjects,2);
%i1_logistic_coefficients_std   = nan(numSubjects,2);
i1_logistic_coefficients_rts_var = nan(numSubjects,3);
i1_logistic_coefficients_rts_var_interact = nan(numSubjects,4);
f_logistic_coefficients_rts   = nan(numSubjects,2);
f_logistic_coefficients_var   = nan(numSubjects,2);
%f_logistic_coefficients_std   = nan(numSubjects,2);
f_logistic_coefficients_rts_var = nan(numSubjects,3);
f_logistic_coefficients_rts_var_interact = nan(numSubjects,4);
% correlate mean of rt sliding window and var for each subject
mrts_var_corr_collector = nan(numSubjects,2);
% postdiction stuff
postdiction_trigger = nan(numSubjects,numAttTrialsTotal);
postdiction_vtc = nan(numSubjects,numAttTrialsTotal);
postdiction_triggerClassic = nan(numSubjects,numAttTrialsTotal);
postdiction_triggerRandomNum = 1000;
postdiction_triggerRandom = nan(numSubjects,numAttTrialsTotal,postdiction_triggerRandomNum);
% cohen's d stuff
numPermutations = 200;
nullBetas_rts = nan(numSubjects, numPermutations);
nullBetas_var = nan(numSubjects, numPermutations);
nullBetas_RTS_var = nan(numSubjects, numPermutations);
nullBetas_rts_VAR = nan(numSubjects, numPermutations);
% for logistic regressions using frequent trials
comp_rts_before_freq = cell(numSubjects,1);
comp_var_before_freq = cell(numSubjects,1);


%% #SF - Set Flags

% Initialise containers for all of our flags (plotting, analysis, and miscellaneous).
pFlags = containers.Map( 'KeyType', 'char' , 'ValueType', 'logical' );
aFlags = containers.Map( 'KeyType', 'char' , 'ValueType', 'logical' );
mFlags = containers.Map( 'KeyType', 'char' , 'ValueType', 'logical' );

%%% Script Flags - Plotting
% plotOverall_precedingTrials_RTs:              Plot results of preceding trials, using raw RTs (i.e. for comparability to deBettencourt et al, 2018, Figure 2a).
% plotOverall_precedingTrials_RTDevs:           Plot results of preceding trials, using RT Deviations (i.e. relating to RT variability).
% plotOverall_precedingTrials_other:            Plot results of preceding trials, using other metrics (e.g. standard deviation of RT windows preceding infrequent trials).
% plotOverall_logistic_RTs:                     Plot all subjects' relationship between performance in attention and in memory, using raw RTs (i.e. for comparability to deBettencourt et al, 2018, Figure 2c).
% plotOverall_logistic_RTDevs:                  Plot all subjects' relationship between performance in attention and in memory, using RT Deviations (i.e. relating to RT variability).
% plotOverall_logistic_other:                   Plot all subjects' relationship between performance in attention and in memory, using other metrics (e.g. standard deviation).
% plotOverall_compareAttentionStates_errors:    Plot commission (false alarm) and omission (miss) errors for good attention states and for bad attention states (i.e. for comparability to Rosenberg et al, 2013, Figure 2b+2c).
% plotOverall_QBA_vigilanceDecrements:          Plot errors and attention metrics to showcase vigilance decrements over the course of the session (i.e. for comparability to Rosenberg et al, 2013, Figure 3).
% plotSubjects_logistic_RTs:                    Plot each subject's relationship between performance in attention and in memory, using raw RTs.
% plotSubjects_logistic_RTDevs:                 Plot each subject's relationship between performance in attention and in memory, using RT Deviations.
% plotSubjects_logistic_other:                  Plot each subject's relationship between performance in attention and in memory, using other metrics.
% plotSubjects_varianceTimeCourse:              Plot each subject's VTC, along with a smoothed version of their VTC (i.e. for replicating Rosenberg et al, 2013, Figure 2a).

pFlags('plotOverall_precedingTrials_RTs')           = true;
pFlags('plotOverall_precedingTrials_RTDevs')        = true;
pFlags('plotOverall_precedingTrials_other')         = true;
pFlags('plotOverall_logistic_RTs')                  = true;
pFlags('plotOverall_logistic_RTDevs')               = true;
pFlags('plotOverall_logistic_other')                = true;
pFlags('plotOverall_compareAttentionStates_errors') = true;
pFlags('plotOverall_QBA_vigilanceDecrements')       = true;
pFlags('plotSubjects_logistic_RTs')                 = false;
pFlags('plotSubjects_logistic_RTDevs')              = false;
pFlags('plotSubjects_logistic_other')               = false;
pFlags('plotSubjects_varianceTimeCourse')           = false;

%%% Script Flags - Analyses
% detrend:           Remove potential time-based effects on attention and/or response times (e.g. vigilance decrement).
%                    (Set to TRUE for comparability to deBettencourt et al, 2018. But for general purposes: recommended to leave as TRUE to extricate vigilance decrement effects from analyses.)
% includeAllTrials:  Include the RTs of infrequent and/or incorrect trials when collating the RTs of trials which precede the target trial of interest.
%                    (Set to TRUE for comparability to deBettencourt et al, 2018. But for general purposes: recommended to leave as FALSE to extricate effects peculiar to errors, e.g. post-error slowing.)
% interpolateRTs:    Replace trials that do not have a response time with an interpolated value. Note that this also means VTC will be interpolated.
%                    (Set to FALSE for comparability to deBettencourt et al, 2018. But for general purposes: recommended to leave as FALSE; this is mainly a fast workaround for functions that hate NaN values.)
% rememberThreshold: Whether logistic regression should exclude subjects who have too low (or too high) remembrance rate. Note that subjects with perfect separation will still be ignored.
%                    (Set to FALSE for comparability to deBettencourt et al, 2018. But for general purposes: recommended to leave as TRUE for more stable logistic regression coefficients.)

aFlags('detrend')            = true;
aFlags('includeAllTrials')   = false;
aFlags('interpolateRTs')     = false;
aFlags('rememberThreshold')  = true;

%%% Script Flags - Miscellaneous
% MASTER_PLOT:       If true, forces all plotting flags to TRUE; otherwise does nothing. Good for laziness, debugging, and lazy debugging.
%                    (Not recommended to run this with MASTER_TEST on, lest you end up with probably thousands of figures.)
% MASTER_TEST:       Goes through every possible combination of analysis flags. Great for generating supplemental stuff (and debugging).
%                    (Each iteration of data analysis takes approximately 15 seconds.)
% SAVE_PLOTS:        Saves all generated figures as both a .fig and a .png to the 'figures' folder.
%                    (Not recommended unless you know precisely which plots will be generated; see above.)
% SAVE_STATS:        Saves all statistical outputs in the form of a .mat data file to the 'mat' folder.
%                    (This is going to save multiple .mat variants if MASTER_TEST is on.)
% @todo: optimise save_stats implementation; seems to be costing a lot of time (around 20s)

mFlags('MASTER_PLOT')        = false;
mFlags('MASTER_TEST')        = false;
mFlags('SAVE_PLOTS')         = false;
mFlags('SAVE_STATS')         = false;


if mFlags('MASTER_PLOT') == true
    
    % Set all flags in pFlags to true.
    pFlags = containers.Map( pFlags.keys , true(1,pFlags.keys) );
    
end


list_flags = aFlags.keys; % get the list of all analysis flags
num_flags = length(list_flags);

if mFlags('MASTER_TEST') == true
    
    % Generate all possible combinations of flags.
    num_combinations = 2 ^ num_flags;
    flag_combinations = sortrows( de2bi( 0:(num_combinations-1) ) , num_flags );
    
elseif mFlags('MASTER_TEST') == false
    
    % Otherwise, generate only one possible combination of flags, which would be the one we specified.
    num_combinations = 1;
    flag_combinations = cell2mat(aFlags.values);
    
end


%% #ISD;#PPSD;#PSD - Subject Analysis
for iFlagConfiguration = 1:num_combinations
    saveSetupPrefix = ''; %#ok<*AGROW> % (this variable is used solely for adding flag-related metadata to saved outputs)
    
    % If mFlags('MASTER_TEST') was true, then this will actively set the analysis-related flags in order to cycle through all their possible combinations.
    % Otherwise, it will still "set" the flags again, but only to the same values that we manually specified above. So it will only echo what we set the flags to be.
    current_flag_combination = flag_combinations(iFlagConfiguration,:);
    
    for i = 1:num_flags
        flagName = list_flags{i}; % get the current analysis flag's name...
        aFlags(flagName) = current_flag_combination(i); % ...and set that flag to its appropriate boolean value
        
        if     aFlags(flagName) == true && flagName == "detrend"
            saveSetupPrefix = [saveSetupPrefix, 'Detrended'];
            
        elseif aFlags(flagName) == true && flagName == "includeAllTrials"
            saveSetupPrefix = [saveSetupPrefix, 'IncludedAll'];
            
        elseif aFlags(flagName) == true && flagName == "interpolateRTs"
            saveSetupPrefix = [saveSetupPrefix, 'InterpolatedRTs'];
        
        elseif aFlags(flagName) == true && flagName == "rememberThreshold"
            saveSetupPrefix = [saveSetupPrefix, 'RemRateThreshold'];
            
        end %@ add more as needed (i.e. if we add in more analysis flags)
        
        
    end
    
    saveSetupPrefix = [saveSetupPrefix, '_'];
    
    % Open up a text file to log output to if 'SAVE_STATS' is toggled.
    if mFlags('SAVE_STATS') == true
        outputLog_name = [saveSetupPrefix, 'SustainedAttentionExperiment_Summary.txt'];
        outputLog = fopen( [projectDirectory_stats, outputLog_name] , 'w' );
    else
        outputLog = 1; % i.e. log to the command window
    end
    
    % Print the currently-toggled analysis flags for this iteration.
    fprintf(outputLog, ['\n', 'SETTINGS: ', saveSetupPrefix, '\n']);
    
for iSubject = 1:length(subjectID_list)
    % Notes on conventions:
    % iSubject is an index for filling subject-related data variables. This is always in order, e.g. [1, 2, 3, ...]
    % subjectID is the actual subject's ID. This is not always in order, e.g. [51, 52, 53, 59, ...]
    % Accuracy is always encoded as either 0 or 1, and not some other value (e.g. NaN)
    
    % Retrieve subject's data folder.
    subjectID = subjectID_list(iSubject);
    if realtime == true
        subjectID_dir = fullfile(projectDirectory_rtdata, num2str(subjectID), '/');
    elseif realtime == false
        subjectID_dir = fullfile(projectDirectory_data  , num2str(subjectID), '/');
    end
    
    % Load subject's attention data.
        file = dir([subjectID_dir, 'attndata_*']);
        assert(numel(file)<=1,['ERROR: More than one attention datafile for subject ID: ' num2str(subjectID)])
        load([subjectID_dir, file.name]); % this will load in the variable called 'attnData'
    % Load subject's memory data.
        file = dir([subjectID_dir, 'memdata_*']);
        assert(numel(file)<=1,['ERROR: More than one memory datafile for subject ID: ' num2str(subjectID)])
        load([subjectID_dir, file.name]); % this will load in the variable called 'memData'
    
    %%% Attention Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % #ISDa
    % Pull relevant data to be analysed.
    att_accs = attnData.accs;
    att_categs = attnData.categs;
    att_RTs = attnData.rts;
    att_trialIDs = attnData.trial;
    
    
    % #PPSDa
    % Replace trials with a reaction time of below 100ms with NaN values (since that is physiologically impossible).
    att_RTs( att_RTs<0.1000 ) = NaN;
    
    % If we want detrended RTs instead (to account for natural attention decrements over time/trials), then
    % convert RTs into RT-residuals (THIS MEANS 'att_RTs' WILL NOW REFER TO THE RESIDUALS AND NOT TO THE RAW RT VALUES).
    if aFlags('detrend') == true
        % Using detrend() directly is not recommended because it doesn't work with NaN values.
        
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
    if aFlags('includeAllTrials') == false
        
        validShift_RTs = att_RTs;
        validShift_RTs( att_accs==0 | att_categs==INFREQUENT ) = NaN;
        
    % If we do want to include all trials, then all RTs are considered "valid": none of them will be replaced with NaN.
    elseif aFlags('includeAllTrials') == true
        
        validShift_RTs = att_RTs;
        
    end
    
    % Calculate the variance time course (VTC) for an alternative measure (the main measure being the raw RTs).
    % (Method of calculating VTC is derived from Rosenberg et al, 2013.)
    rt_mean   = nanmean(validShift_RTs);
    rt_std    = nanstd(validShift_RTs);
    
    % Interpolate RT (and consequently, VTC) values.
    % Note: this affects all data analysis. Leave this on 'false' to use interpolation only for calculating smoothed VTC.
    if aFlags('interpolateRTs') == true
        
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
    vtc_zoneSplit = median(vtc_smooth); % for getting the median split of the VTC in order to determine which values correspond to ITZ and OTZ
    
    
    % #PSDa
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
    
    % #PSDa
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
    
    % Test for each subject as to whether there is a significant difference in RT std values between the accurate and the
    % inaccurate trials; this will establish whether this alternative metric can be used as a measure of attention by testing
    % its ability to distinguish between accurately-responded trials and inaccurately-responded trials.
    % Note that there will be some nan values that are removed first - these are the trials which we couldn't get an stdWindow for.
    [~,std_accVsInacc_ttest2_pValue(iSubject),~,~] = ...
        ttest2( rmmissing(all_std_before_acc) , rmmissing(all_std_before_inacc) );
    
    % Perform postdiction (used only to tweak realtime triggering setup).
    % Note that the values are already detrended and all <100ms values removed.
    numPadTrials = 50;
    [ postdiction_trigger(iSubject,:) , postdiction_vtc(iSubject,:) , postdiction_triggerClassic(iSubject,:) ] = ...
        postdiction( attnData.rts, att_categs, att_accs, FREQUENT, numShifts, numAttTrialsTotal, numPadTrials);
    % issue: linear drift was not based on expanding window; temporary workaround in place(feed raw rts instead)
    % @todo: unscaffold
        trg_vtc(iSubject,:) = [sum(postdiction_trigger(iSubject,:)==-1) , sum(postdiction_trigger(iSubject,:)==1)];
        trg_rts(iSubject,:) = [sum(postdiction_triggerClassic(iSubject,:)==-1) , sum(postdiction_triggerClassic(iSubject,:)==1)];
        [postVsWhole_vtc_r(iSubject), postVsWhole_vtc_p(iSubject)] = corr( postdiction_vtc(iSubject,:)' , vtc' , 'rows','complete' );
    if false
        % trigger counts
        figure;
        hist(postdiction_trigger(iSubject,:));
        % vtc stuff
        figure;
        plot(postdiction_vtc(iSubject,:));
        hold on;
        % treat nans as 0s, due to the way find() works
        [~,asdf_tf] = rmmissing(postdiction_trigger(iSubject,:));
        asdf_nonans = postdiction_trigger(iSubject,:);
        asdf_nonans(asdf_tf) = 0;
        asdf_find = find(asdf_nonans);
%        for iline=1:length(asdf_find)
%            xline( asdf_find(iline) );
%        end
        hold off;
        % correlate simulated vtc with post-session vtc
    end
    % also generate randomly-triggered trials (with same proportion as subject's simulated triggered trials)
    % to have an idea of what our chance/baseline is; note that this only generates random triggers, i.e. no direction
    for iPermutation = 1:postdiction_triggerRandomNum
        num_triggerTrials = sum(  abs(postdiction_trigger(iSubject,:))==1  );
        num_nontriggerTrials = numAttTrialsTotal - (numPadTrials*2) - num_triggerTrials;
        unshuffled_triggers = [ ones(1,num_triggerTrials) , zeros(1,num_nontriggerTrials) ];
        shuffled_triggers = unshuffled_triggers(randperm(length(unshuffled_triggers)));
        postdiction_triggerRandom(iSubject,:,iPermutation) = [nan(1,numPadTrials), shuffled_triggers, nan(1,numPadTrials)];
    end
    
    
    %%% Memory Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % #ISDm
    % Pull relevant data to be analysed.
    mem_accs = memData.accs;
    mem_categs = memData.memCategs; % 1 == frequent+old, 2 == infrequent+old, 3 == frequent+new, 4 == infrequent+new
    mem_ratings = memData.rating; % 1 == definitely new, 2 == maybe new, 3 == maybe old, 4 == definitely old
    mem_trialIDs = memData.trial;
    % Find out how long the session took (in case we're curious).
    % (The +1.00 is because the picture remains for 500ms post-response, followed by a 500ms interstimulus interval.)
    mem_sessionLength(iSubject) = memData.actualTrialOffsets(end) - memData.actualTrialOnsets(1) + 1.00;
    
    
    % #PPSDm
    memFreq_hit         (iSubject) = mean(mem_accs(mem_categs==1)==1);
    memInfreq_hit       (iSubject) = mean(mem_accs(mem_categs==2)==1);
    memFreq_falseAlarm  (iSubject) = mean(mem_accs(mem_categs==3)==0);
    memInfreq_falseAlarm(iSubject) = mean(mem_accs(mem_categs==4)==0);
    % Trials are only considered "remembered" if the responses were also confident:
    memFreq_remembered  (iSubject) = sum( mem_categs==1 & mem_ratings==4 ) / length(mem_accs(mem_categs==1));
    memInfreq_remembered(iSubject) = sum( mem_categs==2 & mem_ratings==4 ) / length(mem_accs(mem_categs==2));
    memFreq_forgotten   (iSubject) = 1.0 - memFreq_remembered  (iSubject);
    memInfreq_forgotten (iSubject) = 1.0 - memInfreq_remembered(iSubject);
    % (Used to calculate A_Prime for memFreq and memInfreq to ensure consistency with calculations used in deBettencourt et al, 2018)
    memFreq_cFalseAlarm  (iSubject) = sum( mem_categs==3 & mem_ratings==4 ) / length(mem_accs(mem_categs==1));
    memInfreq_cFalseAlarm(iSubject) = sum( mem_categs==4 & mem_ratings==4 ) / length(mem_accs(mem_categs==2));
    
    % Get only the infrequent+old trial IDs from the attention task.
    % (Named 'both' because it's both acc and inacc.)
    trial_infreq_both = memData.attOrder(mem_categs==2);
    
    
    % #PSDm
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
    
    % #PSDm
    % Calculate some summary statistics for RTs on the current subject ID.
    rts_before_infreq(iSubject,:) = nanmean(all_rts_before_infreq);   
    comp_rts_before_infreq{iSubject} = nanmean(all_rts_before_infreq,2);
    
    var_before_infreq(iSubject,:) = nanmean(all_var_before_infreq);
    comp_var_before_infreq{iSubject} = nanmean(all_var_before_infreq,2);
    
    std_before_infreq(iSubject) = nanmean(all_std_before_infreq);
    std_before_infreq_all{iSubject} = all_std_before_infreq; % will need this for comparing predicted values against actual values wit logistic regression
    
    
    % #PPSDm
    % Arrange vectors in the correct order for use in logistic regression modeling.
    % We need an inputVariable vector (e.g. RTs), an outputVariable vector (i.e. accuracy), and a vector of ones (because MATLAB is picky).
    % Additionally, all three vectors should have the same length (numInfreqTrials==50) and be column vectors.
    
    % Output: Accuracy Vector (boolean vector)
    accuracyVector = mem_ratings(mem_categs==2); % get only the trials of interest (old+infrequent)
    accuracyVector(accuracyVector~=4)=0; accuracyVector(accuracyVector==4)=1; % only confident responses (4) are considered correct
    accuracyVector=accuracyVector'; % transpose from row to column
    % If 'rememberThreshold' is true, then also check that there is a "good enough mix" of 1s and 0s for the logistic regression outputs to yield feasible values.
    % (Otherwise the beta coefficients could be incredibly high, skewing the results.)
    % If 'rememberThreshold' is false, then at least check that there is not perfect separation or some other problem (i.e. accuracyVector should have only 2 unique values).
    % If there is a threshold and we fail it (i.e. close to perfect separation), or there isn't and we have perfect separation, then all regression-related values will default to NaN for that subject.
    if ( (aFlags('rememberThreshold') == true) && ...
         ((memInfreq_remembered(iSubject) <= 1.0-rememberRateThreshold) && (memInfreq_remembered(iSubject) >= rememberRateThreshold)) ) ...
    ...
    || ( (aFlags('rememberThreshold') == false) && ...
         (length(unique(accuracyVector)) == 2) )
    
    % Input: RT Vector (numeric vector)
    inputVector_rts = comp_rts_before_infreq{iSubject};
    
    % Input: VTC Vector (numeric vector)
    inputVector_var = comp_var_before_infreq{iSubject};
    
    % Input: RT StD Vector (numeric vector)
    inputVector_std = std_before_infreq_all{iSubject};
    
    % Note: All the vectors above already correspond to the correct (i.e. the same) trial IDs, so we can pair the vectors off just fine.
    
    % Auxiliary: "Number of times each value of the input variable is repeated" Vector (numeric vector)
    testedVector = ones(length(accuracyVector),1);
    
    
    % #PSDm
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
    
    % #PSDm
    % Append model results.
    logistic_coefficients_rts(iSubject,:) = logModel_rts_coefficients;
    rts_before_infreq_estimated{iSubject} = logModel_rts_predicted;
    
    logistic_coefficients_var(iSubject,:) = logModel_var_coefficients;
    var_before_infreq_estimated{iSubject} = logModel_var_predicted;
    
    logistic_coefficients_std(iSubject,:) = logModel_std_coefficients;
    std_before_infreq_estimated{iSubject} = logModel_std_predicted;
    
    
    % #AMPMP
    % Plot logistic regression on raw RTs for each individual.
    if pFlags('plotSubjects_logistic_RTs') == true
        figure('Name',[saveSetupPrefix,'logisticRT_ID',num2str(subjectID)],'Position', [10 10 900 600]);
        hold on;
        plot(logModel_rts_input, logModel_rts_output, 'bs', logModel_rts_input, logModel_rts_predicted, 'r-', 'LineWidth', 2)
        set(gca,'fontsize',20);
        title("RT Attention & Memory Performance")
        xlabel("RT"); ylabel("Memory Performance");
        set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);
    end
    
    % Plot logistic regression on var for each individual.
    if pFlags('plotSubjects_logistic_RTDevs') == true
        figure('Name',[saveSetupPrefix,'logisticVar_ID',num2str(subjectID)],'Position', [10 10 900 600]);
        hold on;
        plot(logModel_var_input, logModel_var_output, 'bs', logModel_var_input, logModel_var_predicted, 'r-', 'LineWidth', 2)
        set(gca,'fontsize',20);
        title("RT Var Attention & Memory Performance")
        xlabel("RT Var"); ylabel("Memory Performance");
        set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);
    end
    
    % Plot logistic regression on std for each individual.
    if pFlags('plotSubjects_logistic_other') == true
        figure('Name',[saveSetupPrefix,'logisticSTD_ID',num2str(subjectID)],'Position', [10 10 900 600]);
        hold on;
        plot(logModel_std_input, logModel_std_output, 'bs', logModel_std_input, logModel_std_predicted, 'r-', 'LineWidth', 2)
        set(gca,'fontsize',20);
        title("RT Std Attention & Memory Performance")
        xlabel("RT Std"); ylabel("Memory Performance");
        set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);
    end
    
    
    %% #PSDm - Memory Analysis (Extra)
    
    % Logistic regression, RTs+Var & RTs+Var+Interaction
    % Resort in ascending order for the inputVector (RTs), so we can get the model's estimated values in order (in case we want to plot for some reason).
    resorted = sortrows([inputVector_rts, inputVector_var, accuracyVector, testedVector]);
    % Delete invalid trials (i.e. rows which have NaNs, i.e. rows which have no input value).
    resorted = rmmissing(resorted);
    % Generate arguments for glmfit() and glmval().
    logModel_combined_input_rts = resorted(:,1);
    logModel_combined_input_var = resorted(:,2);
        % (Standardise units to feed into the model with an interaction effect.)
        logModel_combined_input_rts = ( logModel_combined_input_rts - nanmean(logModel_combined_input_rts) ) / nanstd(logModel_combined_input_rts);
        logModel_combined_input_var = ( logModel_combined_input_var - nanmean(logModel_combined_input_var) ) / nanstd(logModel_combined_input_var);
    logModel_combined_output    = resorted(:,3);
    logModel_combined_aux       = resorted(:,4);
    % Generate both logistic regression models.
    logModel_rts_var_coefficients = glmfit( [logModel_combined_input_rts, logModel_combined_input_var], ...
        [logModel_combined_output, logModel_combined_aux], 'binomial', 'link', 'logit');
    logModel_rts_var_interact_coefficients = glmfit( [logModel_combined_input_rts, logModel_combined_input_var, logModel_combined_input_rts .* logModel_combined_input_var], ...
        [logModel_combined_output, logModel_combined_aux], 'binomial', 'link', 'logit');
    % @todo: generate predicted values, in case we want to plot or something
    
    % Append model results.
    logistic_coefficients_rts_var(iSubject,:) = logModel_rts_var_coefficients;
    logistic_coefficients_rts_var_interact(iSubject,:) = logModel_rts_var_interact_coefficients;
    
    
    % Generate null betas for each logistic model type (RTs, Var, RTs+Var).
    for iPermutation = 1:numPermutations
        accuracyVector_shuffled = accuracyVector(randperm(length(accuracyVector)));
        
        
        % Response Times
        
        p_resorted = sortrows( [inputVector_rts, accuracyVector_shuffled, testedVector] );
        % Delete invalid trials (i.e. rows which have NaNs, i.e. rows which have no input value)
        p_resorted = rmmissing(p_resorted);
        % Generate arguments for glmfit()
        p_logModel_rts_input  = p_resorted(:,1);
        p_logModel_rts_output = p_resorted(:,2);
        p_logModel_rts_aux    = p_resorted(:,3);
        % Generate the logistic regression model
        p_logModel_rts_coefficients = glmfit( p_logModel_rts_input, [p_logModel_rts_output, p_logModel_rts_aux], 'binomial', 'link', 'logit');
        
        % Append the newly-generated null beta
        nullBetas_rts(iSubject, iPermutation) = p_logModel_rts_coefficients(2);
        
        
        % Response Time Variance
        
        p_resorted = sortrows( [inputVector_var, accuracyVector_shuffled, testedVector] );
        % Delete invalid trials (i.e. rows which have NaNs, i.e. rows which have no input value)
        p_resorted = rmmissing(p_resorted);
        % Generate arguments for glmfit()
        p_logModel_var_input  = p_resorted(:,1);
        p_logModel_var_output = p_resorted(:,2);
        p_logModel_var_aux    = p_resorted(:,3);
        % Generate the logistic regression model
        p_logModel_var_coefficients = glmfit( p_logModel_var_input, [p_logModel_var_output, p_logModel_var_aux], 'binomial', 'link', 'logit');
        
        % Append the newly-generated null beta
        nullBetas_var(iSubject, iPermutation) = p_logModel_var_coefficients(2);
        
        
        % Response Times + Response Time Variance
        % @todo: this section is currently spitting out warnings (also, the results using these figures seem kind of screwed up)
        
        p_resorted = sortrows( [inputVector_rts, inputVector_var, accuracyVector_shuffled, testedVector] );
        % Delete invalid trials (i.e. rows which have NaNs, i.e. rows which have no input value)
        p_resorted = rmmissing(p_resorted);
        % Generate arguments for glmfit()
        p_logModel_combined_input_rts = p_resorted(:,1);
        p_logModel_combined_input_var = p_resorted(:,2);
        p_logModel_combined_output    = p_resorted(:,3);
        p_logModel_combined_aux       = p_resorted(:,4);
        % Generate the logistic regression model (along with some extra error handling things)
    lastwarn('');
        p_logModel_combined_coefficients = glmfit( [p_logModel_combined_input_rts, p_logModel_combined_input_var], ...
            [p_logModel_combined_output, p_logModel_combined_aux], 'binomial', 'link', 'logit');
    [~,lastWarningID] = lastwarn;
    % If, due to sparsity in the accuracy (output) vector, we can't converge properly, then keep reshuffling it until we can.
    % @todo: implement this workaround for all the other logistic regressions?
    while strcmp(lastWarningID, 'stats:glmfit:IterationLimit')
        fprintf(['Working around iteration limit...', '\n']);
        %disp(['iSubject: ', num2str(iSubject)]);
        %disp(['iPermutation: ', num2str(iPermutation)]);
        lastwarn('');
        p_logModel_combined_output = p_logModel_combined_output(randperm(length(p_logModel_combined_output)));
        p_logModel_combined_coefficients = glmfit( [p_logModel_combined_input_rts, p_logModel_combined_input_var], ...
            [p_logModel_combined_output, p_logModel_combined_aux], 'binomial', 'link', 'logit');
        [~,lastWarningID] = lastwarn;
    end
        
        % Append the newly-generated null beta
        nullBetas_RTS_var(iSubject, iPermutation) = p_logModel_combined_coefficients(2);
        nullBetas_rts_VAR(iSubject, iPermutation) = p_logModel_combined_coefficients(3);
        
    end
    
    
    % Logistic regressions using only the immediately preceding trial.
    i_inputVector_rts = all_rts_before_infreq(:,end);
    i_inputVector_var = all_var_before_infreq(:,end);
    
    i_resorted = sortrows([i_inputVector_rts, accuracyVector, testedVector]);
    i_resorted = rmmissing(i_resorted);
    i_logModel_rts_input  = i_resorted(:,1);
    i_logModel_rts_output = i_resorted(:,2);
    i_logModel_rts_aux    = i_resorted(:,3);
    i_logModel_rts_coefficients = glmfit( i_logModel_rts_input, [i_logModel_rts_output, i_logModel_rts_aux], 'binomial', 'link', 'logit');
    
    i_resorted = sortrows([i_inputVector_var, accuracyVector, testedVector]);
    i_resorted = rmmissing(i_resorted);
    i_logModel_var_input  = i_resorted(:,1);
    i_logModel_var_output = i_resorted(:,2);
    i_logModel_var_aux    = i_resorted(:,3);
    i_logModel_var_coefficients = glmfit( i_logModel_var_input, [i_logModel_var_output, i_logModel_var_aux], 'binomial', 'link', 'logit');
    
    i_resorted = sortrows([i_inputVector_rts, i_inputVector_var, accuracyVector, testedVector]);
    i_resorted = rmmissing(i_resorted);
    i_logModel_combination_input_rts = i_resorted(:,1);
    i_logModel_combination_input_var = i_resorted(:,2);
    i_logModel_combination_output    = i_resorted(:,3);
    i_logModel_combination_aux       = i_resorted(:,4);
    i_logModel_combination_coefficients = glmfit( [i_logModel_combination_input_rts, i_logModel_combination_input_var], ...
        [i_logModel_combination_output, i_logModel_combination_aux], 'binomial', 'link', 'logit');
    i_logModel_combination_interact_coefficients = glmfit( [i_logModel_combination_input_rts, i_logModel_combination_input_var, i_logModel_combination_input_rts .* i_logModel_combination_input_var], ...
        [i_logModel_combination_output, i_logModel_combination_aux], 'binomial', 'link', 'logit');
    
    i1_logistic_coefficients_rts(iSubject,:)     = i_logModel_rts_coefficients;
    i1_logistic_coefficients_var(iSubject,:)     = i_logModel_var_coefficients;
    i1_logistic_coefficients_rts_var(iSubject,:) = i_logModel_combination_coefficients;
    i1_logistic_coefficients_rts_var_interact(iSubject,:) = i_logModel_combination_interact_coefficients;
    
    else % (cf. note above on maintaining a "good enough mix" of values in the accuracy vector)
        if length(unique(accuracyVector)) < 2
            fprintf(['WARNING: ACCURACY VECTOR FOR OLD+INFREQUENT TRIALS CONTAINS ONLY ONE VALUE TYPE', '\n']);
        else
            fprintf(['Subject (#', num2str(iSubject), ') excluded from logistic regression analyses (infrequent trials) due to unstable separation', '\n']);
        end
        % If the accuracy vector is invalid for logistic regression, then we can only leave them as NaN values.
        logistic_coefficients_rts(iSubject,:) = [NaN,NaN];
        rts_before_infreq_estimated{iSubject} = NaN;
        logistic_coefficients_var(iSubject,:) = [NaN,NaN];
        var_before_infreq_estimated{iSubject} = NaN;
        logistic_coefficients_std(iSubject,:) = [NaN,NaN];
        std_before_infreq_estimated{iSubject} = NaN;
        logistic_coefficients_rts_var(iSubject,:) = [NaN,NaN,NaN];
        logistic_coefficients_rts_var_interact(iSubject,:) = [NaN,NaN,NaN,NaN];
        nullBetas_RTS_var(iSubject,:) = nan(1,numPermutations);
        nullBetas_rts_VAR(iSubject,:) = nan(1,numPermutations);
        i1_logistic_coefficients_rts(iSubject,:)     = [NaN,NaN];
        i1_logistic_coefficients_var(iSubject,:)     = [NaN,NaN];
        i1_logistic_coefficients_rts_var(iSubject,:) = [NaN,NaN,NaN];
        i1_logistic_coefficients_rts_var_interact(iSubject,:) = [NaN,NaN,NaN,NaN];
    end
    
    
    % Logistic regressions using frequent trials.
    trial_freq_both = memData.attOrder(mem_categs==1);
    numTrial_freq_both = size(trial_freq_both,2);
    all_rts_before_freq = nan(numTrial_freq_both , numShifts);
    all_var_before_freq = nan(numTrial_freq_both , numShifts);
    
    counter=1;
    for iShift = shifts
        shifted_trials_both = trial_freq_both + iShift;
        whereValid_shiftedTrials_both = shifted_trials_both>=1;
        validIndicesOf_shiftedTrials_both = find(whereValid_shiftedTrials_both);
        valid_shifted_trials_both = shifted_trials_both(whereValid_shiftedTrials_both);
        
        all_rts_before_freq(validIndicesOf_shiftedTrials_both , counter) = validShift_RTs(valid_shifted_trials_both);
        all_var_before_freq(validIndicesOf_shiftedTrials_both , counter) = vtc(valid_shifted_trials_both);
        
        counter=counter+1;
    end
    rts_before_freq(iSubject,:) = nanmean(all_rts_before_freq);   
    comp_rts_before_freq{iSubject} = nanmean(all_rts_before_freq,2);
    var_before_freq(iSubject,:) = nanmean(all_var_before_freq);
    comp_var_before_freq{iSubject} = nanmean(all_var_before_freq,2);
    
    f_accuracyVector = mem_ratings(mem_categs==1); % get only the trials of interest (old+frequent)
    f_accuracyVector(f_accuracyVector~=4)=0; f_accuracyVector(f_accuracyVector==4)=1; % only confident responses (4) are considered correct
    f_accuracyVector=f_accuracyVector'; % transpose from row to column
    f_inputVector_rts = comp_rts_before_freq{iSubject};
    f_inputVector_var = comp_var_before_freq{iSubject};
    f_testedVector = ones(length(f_accuracyVector),1);
    
    if ( (aFlags('rememberThreshold') == true) && ...
         ((memFreq_remembered(iSubject) <= 1.0-rememberRateThreshold) && (memFreq_remembered(iSubject) >= rememberRateThreshold)) ) ...
    ...
    || ( (aFlags('rememberThreshold') == false) && ...
         (length(unique(f_accuracyVector)) == 2) )
    
    f_resorted = sortrows([f_inputVector_rts, f_accuracyVector, f_testedVector]);
    f_resorted = rmmissing(f_resorted);
    f_logModel_rts_input  = f_resorted(:,1);
    f_logModel_rts_output = f_resorted(:,2);
    f_logModel_rts_aux    = f_resorted(:,3);
    f_logModel_rts_coefficients = glmfit( f_logModel_rts_input, [f_logModel_rts_output, f_logModel_rts_aux], 'binomial', 'link', 'logit');
    
    f_resorted = sortrows([f_inputVector_var, f_accuracyVector, f_testedVector]);
    f_resorted = rmmissing(f_resorted);
    f_logModel_var_input  = f_resorted(:,1);
    f_logModel_var_output = f_resorted(:,2);
    f_logModel_var_aux    = f_resorted(:,3);
    f_logModel_var_coefficients = glmfit( f_logModel_var_input, [f_logModel_var_output, f_logModel_var_aux], 'binomial', 'link', 'logit');
    
    f_resorted = sortrows([f_inputVector_rts, f_inputVector_var, f_accuracyVector, f_testedVector]);
    f_resorted = rmmissing(f_resorted);
    f_logModel_combination_input_rts = f_resorted(:,1);
    f_logModel_combination_input_var = f_resorted(:,2);
    f_logModel_combination_output    = f_resorted(:,3);
    f_logModel_combination_aux       = f_resorted(:,4);
    f_logModel_combination_coefficients = glmfit( [f_logModel_combination_input_rts, f_logModel_combination_input_var], ...
        [f_logModel_combination_output, f_logModel_combination_aux], 'binomial', 'link', 'logit');
    f_logModel_combination_interact_coefficients = glmfit( [f_logModel_combination_input_rts, f_logModel_combination_input_var, f_logModel_combination_input_rts .* f_logModel_combination_input_var], ...
        [f_logModel_combination_output, f_logModel_combination_aux], 'binomial', 'link', 'logit');
    
    f_logistic_coefficients_rts(iSubject,:)     = f_logModel_rts_coefficients;
    f_logistic_coefficients_var(iSubject,:)     = f_logModel_var_coefficients;
    f_logistic_coefficients_rts_var(iSubject,:) = f_logModel_combination_coefficients;
    f_logistic_coefficients_rts_var_interact(iSubject,:) = f_logModel_combination_interact_coefficients;
    
    else
        if length(unique(f_accuracyVector)) < 2
            fprintf(['WARNING: ACCURACY VECTOR FOR OLD+FREQUENT TRIALS CONTAINS ONLY ONE VALUE TYPE', '\n']);
        else
            fprintf(['Subject (#', num2str(iSubject), ') excluded from logistic regression analyses (frequent trials) due to unstable separation', '\n']);
        end
        f_logistic_coefficients_rts(iSubject,:)     = [NaN,NaN];
        f_logistic_coefficients_var(iSubject,:)     = [NaN,NaN];
        f_logistic_coefficients_rts_var(iSubject,:) = [NaN,NaN,NaN];
        f_logistic_coefficients_rts_var_interact(iSubject,:) = [NaN,NaN,NaN,NaN];
    end
    
    
    %% #ASC;#QBA - Experimental/Other Analyses
    
    % #ASC
    % Figure 2a from Rosenberg et al 2013 APP paper
    if pFlags('plotSubjects_varianceTimeCourse') == true
        figure('Name',[saveSetupPrefix,'varianceTimeCourse_ID',num2str(subjectID)],'Position', [10 10 900 600]);
        hold on;
        % plot the VTC using a light-gray line
        plot(att_trialIDs, vtc                                    , '-', 'Color', [0.75,0.75,0.75]);
        % plot a smoothed VTC using a red, thick line
        plot(att_trialIDs, vtc_smooth                             , 'r-', 'LineWidth', 3);
        % draw a black horizontal line through the time series, showing the median split of the vtc
        plot(att_trialIDs, vtc_zoneSplit*ones(1,numAttTrialsTotal), '-', 'Color', [0,0,0], 'LineWidth', 3);
        % other plotting stuff
        set(gca,'fontsize',20,'ylim',[-1,6]);
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
    
    % #AMC
    [ rts_var_rho(iSubject) , rts_var_rho_pValue(iSubject) ] = corr( rmmissing(validShift_RTs)', rmmissing(vtc)' );
    
    % #OTHER - trying to disentangle rt, rtVar correlation to find variance uniquely predicted by rtVar
    % calculate correlation between RT sliding window means and VTC this time
    % @todo: this code is absolute shit
    %slidingAverage_rts = movmean(validShift_RTs,[2,0],'omitnan','Endpoints','discard');
    slidingAverage_rts = movmean( validShift_RTs, [4,0], 'omitnan' );
    if false
        figure;
        hold on;
        plot(slidingAverage_rts);
        plot( zscore(interpolated_vtc(3:end))*0.1, 'r-' );
        hold off;
    end
    % using interpolated vtc
%    [someVar1, someVar1_tf] = rmmissing(slidingAverage_rts);
%    someVar2 = interpolated_vtc;
%    someVar2(someVar1_tf) = [];
    % not using interpolated vtc
    [someVar2, someVar2_tf] = rmmissing(vtc);
    someVar1 = slidingAverage_rts;
    someVar1(someVar2_tf) = [];
    [blah,blah_pVal] = corr( someVar1', someVar2' );
    mrts_var_corr_collector(iSubject,:) = [blah,blah_pVal];
    
    % #QBA
    % Compute quartile-based statistics for replicating Figure 3 from Rosenberg et al, 2013
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


%% #PLOTS - Plotting

% #AMC - Histogram of subject-level correlations between RT and RT variability
figure('Name',[saveSetupPrefix,'AMC-correlations']);
histogram(rts_var_rho,10);
title("Subject-level Correlations")
set(gca,'xtick',-0.2:0.1:0.8,'xlim',[-0.2,0.8],'ylim',[0,10])
ylabel('Number of Subjects');
xlabel('Pearson''s r');


% #ASC - Figures 2b and 2c from Rosenberg et al, 2013 (AP&P)
% (Note: MATLAB does not support error bars for histograms)
if pFlags('plotOverall_compareAttentionStates_errors') == true
    
    figure('Name',[saveSetupPrefix,'ASC-errors'],'Position', [10 10 900 600]);
    
    subplot(1,2,1);
    set(gca,'fontsize',20);
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
    set(gca,'fontsize',20);
    hold on;
    errorbar(1,nanmean(att_miss_ITZ),nanstd(att_miss_ITZ)/sqrt(numSubjects),'go','linewidth',4)
    errorbar(2,nanmean(att_miss_OTZ),nanstd(att_miss_OTZ)/sqrt(numSubjects),'ro','linewidth',4)
    set(gca,'xtick',[],'xlim',[0.0,3.0],'ylim',[0.0,0.05],'ytick',[0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10])
    title("Misses: ITZ vs OTZ")
    ylabel('Error Rate (%)');
    xlabel('');

    hold off;
    
    % @todo: use a paired ttest to compare ITZ against OTZ?
    
end


% #QBA - Figure 3 from Rosenberg et al, 2013 (code was largely taken from MdB)
% (Note: QBA calculations do not use detrended RTs, even if detrend is set to true)
if pFlags('plotOverall_QBA_vigilanceDecrements') == true
    
    figure('Name',[saveSetupPrefix,'QBA-vigilance']);
    
    % Figure 3a from Rosenberg, et al, 2013 (AP&P)
    subplot(2,2,1);
    hold on;
    shadedErrorBar(1:4,mean(quartile_commission_errors),std(quartile_commission_errors)/sqrt(numSubjects),'lineprops','k');
    plot(1:4,mean(quartile_commission_errors),'k');
    set(gca,'ylim',[25,40],'ytick',[25 30 35 40]);
    ylabel('Commission errors')
    xlabel('Quartile')
    
    % Figure 3b from Rosenberg, et al, 2013 (AP&P)
    subplot(2,2,2);
    hold on;
    shadedErrorBar(1:4,mean(quartile_rtvar),std(quartile_rtvar)/sqrt(numSubjects),'lineprops','k');
    plot(1:4,mean(quartile_rtvar),'k');
    ylabel('RT variability')
    xlabel('Quartile')
    
    % Figure 3c from Rosenberg, et al, 2013 (AP&P)
    subplot(2,2,3);
    hold on;
    shadedErrorBar(1:4,mean(quartile_ommission_errors),std(quartile_ommission_errors)/sqrt(numSubjects),'lineprops','k');
    plot(1:4,mean(quartile_ommission_errors),'k');
    set(gca,'ylim',[0,5],'ytick',[1 2 3 4 5]);
    ylabel('Omission errors')
    xlabel('Quartile')
    
    % Figure 3d from Rosenberg, et al, 2013 (AP&P)
    subplot(2,2,4);
    hold on;
    shadedErrorBar(1:4,mean(quartile_rt),std(quartile_rt)/sqrt(numSubjects),'lineprops','k');
    plot(1:4,mean(quartile_rt),'k');
    set(gca,'xlim',[.5,4.5],'xtick',[1 2 3 4]);
    ylabel('RT')
    xlabel('Quartile')
    
end


% #AMPMP
% (Note: rmmissing() needed for the actual values (the inputVector) in order to match the size of the estimated values vector)

% Plot logistic regression for all subjects (RTs)
if pFlags('plotOverall_logistic_RTs') == true
    
    figure('Name',[saveSetupPrefix,'logisticRT'],'Position', [10 10 900 600]);
    hold on;
    
    for iSubject = 1:numSubjects
        sorted_actualValues = rmmissing(sort(comp_rts_before_infreq{iSubject}));
        plot(sorted_actualValues, rts_before_infreq_estimated{iSubject}, '-','LineWidth',2);
    end
    
    set(gca,'fontsize',20);
    title("Logistic Regression, all subject RTs");
    xlabel("RT");
    ylabel("Memory Performance");
    set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);
    
end

% Plot logistic regression for all subjects (Var)
if pFlags('plotOverall_logistic_RTDevs') == true
    
    figure('Name',[saveSetupPrefix,'logisticVar'],'Position', [10 10 900 600]);
    hold on;
    
    for iSubject = 1:numSubjects
        sorted_actualValues = rmmissing(sort(comp_var_before_infreq{iSubject}));
        plot(sorted_actualValues, var_before_infreq_estimated{iSubject}, '-','LineWidth',2);
    end
    
    set(gca,'fontsize',20);
    title("Logistic Regression, all subject RT Vars");
    xlabel("RT Vars");
    ylabel("Memory Performance");
    set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);
    
end

% Plot logistic regression for all subjects (StD)
if pFlags('plotOverall_logistic_other') == true
    
    figure('Name',[saveSetupPrefix,'logisticSTD'],'Position', [10 10 900 600]);
    hold on;
    
    for iSubject = 1:numSubjects
        sorted_actualValues = rmmissing(sort(std_before_infreq_all{iSubject}));
        plot(sorted_actualValues, std_before_infreq_estimated{iSubject}, '-','LineWidth',2);
    end
    
    set(gca,'fontsize',20);
    title("Logistic Regression, all subject RT StDs");
    xlabel("RT StDs");
    ylabel("Memory Performance");
    set(gca,'ytick',[0.0,1.0], 'ylim',[0.0,1.0]);
    
end


% #AMPTP
% @todo: clean up this section
%plot across all subjects (note: the difference between detrended and non-detrended plots is just the y-axis scaling)
% plotting variables
if aFlags('detrend') == true
    yLabelPrefix_general = 'Detrended ';
    yLabelSuffix_RTs    = '';
    yLabelSuffix_RTDevs = '';
    yLimit_RTs =    [-0.200 , 0.200];
    yTicks_RTs =    [-0.200 , -0.100 , 0.000 , 0.100 , 0.200];
    yLimit_RTDevs = [ 0.000 , 1.200];
    yTicks_RTDevs = [ 0.000 ,  0.400 , 0.800 , 1.200];
elseif aFlags('detrend') == false
    yLabelPrefix_general = '';
    yLabelSuffix_RTs    = ' [s]';
    yLabelSuffix_RTDevs = '';
    yLimit_RTs =    [-0.200 , 0.700];
    yTicks_RTs =    [ 0.000 , 0.100 ,  0.200 , 0.300 , 0.400 , 0.500 , 0.600 , 0.700];
    yLimit_RTDevs = [ 0.000 , 1.200];
    yTicks_RTDevs = [ 0.000 ,  0.400 , 0.800 , 1.200];
end
yLabelSuffix_general = ', avg across subjects';

title_average = "Average";
title_spaghetti = "All Subjects";
yLabel_RTs    = [ yLabelPrefix_general , 'RTs'     , yLabelSuffix_general , yLabelSuffix_RTs    ];
yLabel_RTDevs = [ yLabelPrefix_general , 'RT Devs' , yLabelSuffix_general , yLabelSuffix_RTDevs ];
xLabel = "Trials before infrequent trial";
legendLabels = ["Accurate","Inaccurate"];

% Plot average preceding RTs
if pFlags('plotOverall_precedingTrials_RTs') == true
    
    figure('Name',[saveSetupPrefix,'precedingRT'],'Position', [10 10 900 600]);
    
    subplot(1,2,1);
    set(gca,'fontsize',20);
    hold on;
    errorbar(shifts, nanmean(rts_before_acc)  , nanstd(rts_before_acc)  / sqrt(numSubjects), 'b', 'linewidth', 4)
    errorbar(shifts, nanmean(rts_before_inacc), nanstd(rts_before_inacc)/ sqrt(numSubjects), 'm', 'linewidth', 4)
    set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',yLimit_RTs,'ytick',yTicks_RTs)
    title( title_average )
    ylabel( yLabel_RTs );
    xlabel( xLabel );
    legend( legendLabels );
    legend boxoff
    
    subplot(1,2,2);
    set(gca,'fontsize',20);
    hold on;
    plot(shifts,rts_before_acc,'b','linewidth',1)
    plot(shifts,rts_before_inacc,'m','linewidth',1)
    plot(shifts,nanmean(rts_before_acc),'b','linewidth',4)
    plot(shifts,nanmean(rts_before_inacc),'m','linewidth',4)
    set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',yLimit_RTs,'ytick',yTicks_RTs)
    title( title_spaghetti )
    xlabel( xLabel );
    
    hold off;
    
end

% Plot average preceding RT Variances
if pFlags('plotOverall_precedingTrials_RTDevs') == true
    
    figure('Name',[saveSetupPrefix,'precedingVar'],'Position', [10 10 900 600]);
    
    subplot(1,2,1);
    set(gca,'fontsize',20);
    hold on;
    errorbar(shifts, nanmean(var_before_acc)  , nanstd(var_before_acc)  /sqrt(numSubjects), 'b', 'linewidth', 4)
    errorbar(shifts, nanmean(var_before_inacc), nanstd(var_before_inacc)/sqrt(numSubjects), 'm', 'linewidth', 4)
    set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',yLimit_RTDevs,'ytick',yTicks_RTDevs)
    title( title_average )
    ylabel( yLabel_RTDevs );
    xlabel( xLabel );
    legend( legendLabels );
    legend boxoff
    
    subplot(1,2,2);
    set(gca,'fontsize',20);
    hold on;
    plot(shifts,var_before_acc,'b','linewidth',1)
    plot(shifts,var_before_inacc,'m','linewidth',1)
    plot(shifts,nanmean(var_before_acc),'b','linewidth',4)
    plot(shifts,nanmean(var_before_inacc),'m','linewidth',4)
    set(gca,'xtick',shifts,'xlim',[shifts(1)-.5,shifts(end)+.5],'ylim',yLimit_RTDevs,'ytick',yTicks_RTDevs)
    title( title_spaghetti )
    xlabel( xLabel );
    
    hold off;
    
end

% Plot other metrics for preceding trials (as of now, it's just standard deviations)
if pFlags('plotOverall_precedingTrials_other') == true
    
    figure('Name',[saveSetupPrefix,'precedingSTD'],'Position', [10 10 450 600]);
    
    %subplot(1,2,1);
    set(gca,'fontsize',20);
    hold on;
    errorbar(1,nanmean(std_before_acc),nanstd(std_before_acc)/sqrt(numSubjects),'b','linewidth',4)
    errorbar(1,nanmean(std_before_inacc),nanstd(std_before_inacc)/sqrt(numSubjects),'m','linewidth',4)
    title(title_average)
    set(gca,'xtick',1,'xlim',[0.5,1.5],'ylim',[0.000,0.200],'ytick',[0.000,0.100,0.200])
    ylabel('RT StDs, avg across subjects [s]');
    xlabel(['StD of Preceding Trials (n=' num2str(stdWindow) ')']);
    legend('Accurate','Inaccurate');
    legend boxoff
    
    hold off;
    
end


%% #STATS - Summary Statistics
% @todo: replicate all of the numbers in deBettencourt et al, 2018

% Create some "data containers" for the output of our analyses.
% (This gives a more organised alternative than having to sift for results through the massive workspace.)
EXPT1_ATTENTION      = struct();
EXPT1_ATTENTION_BOOT = struct();
EXPT1_MEMORY         = struct();
EXPT1_MEMORY_BOOT    = struct();

EXPT2_ATTENTION      = struct();
EXPT2_ATTENTION_BOOT = struct();
EXPT2_MEMORY         = struct();
EXPT2_MEMORY_BOOT    = struct();

bs = subjectID_list_bootstrap; % shorten variable names for convenience
nb = numSamples_bootstrap;
CI95 = [2.5, 97.5];



%%% #STATSA - Attention Analyses (Non-Bootstrapped)

fprintf(outputLog, ['\n\n\n', '--- NON-BOOTSTRAPPED RESULTS [ATTENTION] ---', '\n']);

% Calculate A' for attention.
chance = 0.5;
A_Prime_attention = A_Prime( att_hit, att_falseAlarm, chance );

% 2-sample t-test for whether StDs of preceding RTs are significantly different
% between correct and incorrect trials.
[ ~ , std_ttest2_accVsInacc_p , std_ttest2_accVsInacc_ci , ~ ] = ...
    ttest2( std_before_acc , std_before_inacc );

% 1-sample t-test for whether the subject-level correlation coefficients (between RTs and RT variability) are
% significantly different from 0.
[ ~ , rts_var_rho_ttest_p  , ~ , ~ ] = ...
    ttest( rts_var_rho );


% Print outputs.
EXPT1_ATTENTION.APRIME = nanmean(A_Prime_attention);
fprintf(outputLog, ['CATEGORY JUDGMENT SENSITIVITY: ', ...
    num2str(EXPT1_ATTENTION.APRIME), ...
    '\n']);

EXPT1_ATTENTION.ERROR_INFREQUENT = mean(att_falseAlarm);
fprintf(outputLog, ['ERROR RATE - INFREQUENT: ', ...
    num2str(EXPT1_ATTENTION.ERROR_INFREQUENT), ...
    '\n']);

EXPT1_ATTENTION.ERROR_FREQUENT = mean(att_miss);
fprintf(outputLog, ['ERROR RATE - FREQUENT: ', ...
    num2str(EXPT1_ATTENTION.ERROR_FREQUENT), ...
    '\n']);

EXPT1_ATTENTION.RT_BEFORE_CORRECT = mean(nanmean(rts_before_acc));
fprintf(outputLog, ['OVERALL AVERAGE RT BEFORE CORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION.RT_BEFORE_CORRECT), ...
    ' seconds', '\n']);

EXPT1_ATTENTION.RT_BEFORE_INCORRECT = mean(nanmean(rts_before_inacc));
fprintf(outputLog, ['OVERALL AVERAGE RT BEFORE INCORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION.RT_BEFORE_INCORRECT), ...
    ' seconds', '\n']);

EXPT1_ATTENTION.VAR_BEFORE_CORRECT = mean(nanmean(var_before_acc));
fprintf(outputLog, ['OVERALL AVERAGE RT VARIANCE BEFORE CORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION.VAR_BEFORE_CORRECT), ...
    '\n']);

EXPT1_ATTENTION.VAR_BEFORE_INCORRECT = mean(nanmean(var_before_inacc));
fprintf(outputLog, ['OVERALL AVERAGE RT VARIANCE BEFORE INCORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION.VAR_BEFORE_INCORRECT), ...
    '\n']);

EXPT1_ATTENTION.STD_BEFORE_CORRECT = mean(std_before_acc);
fprintf(outputLog, ['OVERALL AVERAGE RT STANDARD DEVIATION BEFORE CORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION.STD_BEFORE_CORRECT), ...
    '\n']);

EXPT1_ATTENTION.STD_BEFORE_INCORRECT = mean(std_before_inacc);
fprintf(outputLog, ['OVERALL AVERAGE RT STANDARD DEVIATION BEFORE INCORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION.STD_BEFORE_INCORRECT), ...
    '\n']);

EXPT1_ATTENTION.RT_COR_INCOR_DIFF = nanmean(rts_before_acc) - nanmean(rts_before_inacc);
fprintf(outputLog, ['SEPARATE TRIAL ANALYSES, RT DIFFERENCE BETWEEN CORRECT AND INCORRECT, index_distance=[3,2,1]: ', ...
    num2str(EXPT1_ATTENTION.RT_COR_INCOR_DIFF), ...
    ' (seconds)', '\n']);
% intuitively, should be all positive values because accurate is slower

EXPT1_ATTENTION.VAR_COR_INCOR_DIFF = nanmean(var_before_acc) - nanmean(var_before_inacc);
fprintf(outputLog, ['SEPARATE TRIAL ANALYSES, RT VARIANCE DIFFERENCE BETWEEN CORRECT AND INCORRECT, index_distance=[3,2,1]: ', ...
    num2str(EXPT1_ATTENTION.VAR_COR_INCOR_DIFF), ...
    '\n']);
% intuitively, should be all negative values because accurate is lower variability

EXPT1_ATTENTION.STD_COR_INCOR_DIFF = std_ttest2_accVsInacc_p;
fprintf(outputLog, ['PAIRED T-TEST, RT STANDARD DEVIATION DIFFERENCE BETWEEN CORRECT AND INCORRECT, subject-level: ', ...
    'p-value = ' , num2str(EXPT1_ATTENTION.STD_COR_INCOR_DIFF), ...
    '\n']);

EXPT1_ATTENTION.RT_VAR_CORRELATION   = rts_var_rho';
EXPT1_ATTENTION.RT_VAR_CORRELATION_P = rts_var_rho_pValue';
fprintf(outputLog, ['PEARSON CORRELATION BETWEEN RT AND RT VARIANCE, subject-level: ', ...
    '\n' , num2str(EXPT1_ATTENTION.RT_VAR_CORRELATION) ...
    '\n' , '    p-values = ' , num2str(EXPT1_ATTENTION.RT_VAR_CORRELATION_P), ...
    '\n']);

EXPT1_ATTENTION.RT_VAR_CORRELATION_TTEST_P = rts_var_rho_ttest_p;
fprintf(outputLog, ['SIGNIFICANCE OF SUBJECT-LEVEL PEARSON CORRELATION COEFFICIENTS: ', ...
    'p-value (t-test) = ' , num2str(EXPT1_ATTENTION.RT_VAR_CORRELATION_TTEST_P), ...
    '\n']);



%%% #STATSAB - Attention Analyses (Bootstrapped)

fprintf(outputLog, ['\n\n\n', '--- BOOTSTRAPPED RESULTS [ATTENTION] ---', '\n']);

% Bootstrap attention accuracy
att_hit_booted   = att_hit( bs );
att_hit_allBoots = nanmean( att_hit_booted, 2 );
att_hit_avgBoot  = mean( att_hit_allBoots );
att_hit_bootCI   = prctile( att_hit_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.CORRECT_FREQUENT_MEAN = att_hit_avgBoot;
    EXPT1_ATTENTION_BOOT.CORRECT_FREQUENT_CI   = att_hit_bootCI;

att_falseAlarm_booted   = att_falseAlarm( bs );
att_falseAlarm_allBoots = nanmean( att_falseAlarm_booted, 2 );
att_falseAlarm_avgBoot  = mean( att_falseAlarm_allBoots );
att_falseAlarm_bootCI   = prctile( att_falseAlarm_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.ERROR_INFREQUENT_MEAN = att_falseAlarm_avgBoot;
    EXPT1_ATTENTION_BOOT.ERROR_INFREQUENT_CI   = att_falseAlarm_bootCI;

att_miss_booted   = att_miss( bs );
att_miss_allBoots = nanmean( att_miss_booted, 2 );
att_miss_avgBoot  = mean( att_miss_allBoots );
att_miss_bootCI   = prctile( att_miss_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.ERROR_FREQUENT_MEAN = att_miss_avgBoot;
    EXPT1_ATTENTION_BOOT.ERROR_FREQUENT_CI   = att_miss_bootCI;

att_correctReject_booted   = att_correctReject( bs );
att_correctReject_allBoots = nanmean( att_correctReject_booted, 2 );
att_correctReject_avgBoot  = mean( att_correctReject_allBoots );
att_correctReject_bootCI   = prctile( att_correctReject_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.CORRECT_INFREQUENT_MEAN = att_correctReject_avgBoot;
    EXPT1_ATTENTION_BOOT.CORRECT_INFREQUENT_CI   = att_correctReject_bootCI;

A_Prime_attention_booted   = A_Prime_attention( bs );
A_Prime_attention_allBoots = nanmean( A_Prime_attention_booted, 2 );
A_Prime_attention_avgBoot  = mean( A_Prime_attention_allBoots );
A_Prime_attention_bootCI   = prctile( A_Prime_attention_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.APRIME_MEAN = A_Prime_attention_avgBoot;
    EXPT1_ATTENTION_BOOT.APRIME_CI   = A_Prime_attention_bootCI;

% Bootstrap attention metrics
rts_before_acc_booted   = nan(nb,numSubjects,numShifts); % e.g. 100000x32x3 matrix
rts_before_inacc_booted = nan(nb,numSubjects,numShifts);
var_before_acc_booted   = nan(nb,numSubjects,numShifts);
var_before_inacc_booted = nan(nb,numSubjects,numShifts);
for iShift=1:numShifts
    
    rts_before_acc_iShift               = rts_before_acc(:,iShift);
    rts_before_acc_booted(:,:,iShift)   = rts_before_acc_iShift( bs );
    
    rts_before_inacc_iShift             = rts_before_inacc(:,iShift);
    rts_before_inacc_booted(:,:,iShift) = rts_before_inacc_iShift( bs );
    
    var_before_acc_iShift               = var_before_acc(:,iShift);
    var_before_acc_booted(:,:,iShift)   = var_before_acc_iShift( bs );
    
    var_before_inacc_iShift             = var_before_inacc(:,iShift);
    var_before_inacc_booted(:,:,iShift) = var_before_inacc_iShift( bs );
    
end

rts_before_acc_allBoots = reshape( nanmean( rts_before_acc_booted, 2 ) , nb,numShifts );
rts_before_acc_avgBoot  = reshape( mean( rts_before_acc_allBoots, 1 ) , 1,numShifts );
rts_before_acc_bootCI   = prctile( rts_before_acc_allBoots, CI95, 1 );
    EXPT1_ATTENTION_BOOT.RT_BEFORE_CORRECT_MEAN = rts_before_acc_avgBoot;
    EXPT1_ATTENTION_BOOT.RT_BEFORE_CORRECT_CI   = rts_before_acc_bootCI;
rts_before_acc_overall_allBoots = reshape( mean( rts_before_acc_allBoots, 2 ) , nb,1 );
rts_before_acc_overall_avgBoot  = mean(rts_before_acc_overall_allBoots);
rts_before_acc_overall_bootCI   = prctile( rts_before_acc_overall_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.RT_BEFORE_CORRECT_OVERALL_MEAN = rts_before_acc_overall_avgBoot;
    EXPT1_ATTENTION_BOOT.RT_BEFORE_CORRECT_OVERALL_CI   = rts_before_acc_overall_bootCI;

rts_before_inacc_allBoots = reshape( nanmean( rts_before_inacc_booted, 2 ) , nb,numShifts );
rts_before_inacc_avgBoot  = reshape( mean( rts_before_inacc_allBoots, 1 ) , 1,numShifts );
rts_before_inacc_bootCI   = prctile( rts_before_inacc_allBoots, CI95, 1 );
    EXPT1_ATTENTION_BOOT.RT_BEFORE_INCORRECT_MEAN = rts_before_inacc_avgBoot;
    EXPT1_ATTENTION_BOOT.RT_BEFORE_INCORRECT_CI   = rts_before_inacc_bootCI;
rts_before_inacc_overall_allBoots = reshape( mean( rts_before_inacc_allBoots, 2 ) , nb,1 );
rts_before_inacc_overall_avgBoot  = mean(rts_before_inacc_overall_allBoots);
rts_before_inacc_overall_bootCI   = prctile( rts_before_inacc_overall_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.RT_BEFORE_INCORRECT_OVERALL_MEAN = rts_before_inacc_overall_avgBoot;
    EXPT1_ATTENTION_BOOT.RT_BEFORE_INCORRECT_OVERALL_CI   = rts_before_inacc_overall_bootCI;

var_before_acc_allBoots = reshape( nanmean( var_before_acc_booted, 2 ) , nb,numShifts );
var_before_acc_avgBoot  = reshape( mean( var_before_acc_allBoots, 1 ) , 1,numShifts );
var_before_acc_bootCI   = prctile( var_before_acc_allBoots, CI95, 1 );
    EXPT1_ATTENTION_BOOT.VAR_BEFORE_CORRECT_MEAN = var_before_acc_avgBoot;
    EXPT1_ATTENTION_BOOT.VAR_BEFORE_CORRECT_CI   = var_before_acc_bootCI;
var_before_acc_overall_allBoots = reshape( mean( var_before_acc_allBoots, 2 ) , nb,1 );
var_before_acc_overall_avgBoot = mean(var_before_acc_overall_allBoots);
var_before_acc_overall_bootCI  = prctile( var_before_acc_overall_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.VAR_BEFORE_CORRECT_OVERALL_MEAN = var_before_acc_overall_avgBoot;
    EXPT1_ATTENTION_BOOT.VAR_BEFORE_CORRECT_OVERALL_CI   = var_before_acc_overall_bootCI;

var_before_inacc_allBoots = reshape( nanmean( var_before_inacc_booted, 2 ) , nb,numShifts );
var_before_inacc_avgBoot  = reshape( mean( var_before_inacc_allBoots, 1 ) , 1,numShifts );
var_before_inacc_bootCI   = prctile( var_before_inacc_allBoots, CI95, 1 );
    EXPT1_ATTENTION_BOOT.VAR_BEFORE_INCORRECT_MEAN = var_before_inacc_avgBoot;
    EXPT1_ATTENTION_BOOT.VAR_BEFORE_INCORRECT_CI   = var_before_inacc_bootCI;
var_before_inacc_overall_allBoots = reshape( mean( var_before_inacc_allBoots, 2 ) , nb,1 );
var_before_inacc_overall_avgBoot = mean(var_before_inacc_overall_allBoots);
var_before_inacc_overall_bootCI  = prctile( var_before_inacc_overall_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.VAR_BEFORE_INCORRECT_OVERALL_MEAN = var_before_inacc_overall_avgBoot;
    EXPT1_ATTENTION_BOOT.VAR_BEFORE_INCORRECT_OVERALL_CI   = var_before_inacc_overall_bootCI;

std_before_acc_booted   = std_before_acc( bs );
std_before_acc_allBoots = nanmean( std_before_acc_booted, 2 );
std_before_acc_avgBoot  = mean( std_before_acc_allBoots );
std_before_acc_bootCI   = prctile( std_before_acc_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.STD_BEFORE_CORRECT_MEAN = std_before_acc_avgBoot;
    EXPT1_ATTENTION_BOOT.STD_BEFORE_CORRECT_CI   = std_before_acc_bootCI;

std_before_inacc_booted   = std_before_inacc( bs );
std_before_inacc_allBoots = nanmean( std_before_inacc_booted, 2 );
std_before_inacc_avgBoot  = mean( std_before_inacc_allBoots );
std_before_inacc_bootCI   = prctile( std_before_inacc_allBoots, CI95 );
    EXPT1_ATTENTION_BOOT.STD_BEFORE_INCORRECT_MEAN = std_before_inacc_avgBoot;
    EXPT1_ATTENTION_BOOT.STD_BEFORE_INCORRECT_CI   = std_before_inacc_bootCI;

% Response sensitivity greater than chance?
A_Prime_attention_bootp = bootp( A_Prime_attention_allBoots, 0.5, nb, 'below' );
    EXPT1_ATTENTION_BOOT.APRIME_P = A_Prime_attention_bootp;

% Error rate for infrequent higher than error rate for frequent? (I.e. is the difference greater than 0.0?)
errorRate_difference_booted = att_falseAlarm_booted - att_miss_booted;
errorRate_difference_allBoots = nanmean( errorRate_difference_booted, 2 );
errorRate_betweenCategories_bootp = bootp( errorRate_difference_allBoots, 0.0, nb , 'below' );
    EXPT1_ATTENTION_BOOT.ERR_INFREQ_FREQ_P = errorRate_betweenCategories_bootp;

% Preceding response times (overall) are higher for correct responses than incorrect ones?
responseTimes_overallDifference_booted = rts_before_acc_booted - rts_before_inacc_booted;
responseTimes_overallDifference_allBoots = reshape( nanmean( responseTimes_overallDifference_booted, 2 ) , nb,numShifts );
responseTimes_overallDifference_allBoots = nanmean( responseTimes_overallDifference_allBoots, 2 );
responseTimes_overallBetweenPerformance_bootp = bootp( responseTimes_overallDifference_allBoots, 0.0, nb, 'below' );
    EXPT1_ATTENTION_BOOT.RT_COR_INCOR_SUMDIFF_P = responseTimes_overallBetweenPerformance_bootp;

% Preceding response time variances (overall) are higher for incorrect responses than correct ones?
rtVars_overallDifference_booted = var_before_inacc_booted - var_before_acc_booted;
rtVars_overallDifference_allBoots = reshape( nanmean( rtVars_overallDifference_booted, 2 ) , nb,numShifts );
rtVars_overallDifference_allBoots = nanmean( rtVars_overallDifference_allBoots, 2 );
rtVars_overallBetweenPerformance_bootp = bootp( rtVars_overallDifference_allBoots, 0.0, nb, 'below' );
    EXPT1_ATTENTION_BOOT.VAR_COR_INCOR_SUMDIFF_P = rtVars_overallBetweenPerformance_bootp;

% Perform the same analyses, but this time do it for each iShift individually.
responseTimes_difference_booted = rts_before_acc_booted - rts_before_inacc_booted;
responseTimes_difference_allBoots = reshape( nanmean( responseTimes_difference_booted, 2 ) , nb,numShifts );
responseTimes_betweenPerformance_bootp = nan(1,numShifts);
rtVars_difference_booted = var_before_inacc_booted - var_before_acc_booted;
rtVars_difference_allBoots = reshape( nanmean( rtVars_difference_booted, 2 ) , nb,numShifts );
rtVars_betweenPerformance_bootp = nan(1,numShifts);
for iShift=1:numShifts
    
    % Preceding response times (for each shift) are higher for correct responses than incorrect ones?
    responseTimes_betweenPerformance_bootp(iShift) = bootp( responseTimes_difference_allBoots(:,iShift), 0.0, nb, 'below' );
    
    % Preceding response time variances (for each shift) are higher for incorrect responses than correct ones?
    rtVars_betweenPerformance_bootp(iShift) = bootp( rtVars_difference_allBoots(:,iShift), 0.0, nb, 'below' );
    
end
    EXPT1_ATTENTION_BOOT.RT_COR_INCOR_DIFF_MEAN = mean(responseTimes_difference_allBoots);
    EXPT1_ATTENTION_BOOT.VAR_COR_INCOR_DIFF_MEAN = mean(rtVars_difference_allBoots);
    EXPT1_ATTENTION_BOOT.RT_COR_INCOR_DIFF_P = responseTimes_betweenPerformance_bootp;
    EXPT1_ATTENTION_BOOT.VAR_COR_INCOR_DIFF_P = rtVars_betweenPerformance_bootp;

% Preceding response time standard deviations are higher for incorrect responses than correct ones?
rtSTDs_difference_booted = std_before_inacc_booted - std_before_acc_booted;
rtSTDs_difference_allBoots = nanmean( rtSTDs_difference_booted );
rtSTDs_betweenPerformance_bootp = bootp( rtSTDs_difference_allBoots, 0.0, nb, 'below' );
    EXPT1_ATTENTION_BOOT.STD_COR_INCOR_DIFF_P = rtSTDs_betweenPerformance_bootp;


% Print outputs.
fprintf(outputLog, ['CATEGORY JUDGMENT SENSITIVITY: ', ...
    num2str(EXPT1_ATTENTION_BOOT.APRIME_MEAN), ...
    ' (', num2str(EXPT1_ATTENTION_BOOT.APRIME_CI(1)), ', ', num2str(EXPT1_ATTENTION_BOOT.APRIME_CI(2)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_ATTENTION_BOOT.APRIME_P), ...
    '\n']);

fprintf(outputLog, ['ERROR RATE - INFREQUENT: ', ...
    num2str(EXPT1_ATTENTION_BOOT.ERROR_INFREQUENT_MEAN), ...
    ' (', num2str(EXPT1_ATTENTION_BOOT.ERROR_INFREQUENT_CI(1)), ', ', num2str(EXPT1_ATTENTION_BOOT.ERROR_INFREQUENT_CI(2)), ')', ...
    '\n']);

fprintf(outputLog, ['ERROR RATE - FREQUENT: ', ...
    num2str(EXPT1_ATTENTION_BOOT.ERROR_FREQUENT_MEAN), ...
    ' (', num2str(EXPT1_ATTENTION_BOOT.ERROR_FREQUENT_CI(1)), ', ', num2str(EXPT1_ATTENTION_BOOT.ERROR_FREQUENT_CI(2)), ')', ...
    '\n']);

fprintf(outputLog, ['ERROR RATE - SIGNIFICANCE (NON-PARAMETRIC) OF DIFFERENCE BETWEEN INFREQUENT AND FREQUENT: ', ...
    'p-value = ' , num2str(EXPT1_ATTENTION_BOOT.ERR_INFREQ_FREQ_P), ...
    '\n']);

fprintf(outputLog, ['OVERALL AVERAGE RT BEFORE CORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION_BOOT.RT_BEFORE_CORRECT_OVERALL_MEAN), ...
    ' (', num2str(EXPT1_ATTENTION_BOOT.RT_BEFORE_CORRECT_OVERALL_CI(1)), ', ', num2str(EXPT1_ATTENTION_BOOT.RT_BEFORE_CORRECT_OVERALL_CI(2)), ')', ...
    ' seconds', '\n']);

fprintf(outputLog, ['OVERALL AVERAGE RT BEFORE INCORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION_BOOT.RT_BEFORE_INCORRECT_OVERALL_MEAN), ...
    ' (', num2str(EXPT1_ATTENTION_BOOT.RT_BEFORE_INCORRECT_OVERALL_CI(1)), ', ', num2str(EXPT1_ATTENTION_BOOT.RT_BEFORE_INCORRECT_OVERALL_CI(2)), ')', ...
    ' seconds', '\n']);

fprintf(outputLog, ['AVERAGE PRECEDING RESPONSE TIMES - SIGNIFICANCE (NON-PARAMETRIC) OF DIFFERENCE BETWEEN CORRECT AND INCORRECT: ', ...
    'p-value = ' , num2str(EXPT1_ATTENTION_BOOT.RT_COR_INCOR_SUMDIFF_P), ...
    '\n']);

fprintf(outputLog, ['OVERALL AVERAGE RT VARIANCE BEFORE CORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION_BOOT.VAR_BEFORE_CORRECT_OVERALL_MEAN), ...
    ' (', num2str(EXPT1_ATTENTION_BOOT.VAR_BEFORE_CORRECT_OVERALL_CI(1)), ', ', num2str(EXPT1_ATTENTION_BOOT.VAR_BEFORE_CORRECT_OVERALL_CI(2)), ')', ...
    '\n']);

fprintf(outputLog, ['OVERALL AVERAGE RT VARIANCE BEFORE INCORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION_BOOT.VAR_BEFORE_INCORRECT_OVERALL_MEAN), ...
    ' (', num2str(EXPT1_ATTENTION_BOOT.VAR_BEFORE_INCORRECT_OVERALL_CI(1)), ', ', num2str(EXPT1_ATTENTION_BOOT.VAR_BEFORE_INCORRECT_OVERALL_CI(2)), ')', ...
    '\n']);

fprintf(outputLog, ['AVERAGE PRECEDING RT VARIANCES - SIGNIFICANCE (NON-PARAMETRIC) OF DIFFERENCE BETWEEN CORRECT AND INCORRECT: ', ...
    'p-value = ' , num2str(EXPT1_ATTENTION_BOOT.VAR_COR_INCOR_SUMDIFF_P), ...
    '\n']);

fprintf(outputLog, ['OVERALL AVERAGE RT STANDARD DEVIATION BEFORE CORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION_BOOT.STD_BEFORE_CORRECT_MEAN), ...
    ' (', num2str(EXPT1_ATTENTION_BOOT.STD_BEFORE_CORRECT_CI(1)), ', ', num2str(EXPT1_ATTENTION_BOOT.STD_BEFORE_CORRECT_CI(2)), ')', ...
    '\n']);

fprintf(outputLog, ['OVERALL AVERAGE RT STANDARD DEVIATION BEFORE INCORRECT RESPONSE: ', ...
    num2str(EXPT1_ATTENTION_BOOT.STD_BEFORE_INCORRECT_MEAN), ...
    ' (', num2str(EXPT1_ATTENTION_BOOT.STD_BEFORE_INCORRECT_CI(1)), ', ', num2str(EXPT1_ATTENTION_BOOT.STD_BEFORE_INCORRECT_CI(2)), ')', ...
    '\n']);

fprintf(outputLog, ['AVERAGE PRECEDING RT STD - SIGNIFICANCE (NON-PARAMETRIC) OF DIFFERENCE BETWEEN CORRECT AND INCORRECT: ', ...
    'p-value = ' , num2str(EXPT1_ATTENTION_BOOT.STD_COR_INCOR_DIFF_P), ...
    '\n']);

fprintf(outputLog, ['SEPARATE TRIAL ANALYSES, RT DIFFERENCE BETWEEN CORRECT AND INCORRECT, index_distance=[3,2,1]: ', ...
    num2str(EXPT1_ATTENTION_BOOT.RT_COR_INCOR_DIFF_MEAN), ...
    ' (seconds)', '\n', '    p-values = ', num2str(EXPT1_ATTENTION_BOOT.RT_COR_INCOR_DIFF_P), ...
    '\n']);
% (Intuitively, should be all positive values because accurate is slower)

fprintf(outputLog, ['SEPARATE TRIAL ANALYSES, RT VARIANCE DIFFERENCE BETWEEN CORRECT AND INCORRECT, index_distance=[3,2,1]: ', ...
    num2str(EXPT1_ATTENTION_BOOT.VAR_COR_INCOR_DIFF_MEAN), ...
    '\n', '    p-values = ', num2str(EXPT1_ATTENTION_BOOT.VAR_COR_INCOR_DIFF_P), ...
    '\n']);
% (Intuitively, should be all negative values because accurate is lower variability)



%%% #STATSM - Memory Analyses (Non-Bootstrapped)

fprintf(outputLog, ['\n\n\n', '--- NON-BOOTSTRAPPED RESULTS [MEMORY] ---', '\n']);

% Calculate A' for memory (both categories).
chance = 0.5;
% Note: mem_hit is not the same as mem_remembered; the latter has a more stringent requirement where responses also have to be confident.
% The same applies as with mem_falseAlarm versus mem_cFalseAlarm. (The "c" stands for "confident")
A_Prime_memoryFrequent   = A_Prime( memFreq_remembered   , memFreq_cFalseAlarm   , chance );
A_Prime_memoryInfrequent = A_Prime( memInfreq_remembered , memInfreq_cFalseAlarm , chance );

% 1-sample t-test on beta coefficients to compare significance of different attention metrics.
% Note: rmmissing() is applied to the beta coefficients in case one of the subjects had an invalid accuracy vector.
% Model with RTs:
[ ~ , logistic_coefficients_rts_pValue , logistic_coefficients_rts_CI , ~ ] = ...
    ttest( rmmissing(logistic_coefficients_rts(:,2)) );
% Model with RT Vars:
[ ~ , logistic_coefficients_var_pValue , logistic_coefficients_var_CI , ~ ] = ...
    ttest( rmmissing(logistic_coefficients_var(:,2)) );
% Model with RT StDs, for completeness:
[ ~ , logistic_coefficients_std_pValue , logistic_coefficients_std_CI , ~ ] = ...
    ttest( rmmissing(logistic_coefficients_std(:,2)) );
% Model with RTs and RT Vars:
[ ~ , logistic_coefficients_RTS_var_pValue , logistic_coefficients_RTS_var_CI , ~ ] = ...
    ttest( rmmissing(logistic_coefficients_rts_var(:,2)) );
[ ~ , logistic_coefficients_rts_VAR_pValue , logistic_coefficients_rts_VAR_CI , ~ ] = ...
    ttest( rmmissing(logistic_coefficients_rts_var(:,3)) );
% Model with RTs, RT Vars, and interaction term:
[ ~ , logistic_coefficients_RTS_var_interact_pValue , logistic_coefficients_RTS_var_interact_CI , ~ ] = ...
    ttest( rmmissing(logistic_coefficients_rts_var_interact(:,2)) );
[ ~ , logistic_coefficients_rts_VAR_interact_pValue , logistic_coefficients_rts_VAR_interact_CI , ~ ] = ...
    ttest( rmmissing(logistic_coefficients_rts_var_interact(:,3)) );
[ ~ , logistic_coefficients_rts_var_INTERACT_pValue , logistic_coefficients_rts_var_INTERACT_CI , ~ ] = ...
    ttest( rmmissing(logistic_coefficients_rts_var_interact(:,4)) );

% Alternative way to get preceding attention metrics; needed for a calculation below.
% (Proof of concept that direction of means matters)
% (Compare the values of 'prevRTsInfreq_ShiftMean_Cor' and 'prevRTsInfreq_ShiftMean_Err'
%  against 'rts_before_acc' and 'rts_before_inacc')
prevRTsInfreq_ShiftMean_Cor = [];
prevRTsInfreq_ShiftMean_Err = [];
prevVarsInfreq_ShiftMean_Cor = [];
prevVarsInfreq_ShiftMean_Err = [];
for i = 1:numSubjects
    prevRTsInfreq_ShiftMean_Cor  = cat(2, prevRTsInfreq_ShiftMean_Cor, nanmean(comp_rts_before_acc{i}));
    prevRTsInfreq_ShiftMean_Err  = cat(2, prevRTsInfreq_ShiftMean_Err, nanmean(comp_rts_before_inacc{i}));
    prevVarsInfreq_ShiftMean_Cor = cat(2, prevVarsInfreq_ShiftMean_Cor, nanmean(comp_var_before_acc{i}));
    prevVarsInfreq_ShiftMean_Err = cat(2, prevVarsInfreq_ShiftMean_Err, nanmean(comp_var_before_inacc{i}));
end

% Compare each subject's overall memory performance (A_Prime) with each subject's attentional metrics (RT and RT variance).
% If they are significantly correlated, then the relationship between memory and attention is more trait-like, rather than
% good memory being a result of paying attention at that time.
% @todo: what is the bootstrapped version of this?
[ APrime_rts_rho , APrime_rts_rho_pValue ] = corr( A_Prime_memoryInfrequent' , prevRTsInfreq_ShiftMean_Err'  );
[ APrime_var_rho , APrime_var_rho_pValue ] = corr( A_Prime_memoryInfrequent' , prevVarsInfreq_ShiftMean_Err' );

%%%% correlation between RTs and RT_Vars %%%%
% corr_coefficients = [];
% for i = 1:length(rt_rtVar_corr_collector)
%     corr_coefficients = cat(2, corr_coefficients, rt_rtVar_corr_collector{i}.p1);
% end
% corr_coefficients_mean = mean(corr_coefficients)
% bootci(numSamples_bootstrap, @mean, corr_coefficients)


% Print outputs.
% @todo: store t-test outputs to struct output variables? (we also don't do t-tests on the single-term models)
EXPT1_MEMORY.APRIME_FREQUENT = nanmean(A_Prime_memoryFrequent);
fprintf(outputLog, ['MEMORY JUDGMENT SENSITIVITY (FREQUENT): ', ...
    num2str(EXPT1_MEMORY.APRIME_FREQUENT), ...
    '\n']);

EXPT1_MEMORY.APRIME_INFREQUENT = nanmean(A_Prime_memoryInfrequent);
fprintf(outputLog, ['MEMORY JUDGMENT SENSITIVITY (INFREQUENT): ', ...
    num2str(EXPT1_MEMORY.APRIME_INFREQUENT), ...
    '\n']);

EXPT1_MEMORY.APRIME_DIFFERENCE = nanmean(A_Prime_memoryInfrequent - A_Prime_memoryFrequent);
fprintf(outputLog, ['MEMORY JUDGMENT SENSITIVITY (DIFFERENCE BETWEEN INFREQUENT & FREQUENT): ', ...
    num2str(EXPT1_MEMORY.APRIME_DIFFERENCE), ...
    '\n']);

EXPT1_MEMORY.REMEMBERED_INFREQUENT = nanmean(memInfreq_remembered);
fprintf(outputLog, ['PROPORTION REMEMBERED (INFREQUENT): ', ...
    num2str(EXPT1_MEMORY.REMEMBERED_INFREQUENT), ...
    '\n']);

EXPT1_MEMORY.FORGOTTEN_INFREQUENT = nanmean(memInfreq_forgotten);
fprintf(outputLog, ['PROPORTION FORGOTTEN (INFREQUENT): ', ...
    num2str(EXPT1_MEMORY.FORGOTTEN_INFREQUENT), ...
    '\n']);

EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_RTS = nanmean(logistic_coefficients_rts(:,2));
fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT (RT): ', ...
    num2str(EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_RTS), ...
    '\n', '    p-value (t-test) = ', num2str(logistic_coefficients_rts_pValue), ...    
    '\n']);

EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_VAR = nanmean(logistic_coefficients_var(:,2));
fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT (VAR): ', ...
    num2str(EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_VAR), ...
    '\n', '    p-value (t-test) = ', num2str(logistic_coefficients_var_pValue), ...
    '\n']);

EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_STD = nanmean(logistic_coefficients_std(:,2));
fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT (STD): ', ...
    num2str(EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_STD), ...
    '\n', '    p-value (t-test) = ', num2str(logistic_coefficients_std_pValue), ...
    '\n']);

EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_2TERMS_RTS = nanmean(logistic_coefficients_rts_var(:,2));
fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT, RT+VAR MODEL (RT): ', ...
    num2str(EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_2TERMS_RTS), ...
    '\n', '    p-value (t-test) = ', num2str(logistic_coefficients_RTS_var_pValue), ...
    '\n']);

EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_2TERMS_VAR = nanmean(logistic_coefficients_rts_var(:,3));
fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT, RT+VAR MODEL (VAR): ', ...
    num2str(EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_2TERMS_VAR), ...
    '\n', '    p-value (t-test) = ', num2str(logistic_coefficients_rts_VAR_pValue), ...
    '\n']);

EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_3TERMS_RTS = nanmean(logistic_coefficients_rts_var_interact(:,2));
fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT, RT+VAR+INTERACT MODEL (RT): ', ...
    num2str(EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_3TERMS_RTS), ...
    '\n', '    p-value (t-test) = ', num2str(logistic_coefficients_RTS_var_interact_pValue), ...
    '\n']);

EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_3TERMS_VAR = nanmean(logistic_coefficients_rts_var_interact(:,3));
fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT, RT+VAR+INTERACT MODEL (VAR): ', ...
    num2str(EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_3TERMS_VAR), ...
    '\n', '    p-value (t-test) = ', num2str(logistic_coefficients_rts_VAR_interact_pValue), ...
    '\n']);

EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_3TERMS_INTERACT = nanmean(logistic_coefficients_rts_var_interact(:,4));
fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT, RT+VAR+INTERACT MODEL (INTERACT): ', ...
    num2str(EXPT1_MEMORY.LOGISTIC_COEFFICIENTS_3TERMS_INTERACT), ...
    '\n', '    p-value (t-test) = ', num2str(logistic_coefficients_rts_var_INTERACT_pValue), ...
    '\n']);

EXPT1_MEMORY.APRIME_RT_CORRELATION   = APrime_rts_rho;
EXPT1_MEMORY.APRIME_RT_CORRELATION_P = APrime_rts_rho_pValue;
fprintf(outputLog, ['PEARSON CORRELATION BETWEEN A-PRIME (INFREQUENT) AND RT: ', ...
    num2str(EXPT1_MEMORY.APRIME_RT_CORRELATION), ...
    '\n', '    p-value = ', num2str(EXPT1_MEMORY.APRIME_RT_CORRELATION_P), ...
    '\n']);

EXPT1_MEMORY.APRIME_VAR_CORRELATION   = APrime_var_rho;
EXPT1_MEMORY.APRIME_VAR_CORRELATION_P = APrime_var_rho_pValue;
fprintf(outputLog, ['PEARSON CORRELATION BETWEEN A-PRIME (INFREQUENT) AND VAR: ', ...
    num2str(EXPT1_MEMORY.APRIME_VAR_CORRELATION), ...
    '\n', '    p-value = ', num2str(EXPT1_MEMORY.APRIME_VAR_CORRELATION_P), ...
    '\n']);



%%% #STATSMB - Memory Analyses (Bootstrapped)

fprintf(outputLog, ['\n\n\n', '--- BOOTSTRAPPED RESULTS [MEMORY] ---', '\n']);

% Bootstrap memory performance.
% (In lieu of hit/miss/FA/CR, only remembered/forgotten will generally be reported.)
memFreq_remembered_booted   = memFreq_remembered( bs );
memFreq_remembered_allBoots = nanmean( memFreq_remembered_booted, 2 );
memFreq_remembered_avgBoot  = mean( memFreq_remembered_allBoots );
memFreq_remembered_bootCI   = prctile( memFreq_remembered_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.REMEMBERED_FREQUENT_MEAN = memFreq_remembered_avgBoot;
    EXPT1_MEMORY_BOOT.REMEMBERED_FREQUENT_CI   = memFreq_remembered_bootCI;

memInfreq_remembered_booted   = memInfreq_remembered( bs );
memInfreq_remembered_allBoots = nanmean( memInfreq_remembered_booted, 2 );
memInfreq_remembered_avgBoot  = mean( memInfreq_remembered_allBoots );
memInfreq_remembered_bootCI   = prctile( memInfreq_remembered_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.REMEMBERED_INFREQUENT_MEAN = memInfreq_remembered_avgBoot;
    EXPT1_MEMORY_BOOT.REMEMBERED_INFREQUENT_CI   = memInfreq_remembered_bootCI;

memFreq_forgotten_booted   = memFreq_forgotten( bs );
memFreq_forgotten_allBoots = nanmean( memFreq_forgotten_booted, 2 );
memFreq_forgotten_avgBoot  = mean( memFreq_forgotten_allBoots );
memFreq_forgotten_bootCI   = prctile( memFreq_forgotten_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.FORGOTTEN_FREQUENT_MEAN = memFreq_forgotten_avgBoot;
    EXPT1_MEMORY_BOOT.FORGOTTEN_FREQUENT_CI   = memFreq_forgotten_bootCI;

memInfreq_forgotten_booted   = memInfreq_forgotten( bs );
memInfreq_forgotten_allBoots = nanmean( memInfreq_forgotten_booted, 2 );
memInfreq_forgotten_avgBoot  = mean( memInfreq_forgotten_allBoots );
memInfreq_forgotten_bootCI   = prctile( memInfreq_forgotten_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.FORGOTTEN_INFREQUENT_MEAN = memInfreq_forgotten_avgBoot;
    EXPT1_MEMORY_BOOT.FORGOTTEN_INFREQUENT_CI   = memInfreq_forgotten_bootCI;

A_Prime_memoryFrequent_booted   = A_Prime_memoryFrequent( bs );
A_Prime_memoryFrequent_allBoots = nanmean( A_Prime_memoryFrequent_booted, 2 );
A_Prime_memoryFrequent_avgBoot  = mean( A_Prime_memoryFrequent_allBoots );
A_Prime_memoryFrequent_bootCI   = prctile( A_Prime_memoryFrequent_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.APRIME_FREQUENT_MEAN = A_Prime_memoryFrequent_avgBoot;
    EXPT1_MEMORY_BOOT.APRIME_FREQUENT_CI   = A_Prime_memoryFrequent_bootCI;

A_Prime_memoryInfrequent_booted   = A_Prime_memoryInfrequent( bs );
A_Prime_memoryInfrequent_allBoots = nanmean( A_Prime_memoryInfrequent_booted, 2 );
A_Prime_memoryInfrequent_avgBoot  = mean( A_Prime_memoryInfrequent_allBoots );
A_Prime_memoryInfrequent_bootCI   = prctile( A_Prime_memoryInfrequent_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.APRIME_INFREQUENT_MEAN = A_Prime_memoryInfrequent_avgBoot;
    EXPT1_MEMORY_BOOT.APRIME_INFREQUENT_CI   = A_Prime_memoryInfrequent_bootCI;

% Memory sensitivity significantly above chance (both categories)?
A_Prime_memoryFrequent_bootp   = bootp( A_Prime_memoryFrequent_allBoots  , 0.5, nb, 'below' );
    EXPT1_MEMORY_BOOT.APRIME_FREQUENT_P   = A_Prime_memoryFrequent_bootp;
A_Prime_memoryInfrequent_bootp = bootp( A_Prime_memoryInfrequent_allBoots, 0.5, nb, 'below' );
    EXPT1_MEMORY_BOOT.APRIME_INFREQUENT_P = A_Prime_memoryInfrequent_bootp;

% Memory sensitivity of infrequent significantly above frequent?
A_Prime_memoryDifference          = A_Prime_memoryInfrequent - A_Prime_memoryFrequent;
A_Prime_memoryDifference_booted   = A_Prime_memoryDifference( bs );
A_Prime_memoryDifference_allBoots = nanmean( A_Prime_memoryDifference_booted, 2 );
A_Prime_memoryDifference_avgBoot  = mean( A_Prime_memoryDifference_allBoots );
A_Prime_memoryDifference_bootCI   = prctile( A_Prime_memoryDifference_allBoots, CI95 );
A_Prime_memoryDifference_bootp    = bootp( A_Prime_memoryDifference_allBoots, 0.0, nb, 'below' );
    EXPT1_MEMORY_BOOT.APRIME_DIFFERENCE_MEAN = A_Prime_memoryDifference_avgBoot;
    EXPT1_MEMORY_BOOT.APRIME_DIFFERENCE_CI   = A_Prime_memoryDifference_bootCI;
    EXPT1_MEMORY_BOOT.APRIME_DIFFERENCE_P    = A_Prime_memoryDifference_bootp;

% Bootstrap logistic regression coefficients.
% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_rts_trimmed  = rmmissing(logistic_coefficients_rts( :,2 ));
logistic_coefficients_rts_trimmed=logistic_coefficients_rts(:,2);
logistic_coefficients_rts_booted   = logistic_coefficients_rts_trimmed( bs );
logistic_coefficients_rts_allBoots = nanmean( logistic_coefficients_rts_booted, 2 );
logistic_coefficients_rts_avgBoot  = mean( logistic_coefficients_rts_allBoots );
logistic_coefficients_rts_bootCI   = prctile( logistic_coefficients_rts_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_MEAN = logistic_coefficients_rts_avgBoot;
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_CI   = logistic_coefficients_rts_bootCI;

% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_var_trimmed  = rmmissing(logistic_coefficients_var( :,2 ));
logistic_coefficients_var_trimmed=logistic_coefficients_var(:,2);
logistic_coefficients_var_booted   = logistic_coefficients_var_trimmed( bs );
logistic_coefficients_var_allBoots = nanmean( logistic_coefficients_var_booted, 2 );
logistic_coefficients_var_avgBoot  = mean( logistic_coefficients_var_allBoots );
logistic_coefficients_var_bootCI   = prctile( logistic_coefficients_var_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_VAR_MEAN = logistic_coefficients_var_avgBoot;
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_VAR_CI   = logistic_coefficients_var_bootCI;

% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_std_trimmed  = rmmissing(logistic_coefficients_std( :,2 ));
logistic_coefficients_std_trimmed=logistic_coefficients_std(:,2);
logistic_coefficients_std_booted   = logistic_coefficients_std_trimmed( bs );
logistic_coefficients_std_allBoots = nanmean( logistic_coefficients_std_booted, 2 );
logistic_coefficients_std_avgBoot  = mean( logistic_coefficients_std_allBoots );
logistic_coefficients_std_bootCI   = prctile( logistic_coefficients_std_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_STD_MEAN = logistic_coefficients_std_avgBoot;
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_STD_CI   = logistic_coefficients_std_bootCI;

% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_RTS_var_trimmed  = rmmissing(logistic_coefficients_rts_var( :,2 ));
%logistic_coefficients_rts_VAR_trimmed  = rmmissing(logistic_coefficients_rts_var( :,3 ));
logistic_coefficients_RTS_var_trimmed=logistic_coefficients_rts_var(:,2);
logistic_coefficients_rts_VAR_trimmed=logistic_coefficients_rts_var(:,3);
logistic_coefficients_RTS_var_booted   = logistic_coefficients_RTS_var_trimmed( bs );
logistic_coefficients_rts_VAR_booted   = logistic_coefficients_rts_VAR_trimmed( bs );
logistic_coefficients_RTS_var_allBoots = nanmean( logistic_coefficients_RTS_var_booted, 2 );
logistic_coefficients_rts_VAR_allBoots = nanmean( logistic_coefficients_rts_VAR_booted, 2 );
logistic_coefficients_RTS_var_avgBoot  = mean( logistic_coefficients_RTS_var_allBoots );
logistic_coefficients_rts_VAR_avgBoot  = mean( logistic_coefficients_rts_VAR_allBoots );
logistic_coefficients_RTS_var_bootCI   = prctile( logistic_coefficients_RTS_var_allBoots, CI95 );
logistic_coefficients_rts_VAR_bootCI   = prctile( logistic_coefficients_rts_VAR_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_MEAN = ...
        [logistic_coefficients_RTS_var_avgBoot, logistic_coefficients_rts_VAR_avgBoot];
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_CI   = ...
        [logistic_coefficients_RTS_var_bootCI , logistic_coefficients_rts_VAR_bootCI ];

% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_RTS_var_interact_trimmed  = rmmissing(logistic_coefficients_rts_var_interact( :,2 ));
%logistic_coefficients_rts_VAR_interact_trimmed  = rmmissing(logistic_coefficients_rts_var_interact( :,3 ));
%logistic_coefficients_rts_var_INTERACT_trimmed  = rmmissing(logistic_coefficients_rts_var_interact( :,4 ));
logistic_coefficients_RTS_var_interact_trimmed=logistic_coefficients_rts_var_interact(:,2);
logistic_coefficients_rts_VAR_interact_trimmed=logistic_coefficients_rts_var_interact(:,3);
logistic_coefficients_rts_var_INTERACT_trimmed=logistic_coefficients_rts_var_interact(:,4);
logistic_coefficients_RTS_var_interact_booted   = logistic_coefficients_RTS_var_interact_trimmed( bs );
logistic_coefficients_rts_VAR_interact_booted   = logistic_coefficients_rts_VAR_interact_trimmed( bs );
logistic_coefficients_rts_var_INTERACT_booted   = logistic_coefficients_rts_var_INTERACT_trimmed( bs );
logistic_coefficients_RTS_var_interact_allBoots = nanmean( logistic_coefficients_RTS_var_interact_booted, 2 );
logistic_coefficients_rts_VAR_interact_allBoots = nanmean( logistic_coefficients_rts_VAR_interact_booted, 2 );
logistic_coefficients_rts_var_INTERACT_allBoots = nanmean( logistic_coefficients_rts_var_INTERACT_booted, 2 );
logistic_coefficients_RTS_var_interact_avgBoot  = mean( logistic_coefficients_RTS_var_interact_allBoots );
logistic_coefficients_rts_VAR_interact_avgBoot  = mean( logistic_coefficients_rts_VAR_interact_allBoots );
logistic_coefficients_rts_var_INTERACT_avgBoot  = mean( logistic_coefficients_rts_var_INTERACT_allBoots );
logistic_coefficients_RTS_var_interact_bootCI   = prctile( logistic_coefficients_RTS_var_interact_allBoots, CI95 );
logistic_coefficients_rts_VAR_interact_bootCI   = prctile( logistic_coefficients_rts_VAR_interact_allBoots, CI95 );
logistic_coefficients_rts_var_INTERACT_bootCI   = prctile( logistic_coefficients_rts_var_INTERACT_allBoots, CI95 );
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_MEAN = ...
        [logistic_coefficients_RTS_var_interact_avgBoot, logistic_coefficients_rts_VAR_interact_avgBoot, logistic_coefficients_rts_var_INTERACT_avgBoot];
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_CI   = ...
        [logistic_coefficients_RTS_var_interact_bootCI , logistic_coefficients_rts_VAR_interact_bootCI , logistic_coefficients_rts_var_INTERACT_bootCI ];

% Are the beta coefficients significantly different from 0?
logistic_coefficients_rts_bootp = bootp( logistic_coefficients_rts_allBoots, 0.0, nb, 'below' );
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_P = logistic_coefficients_rts_bootp;
logistic_coefficients_var_bootp = bootp( logistic_coefficients_var_allBoots, 0.0, nb, 'above' );
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_VAR_P = logistic_coefficients_var_bootp;
logistic_coefficients_std_bootp = bootp( logistic_coefficients_std_allBoots, 0.0, nb, 'above' );
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_STD_P = logistic_coefficients_std_bootp;

logistic_coefficients_RTS_var_bootp = bootp( logistic_coefficients_RTS_var_allBoots, 0.0, nb, 'below' );
logistic_coefficients_rts_VAR_bootp = bootp( logistic_coefficients_rts_VAR_allBoots, 0.0, nb, 'above' );
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_P = ...
        [logistic_coefficients_RTS_var_bootp, logistic_coefficients_rts_VAR_bootp];

logistic_coefficients_RTS_var_interact_bootp = bootp( logistic_coefficients_RTS_var_interact_allBoots, 0.0, nb, 'below' );
logistic_coefficients_rts_VAR_interact_bootp = bootp( logistic_coefficients_rts_VAR_interact_allBoots, 0.0, nb, 'above' );
logistic_coefficients_rts_var_INTERACT_bootp = bootp( logistic_coefficients_rts_var_INTERACT_allBoots, 0.0, nb, 'below' );
    EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_P = ...
        [logistic_coefficients_RTS_var_interact_bootp, logistic_coefficients_rts_VAR_interact_bootp, logistic_coefficients_rts_var_INTERACT_bootp];


% Print outputs.
% @todo: this section was outputting broken numbers(?)
fprintf(outputLog, ['MEMORY JUDGMENT SENSITIVITY (FREQUENT): ', ...
    num2str(EXPT1_MEMORY_BOOT.APRIME_FREQUENT_MEAN), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.APRIME_FREQUENT_CI(1)), ', ', num2str(EXPT1_MEMORY_BOOT.APRIME_FREQUENT_CI(2)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.APRIME_FREQUENT_P), ...
    '\n']);

fprintf(outputLog, ['MEMORY JUDGMENT SENSITIVITY (INFREQUENT): ', ...
    num2str(EXPT1_MEMORY_BOOT.APRIME_INFREQUENT_MEAN), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.APRIME_INFREQUENT_CI(1)), ', ', num2str(EXPT1_MEMORY_BOOT.APRIME_INFREQUENT_CI(2)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.APRIME_INFREQUENT_P), ...
    '\n']);

fprintf(outputLog, ['MEMORY JUDGMENT SENSITIVITY (DIFFERENCE BETWEEN INFREQUENT & FREQUENT): ', ...
    num2str(EXPT1_MEMORY_BOOT.APRIME_DIFFERENCE_MEAN), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.APRIME_DIFFERENCE_CI(1)), ', ', num2str(EXPT1_MEMORY_BOOT.APRIME_DIFFERENCE_CI(2)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.APRIME_DIFFERENCE_P), ...
    '\n']);
% @todo: in the original python script, the mean seems to be calculated not from the bootstrap figures?
fprintf(outputLog, ['PROPORTION REMEMBERED (INFREQUENT): ', ...
    num2str(EXPT1_MEMORY_BOOT.REMEMBERED_INFREQUENT_MEAN), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.REMEMBERED_INFREQUENT_CI(1)), ', ', num2str(EXPT1_MEMORY_BOOT.REMEMBERED_INFREQUENT_CI(2)), ')', ...
    '\n']);

fprintf(outputLog, ['PROPORTION FORGOTTEN (INFREQUENT): ', ...
    num2str(EXPT1_MEMORY_BOOT.FORGOTTEN_INFREQUENT_MEAN), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.FORGOTTEN_INFREQUENT_CI(1)), ', ', num2str(EXPT1_MEMORY_BOOT.FORGOTTEN_INFREQUENT_CI(2)), ')', ...
    '\n']);

fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT (RT): ', ...
    num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_MEAN), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_CI(1)), ', ', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_CI(2)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_P), ...
    '\n']);

fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT (VAR): ', ...
    num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_VAR_MEAN), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_VAR_CI(1)), ', ', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_VAR_CI(2)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_VAR_P), ...
    '\n']);

fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT (STD): ', ...
    num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_STD_MEAN), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_STD_CI(1)), ', ', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_STD_CI(2)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_STD_P), ...
    '\n']);

fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT, RT+VAR MODEL (RT): ', ...
    num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_MEAN(1)), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_CI(1)), ', ', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_CI(2)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_P(1)), ...
    '\n']);

fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT, RT+VAR MODEL (VAR): ', ...
    num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_MEAN(2)), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_CI(3)), ', ', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_CI(4)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_P(2)), ...
    '\n']);

fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT, RT+VAR+INTERACT MODEL (RT): ', ...
    num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_MEAN(1)), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_CI(1)), ', ', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_CI(2)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_P(1)), ...
    '\n']);

fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT, RT+VAR+INTERACT MODEL (VAR): ', ...
    num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_MEAN(2)), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_CI(3)), ', ', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_CI(4)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_P(2)), ...
    '\n']);

fprintf(outputLog, ['MEAN LOGISTIC COEFFICIENT, RT+VAR+INTERACT MODEL (INTERACT): ', ...
    num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_MEAN(3)), ...
    ' (', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_CI(5)), ', ', num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_CI(6)), ')', ...
    '\n    ', 'p-value = ' , num2str(EXPT1_MEMORY_BOOT.LOGISTIC_COEFFICIENT_RTS_VAR_INTERACT_P(3)), ...
    '\n']);





%% Unfinished/Extra Things

% see if our alternative metric for RT variability is any good for telling apart acc from inacc trials on the basis of RT
% (this gets all "significant" p-values from each subject)
std_accVsInacc_ttest2_pValue(std_accVsInacc_ttest2_pValue<0.05)

% compare amount of overlapping triggers for the two triggering types
% convert 0 to nan and -1 to 1 so we have an easy way to distinguish non-trigger trials from trigger trials
successOnly_triggers = postdiction_trigger;
successOnly_triggers(postdiction_trigger==0) = nan;
successOnly_triggers = abs(successOnly_triggers); % now, trigger == 1, and no_trigger == nan
successOnly_triggersClassic = postdiction_triggerClassic;
successOnly_triggersClassic(postdiction_triggerClassic==0) = nan;
successOnly_triggersClassic = abs(successOnly_triggersClassic);
whichOverlapping = successOnly_triggers==1 & successOnly_triggersClassic==1;

%countOverlapping_all = sum(whichOverlapping,'all');
countOverlapping_all = sum(sum(whichOverlapping));
%countVTC_all = sum(successOnly_triggers==1,'all');
countVTC_all = sum(sum(successOnly_triggers==1));
overlapRate_all = countOverlapping_all/countVTC_all;

countOverlapping_subject = sum(whichOverlapping,2);
countVTC_subject = sum(successOnly_triggers==1,2);
overlapRate_subject = countOverlapping_subject./countVTC_subject;
min(overlapRate_subject)
max(overlapRate_subject)
range(overlapRate_subject)
mean(overlapRate_subject)
% now compare overlap between each metric and chance
for iSubject=1:numSubjects
    random_overlap_vtc(iSubject,:,:) = successOnly_triggers(iSubject,:)==1 & postdiction_triggerRandom(iSubject,:,:)==1;
    random_overlap_rts(iSubject,:,:) = successOnly_triggersClassic(iSubject,:)==1 & postdiction_triggerRandom(iSubject,:,:)==1;
end
countRTS_subject = sum(successOnly_triggersClassic==1,2);
countOverlapping_subject_random_vtc = sum(random_overlap_vtc,2);
countOverlapping_subject_random_rts = sum(random_overlap_rts,2);
overlapRate_subject_random_vtc = countOverlapping_subject_random_vtc ./ countVTC_subject;
overlapRate_subject_random_rts = countOverlapping_subject_random_rts ./ countRTS_subject;
overlapRate_subject_random_vtc = mean(overlapRate_subject_random_vtc,3);
overlapRate_subject_random_rts = mean(overlapRate_subject_random_rts,3);
%mean(overlapRate_subject_random_vtc)
%mean(overlapRate_subject_random_rts)
% is the amount of overlap between rt-triggered and vtc-triggered trials (across subjects) significantly different from just triggering at random?
%ttest(overlapRate_subject,overlapRate_subject_random_vtc)
%ttest(overlapRate_subject,overlapRate_subject_random_rts)

% probably to do next is to plot a representative subject showing at which trials one trigger type would have kicked in rather than another
% helps to show the spatial distribution of trigger types (maybe all of one trigger type happened at the beginning???)

successAlt_triggers = postdiction_trigger;
successAlt_triggers(postdiction_trigger==0) = nan;
successAlt_triggersClassic = postdiction_triggerClassic;
successAlt_triggersClassic(postdiction_triggerClassic==0) = nan;
whichOverlapping_above = successAlt_triggers==1 & successAlt_triggersClassic==1;
whichOverlapping_below = successAlt_triggers==-1 & successAlt_triggersClassic==-1;
whichOverlapping_contra = (successAlt_triggers==1 & successAlt_triggersClassic==-1) | (successAlt_triggers==-1 & successAlt_triggersClassic==1);

%countOverlapping_above = sum(whichOverlapping_above,'all');
countOverlapping_above = sum(sum(whichOverlapping_above));
%countOverlapping_below = sum(whichOverlapping_below,'all');
countOverlapping_below = sum(sum(whichOverlapping_below));
%countOverlapping_contra = sum(whichOverlapping_contra,'all');
countOverlapping_contra = sum(sum(whichOverlapping_contra));
%countVTC_above = sum(successAlt_triggers==1,'all');
countVTC_above = sum(sum(successAlt_triggers==1));
%countVTC_below = sum(successAlt_triggers==-1,'all');
countVTC_below = sum(sum(successAlt_triggers==-1));
% these are inset
overlapRate_above = countOverlapping_above / countVTC_above;
overlapRate_below = countOverlapping_below / countVTC_below;
% not inset
overlapRate_contra = countOverlapping_contra / countVTC_all;

countOverlapping_subjects_above = sum(whichOverlapping_above,2);
countOverlapping_subjects_below = sum(whichOverlapping_below,2);
countVTC_subjects_above = sum(successAlt_triggers==1,2);
countVTC_subjects_below = sum(successAlt_triggers==-1,2);
overlapRate_subjects_above = countOverlapping_subjects_above ./ countVTC_subjects_above;
overlapRate_subjects_below = countOverlapping_subjects_below ./ countVTC_subjects_below;
[min(overlapRate_subjects_above), max(overlapRate_subjects_above), range(overlapRate_subjects_above), mean(overlapRate_subjects_above)]
[min(overlapRate_subjects_below), max(overlapRate_subjects_below), range(overlapRate_subjects_below), mean(overlapRate_subjects_below)]


% log models on other things (like single infrequent trial, or frequent trials)
% @todo: these are currently spitting out problematic numbers and/or warnings


% cohen's d
cohenD_rts     = ( nanmean(nanmean(nullBetas_rts,2)) - nanmean(logistic_coefficients_rts(:,2)) ) / ( (nanvar(nanmean(nullBetas_rts,2)) + nanvar(logistic_coefficients_rts(:,2)))/2 );
cohenD_var     = ( nanmean(nanmean(nullBetas_var,2)) - nanmean(logistic_coefficients_var(:,2)) ) / ( (nanvar(nanmean(nullBetas_var,2)) + nanvar(logistic_coefficients_var(:,2)))/2 );
cohenD_RTS_var = ( nanmean(nanmean(nullBetas_RTS_var,2)) - nanmean(logistic_coefficients_rts_var(:,2)) ) / ( (nanvar(nanmean(nullBetas_RTS_var,2)) + nanvar(logistic_coefficients_rts_var(:,2)))/2 );
cohenD_rts_VAR = ( nanmean(nanmean(nullBetas_rts_VAR,2)) - nanmean(logistic_coefficients_rts_var(:,3)) ) / ( (nanvar(nanmean(nullBetas_rts_VAR,2)) + nanvar(logistic_coefficients_rts_var(:,3)))/2 );
disp('cohen d stuff');
disp(num2str(cohenD_rts));
disp(num2str(cohenD_var));
disp(num2str(cohenD_RTS_var));
disp(num2str(cohenD_rts_VAR));

% logistic regression bootstrapping but for the immediately preceding trial
% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_rts_trimmed  = rmmissing(logistic_coefficients_rts( :,2 ));
i1_logistic_coefficients_rts_trimmed=i1_logistic_coefficients_rts(:,2);
i1_logistic_coefficients_rts_booted   = i1_logistic_coefficients_rts_trimmed( bs );
i1_logistic_coefficients_rts_allBoots = nanmean( i1_logistic_coefficients_rts_booted, 2 );
i1_logistic_coefficients_rts_avgBoot  = mean( i1_logistic_coefficients_rts_allBoots );
i1_logistic_coefficients_rts_bootCI   = prctile( i1_logistic_coefficients_rts_allBoots, CI95 );

% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_var_trimmed  = rmmissing(logistic_coefficients_var( :,2 ));
i1_logistic_coefficients_var_trimmed=i1_logistic_coefficients_var(:,2);
i1_logistic_coefficients_var_booted   = i1_logistic_coefficients_var_trimmed( bs );
i1_logistic_coefficients_var_allBoots = nanmean( i1_logistic_coefficients_var_booted, 2 );
i1_logistic_coefficients_var_avgBoot  = mean( i1_logistic_coefficients_var_allBoots );
i1_logistic_coefficients_var_bootCI   = prctile( i1_logistic_coefficients_var_allBoots, CI95 );

% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_RTS_var_trimmed  = rmmissing(logistic_coefficients_rts_var( :,2 ));
%logistic_coefficients_rts_VAR_trimmed  = rmmissing(logistic_coefficients_rts_var( :,3 ));
i1_logistic_coefficients_RTS_var_trimmed=i1_logistic_coefficients_rts_var(:,2);
i1_logistic_coefficients_rts_VAR_trimmed=i1_logistic_coefficients_rts_var(:,3);
i1_logistic_coefficients_RTS_var_booted   = i1_logistic_coefficients_RTS_var_trimmed( bs );
i1_logistic_coefficients_rts_VAR_booted   = i1_logistic_coefficients_rts_VAR_trimmed( bs );
i1_logistic_coefficients_RTS_var_allBoots = nanmean( i1_logistic_coefficients_RTS_var_booted, 2 );
i1_logistic_coefficients_rts_VAR_allBoots = nanmean( i1_logistic_coefficients_rts_VAR_booted, 2 );
i1_logistic_coefficients_RTS_var_avgBoot  = mean( i1_logistic_coefficients_RTS_var_allBoots );
i1_logistic_coefficients_rts_VAR_avgBoot  = mean( i1_logistic_coefficients_rts_VAR_allBoots );
i1_logistic_coefficients_RTS_var_bootCI   = prctile( i1_logistic_coefficients_RTS_var_allBoots, CI95 );
i1_logistic_coefficients_rts_VAR_bootCI   = prctile( i1_logistic_coefficients_rts_VAR_allBoots, CI95 );

% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_RTS_var_interact_trimmed  = rmmissing(logistic_coefficients_rts_var_interact( :,2 ));
%logistic_coefficients_rts_VAR_interact_trimmed  = rmmissing(logistic_coefficients_rts_var_interact( :,3 ));
%logistic_coefficients_rts_var_INTERACT_trimmed  = rmmissing(logistic_coefficients_rts_var_interact( :,4 ));
i1_logistic_coefficients_RTS_var_interact_trimmed=i1_logistic_coefficients_rts_var_interact(:,2);
i1_logistic_coefficients_rts_VAR_interact_trimmed=i1_logistic_coefficients_rts_var_interact(:,3);
i1_logistic_coefficients_rts_var_INTERACT_trimmed=i1_logistic_coefficients_rts_var_interact(:,4);
i1_logistic_coefficients_RTS_var_interact_booted   = i1_logistic_coefficients_RTS_var_interact_trimmed( bs );
i1_logistic_coefficients_rts_VAR_interact_booted   = i1_logistic_coefficients_rts_VAR_interact_trimmed( bs );
i1_logistic_coefficients_rts_var_INTERACT_booted   = i1_logistic_coefficients_rts_var_INTERACT_trimmed( bs );
i1_logistic_coefficients_RTS_var_interact_allBoots = nanmean( i1_logistic_coefficients_RTS_var_interact_booted, 2 );
i1_logistic_coefficients_rts_VAR_interact_allBoots = nanmean( i1_logistic_coefficients_rts_VAR_interact_booted, 2 );
i1_logistic_coefficients_rts_var_INTERACT_allBoots = nanmean( i1_logistic_coefficients_rts_var_INTERACT_booted, 2 );
i1_logistic_coefficients_RTS_var_interact_avgBoot  = mean( i1_logistic_coefficients_RTS_var_interact_allBoots );
i1_logistic_coefficients_rts_VAR_interact_avgBoot  = mean( i1_logistic_coefficients_rts_VAR_interact_allBoots );
i1_logistic_coefficients_rts_var_INTERACT_avgBoot  = mean( i1_logistic_coefficients_rts_var_INTERACT_allBoots );
i1_logistic_coefficients_RTS_var_interact_bootCI   = prctile( i1_logistic_coefficients_RTS_var_interact_allBoots, CI95 );
i1_logistic_coefficients_rts_VAR_interact_bootCI   = prctile( i1_logistic_coefficients_rts_VAR_interact_allBoots, CI95 );
i1_logistic_coefficients_rts_var_INTERACT_bootCI   = prctile( i1_logistic_coefficients_rts_var_INTERACT_allBoots, CI95 );

% Are the beta coefficients significantly different from 0?
i1_logistic_coefficients_rts_bootp = bootp( i1_logistic_coefficients_rts_allBoots, 0.0, nb, 'below' );
i1_logistic_coefficients_var_bootp = bootp( i1_logistic_coefficients_var_allBoots, 0.0, nb, 'above' );

i1_logistic_coefficients_RTS_var_bootp = bootp( i1_logistic_coefficients_RTS_var_allBoots, 0.0, nb, 'below' );
i1_logistic_coefficients_rts_VAR_bootp = bootp( i1_logistic_coefficients_rts_VAR_allBoots, 0.0, nb, 'above' );

i1_logistic_coefficients_RTS_var_interact_bootp = bootp( i1_logistic_coefficients_RTS_var_interact_allBoots, 0.0, nb, 'below' );
i1_logistic_coefficients_rts_VAR_interact_bootp = bootp( i1_logistic_coefficients_rts_VAR_interact_allBoots, 0.0, nb, 'above' );
i1_logistic_coefficients_rts_var_INTERACT_bootp = bootp( i1_logistic_coefficients_rts_var_INTERACT_allBoots, 0.0, nb, 'below' );

% logistic regression bootstrapping but for frequent trials
% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_rts_trimmed  = rmmissing(logistic_coefficients_rts( :,2 ));
f_logistic_coefficients_rts_trimmed=f_logistic_coefficients_rts(:,2);
f_logistic_coefficients_rts_booted   = f_logistic_coefficients_rts_trimmed( bs );
f_logistic_coefficients_rts_allBoots = nanmean( f_logistic_coefficients_rts_booted, 2 );
f_logistic_coefficients_rts_avgBoot  = mean( f_logistic_coefficients_rts_allBoots );
f_logistic_coefficients_rts_bootCI   = prctile( f_logistic_coefficients_rts_allBoots, CI95 );

% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_var_trimmed  = rmmissing(logistic_coefficients_var( :,2 ));
f_logistic_coefficients_var_trimmed=f_logistic_coefficients_var(:,2);
f_logistic_coefficients_var_booted   = f_logistic_coefficients_var_trimmed( bs );
f_logistic_coefficients_var_allBoots = nanmean( f_logistic_coefficients_var_booted, 2 );
f_logistic_coefficients_var_avgBoot  = mean( f_logistic_coefficients_var_allBoots );
f_logistic_coefficients_var_bootCI   = prctile( f_logistic_coefficients_var_allBoots, CI95 );

% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_RTS_var_trimmed  = rmmissing(logistic_coefficients_rts_var( :,2 ));
%logistic_coefficients_rts_VAR_trimmed  = rmmissing(logistic_coefficients_rts_var( :,3 ));
f_logistic_coefficients_RTS_var_trimmed=f_logistic_coefficients_rts_var(:,2);
f_logistic_coefficients_rts_VAR_trimmed=f_logistic_coefficients_rts_var(:,3);
f_logistic_coefficients_RTS_var_booted   = f_logistic_coefficients_RTS_var_trimmed( bs );
f_logistic_coefficients_rts_VAR_booted   = f_logistic_coefficients_rts_VAR_trimmed( bs );
f_logistic_coefficients_RTS_var_allBoots = nanmean( f_logistic_coefficients_RTS_var_booted, 2 );
f_logistic_coefficients_rts_VAR_allBoots = nanmean( f_logistic_coefficients_rts_VAR_booted, 2 );
f_logistic_coefficients_RTS_var_avgBoot  = mean( f_logistic_coefficients_RTS_var_allBoots );
f_logistic_coefficients_rts_VAR_avgBoot  = mean( f_logistic_coefficients_rts_VAR_allBoots );
f_logistic_coefficients_RTS_var_bootCI   = prctile( f_logistic_coefficients_RTS_var_allBoots, CI95 );
f_logistic_coefficients_rts_VAR_bootCI   = prctile( f_logistic_coefficients_rts_VAR_allBoots, CI95 );

% @todo: figure out what to do with subjects that have a nan value; include them or leave them out of bootstrapping?
%logistic_coefficients_RTS_var_interact_trimmed  = rmmissing(logistic_coefficients_rts_var_interact( :,2 ));
%logistic_coefficients_rts_VAR_interact_trimmed  = rmmissing(logistic_coefficients_rts_var_interact( :,3 ));
%logistic_coefficients_rts_var_INTERACT_trimmed  = rmmissing(logistic_coefficients_rts_var_interact( :,4 ));
f_logistic_coefficients_RTS_var_interact_trimmed=f_logistic_coefficients_rts_var_interact(:,2);
f_logistic_coefficients_rts_VAR_interact_trimmed=f_logistic_coefficients_rts_var_interact(:,3);
f_logistic_coefficients_rts_var_INTERACT_trimmed=f_logistic_coefficients_rts_var_interact(:,4);
f_logistic_coefficients_RTS_var_interact_booted   = f_logistic_coefficients_RTS_var_interact_trimmed( bs );
f_logistic_coefficients_rts_VAR_interact_booted   = f_logistic_coefficients_rts_VAR_interact_trimmed( bs );
f_logistic_coefficients_rts_var_INTERACT_booted   = f_logistic_coefficients_rts_var_INTERACT_trimmed( bs );
f_logistic_coefficients_RTS_var_interact_allBoots = nanmean( f_logistic_coefficients_RTS_var_interact_booted, 2 );
f_logistic_coefficients_rts_VAR_interact_allBoots = nanmean( f_logistic_coefficients_rts_VAR_interact_booted, 2 );
f_logistic_coefficients_rts_var_INTERACT_allBoots = nanmean( f_logistic_coefficients_rts_var_INTERACT_booted, 2 );
f_logistic_coefficients_RTS_var_interact_avgBoot  = mean( f_logistic_coefficients_RTS_var_interact_allBoots );
f_logistic_coefficients_rts_VAR_interact_avgBoot  = mean( f_logistic_coefficients_rts_VAR_interact_allBoots );
f_logistic_coefficients_rts_var_INTERACT_avgBoot  = mean( f_logistic_coefficients_rts_var_INTERACT_allBoots );
f_logistic_coefficients_RTS_var_interact_bootCI   = prctile( f_logistic_coefficients_RTS_var_interact_allBoots, CI95 );
f_logistic_coefficients_rts_VAR_interact_bootCI   = prctile( f_logistic_coefficients_rts_VAR_interact_allBoots, CI95 );
f_logistic_coefficients_rts_var_INTERACT_bootCI   = prctile( f_logistic_coefficients_rts_var_INTERACT_allBoots, CI95 );

% Are the beta coefficients significantly different from 0?
f_logistic_coefficients_rts_bootp = bootp( f_logistic_coefficients_rts_allBoots, 0.0, nb, 'below' );
f_logistic_coefficients_var_bootp = bootp( f_logistic_coefficients_var_allBoots, 0.0, nb, 'above' );

f_logistic_coefficients_RTS_var_bootp = bootp( f_logistic_coefficients_RTS_var_allBoots, 0.0, nb, 'below' );
f_logistic_coefficients_rts_VAR_bootp = bootp( f_logistic_coefficients_rts_VAR_allBoots, 0.0, nb, 'above' );

f_logistic_coefficients_RTS_var_interact_bootp = bootp( f_logistic_coefficients_RTS_var_interact_allBoots, 0.0, nb, 'below' );
f_logistic_coefficients_rts_VAR_interact_bootp = bootp( f_logistic_coefficients_rts_VAR_interact_allBoots, 0.0, nb, 'above' );
f_logistic_coefficients_rts_var_INTERACT_bootp = bootp( f_logistic_coefficients_rts_var_INTERACT_allBoots, 0.0, nb, 'below' );

% print the things
%
disp('i1, logistic coefficient, rts');
disp(i1_logistic_coefficients_rts_avgBoot);
disp(i1_logistic_coefficients_rts_bootCI );

disp('i1, logistic coefficient, var');
disp(i1_logistic_coefficients_var_avgBoot);
disp(i1_logistic_coefficients_var_bootCI );

disp('i1, logistic coefficient, RTS+var');
disp(i1_logistic_coefficients_RTS_var_avgBoot);
disp(i1_logistic_coefficients_RTS_var_bootCI );
disp('i1, logistic coefficient, rts+VAR');
disp(i1_logistic_coefficients_rts_VAR_avgBoot);
disp(i1_logistic_coefficients_rts_VAR_bootCI );

disp('i1, logistic coefficient, RTS+var+interact');
disp(i1_logistic_coefficients_RTS_var_interact_avgBoot);
disp(i1_logistic_coefficients_RTS_var_interact_bootCI );
disp('i1, logistic coefficient, rts+VAR+interact');
disp(i1_logistic_coefficients_rts_VAR_interact_avgBoot);
disp(i1_logistic_coefficients_rts_VAR_interact_bootCI );
disp('i1, logistic coefficient, rts+var+INTERACT');
disp(i1_logistic_coefficients_rts_var_INTERACT_avgBoot);
disp(i1_logistic_coefficients_rts_var_INTERACT_bootCI );

% Are the beta coefficients significantly different from 0?
disp('i1, logistic coefficient, nonparametric p, rts');
disp(i1_logistic_coefficients_rts_bootp);
[~,asdfp,~,~] = ttest(i1_logistic_coefficients_rts(:,2))
disp('i1, logistic coefficient, nonparametric p, var');
disp(i1_logistic_coefficients_var_bootp);
[~,asdfp,~,~] = ttest(i1_logistic_coefficients_var(:,2))

disp('i1, logistic coefficient, nonparametric p, RTS+var');
disp(i1_logistic_coefficients_RTS_var_bootp);
[~,asdfp,~,~] = ttest(i1_logistic_coefficients_rts_var(:,2))
disp('i1, logistic coefficient, nonparametric p, rts+VAR');
disp(i1_logistic_coefficients_rts_VAR_bootp);
[~,asdfp,~,~] = ttest(i1_logistic_coefficients_rts_var(:,3))

disp('i1, logistic coefficient, nonparametric p, RTS+var+interact');
disp(i1_logistic_coefficients_RTS_var_interact_bootp);
[~,asdfp,~,~] = ttest(i1_logistic_coefficients_rts_var_interact(:,2))
disp('i1, logistic coefficient, nonparametric p, rts+VAR+interact');
disp(i1_logistic_coefficients_rts_VAR_interact_bootp);
[~,asdfp,~,~] = ttest(i1_logistic_coefficients_rts_var_interact(:,3))
disp('i1, logistic coefficient, nonparametric p, rts+var+INTERACT');
disp(i1_logistic_coefficients_rts_var_INTERACT_bootp);
[~,asdfp,~,~] = ttest(i1_logistic_coefficients_rts_var_interact(:,4))

%
disp('frequent, logistic coefficient, rts');
disp(f_logistic_coefficients_rts_avgBoot);
disp(f_logistic_coefficients_rts_bootCI );

disp('frequent, logistic coefficient, var');
disp(f_logistic_coefficients_var_avgBoot);
disp(f_logistic_coefficients_var_bootCI );

disp('frequent, logistic coefficient, RTS+var');
disp(f_logistic_coefficients_RTS_var_avgBoot);
disp(f_logistic_coefficients_RTS_var_bootCI );
disp('frequent, logistic coefficient, rts+VAR');
disp(f_logistic_coefficients_rts_VAR_avgBoot);
disp(f_logistic_coefficients_rts_VAR_bootCI );

disp('frequent, logistic coefficient, RTS+var+interact');
disp(f_logistic_coefficients_RTS_var_interact_avgBoot);
disp(f_logistic_coefficients_RTS_var_interact_bootCI );
disp('frequent, logistic coefficient, rts+VAR+interact');
disp(f_logistic_coefficients_rts_VAR_interact_avgBoot);
disp(f_logistic_coefficients_rts_VAR_interact_bootCI );
disp('frequent, logistic coefficient, rts+var+INTERACT');
disp(f_logistic_coefficients_rts_var_INTERACT_avgBoot);
disp(f_logistic_coefficients_rts_var_INTERACT_bootCI );

% Are the beta coefficients significantly different from 0?
disp('frequent, logistic coefficient, nonparametric p, rts');
disp(f_logistic_coefficients_rts_bootp);
[~,asdfp,~,~] = ttest(f_logistic_coefficients_rts(:,2))
disp('frequent, logistic coefficient, nonparametric p, var');
disp(f_logistic_coefficients_var_bootp);
[~,asdfp,~,~] = ttest(f_logistic_coefficients_var(:,2))

disp('frequent, logistic coefficient, nonparametric p, RTS+var');
disp(f_logistic_coefficients_RTS_var_bootp);
[~,asdfp,~,~] = ttest(f_logistic_coefficients_rts_var(:,2))
disp('frequent, logistic coefficient, nonparametric p, rts+VAR');
disp(f_logistic_coefficients_rts_VAR_bootp);
[~,asdfp,~,~] = ttest(f_logistic_coefficients_rts_var(:,3))

disp('frequent, logistic coefficient, nonparametric p, RTS+var+interact');
disp(f_logistic_coefficients_RTS_var_interact_bootp);
[~,asdfp,~,~] = ttest(f_logistic_coefficients_rts_var_interact(:,2))
disp('frequent, logistic coefficient, nonparametric p, rts+VAR+interact');
disp(f_logistic_coefficients_rts_VAR_interact_bootp);
[~,asdfp,~,~] = ttest(f_logistic_coefficients_rts_var_interact(:,3))
disp('frequent, logistic coefficient, nonparametric p, rts+var+INTERACT');
disp(f_logistic_coefficients_rts_var_INTERACT_bootp);
[~,asdfp,~,~] = ttest(f_logistic_coefficients_rts_var_interact(:,4))





%% #SAVE - Save outputs
% NOTE THAT THIS WILL OVERRIDE DUPLICATES (i.e. files with the same file name)

if mFlags('SAVE_PLOTS') == true
    
    % groot is a hardcoded MATLAB variable (like gcf and gca)
    graphicsHandler = groot;
    % Get the total number of figures that are open...
    numFigures = length(graphicsHandler.Children);
    % ...so we can cycle through each and save all of them.
    for i=1:numFigures
        currentFigure = figure(i);
        saveas( currentFigure , [projectDirectory_plots, currentFigure.Name] , 'png' );
        saveas( currentFigure , [projectDirectory_plots, currentFigure.Name] , 'fig' );
    end
    
    % Bam, all the outputted figures should be automagically saved to the 'figures' folder.
    
    % If we are going to be generating more figures after this, then clear out the old figures that we just saved.
    % This avoids the script from having to save all the old figures again on top of the new figures for each iteration.
    if mFlags('MASTER_TEST') == true
        close all;
    end
    
end
% Clear all graphics-related variables we've instantiated in the workspace.
% (If a graphics variable is saved and then loaded again, it will also generate the figure again, which is unnecessary.)
clear graphicsHandler
clear currentFigure

if mFlags('SAVE_STATS') == true
    
    % Save the workspace to a .mat file. Make sure to add all of the relevant flag-based prefixes.
    filename = [saveSetupPrefix, 'SustainedAttentionExperiment.mat'];
    save( [projectDirectory_stats, filename] );
    
    % Bam, all variables in the workspace should be saved to a "snapshot" now.
    
end


end


