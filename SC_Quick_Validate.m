% This script quickly visualises and validates raw data as a sanity check.
clearvars; clc;
numberOfSubjects = 34;

%% Variable Outline

% attnData Variables:
%
% trialsPerRun       (1  ) = integer; number of trials (default = 500)
% trial              (500) = integer; index list (e.g. [1, 2, 3, 4, ...])
% plannedTrialOnsets (500) = double; timestamps, calculated in advance (i.e. spaced exactly 1s apart)
% actualTrialOnsets  (500) = double; timestamps, measured from screen flips (roughly equal to plannedTrialOnsets)
% rts                (500) = double; ranges from 0.0s to 1.0s
% accs               (500) = boolean; 1 is correct, 0 is incorrect
% resps              (500) = integer; input keyCode (11='h', 13='j')
% corrresps          (500) = integer; hypothetical input keyCode (11='h', 13='j')
% categs             (500) = integer; image classification category (OUTDOOR=1, INDOOR=2)
% images             (500) = integer; ID of the image, ranging from 1 to ???
% imageName          {500} = char; filename of the image

% memData Variables:
%
% trial              (200) = integer; index list (e.g. [1, 2, 3, 4, ...])
% attCategs          (200) = integer; frequency-type of image from attention task (1 or 2 == prepotent or counter == indoor or outdoor)
% memCategs          (200) = integer; "punnett square" classification between both indoor/outdoor and old/new (1 and 3 == prepotent; 1 and 2 == old; etc)
% images             (200) = integer; ID of the image, ranging from 1 to ???
% attOrder           (200) = integer/nan; trial number of the image from previous task (applicable only to old images, otherwise == nan)
% rts                (200) = double; no range limit
% resps              (200) = integer; input keyCode (30='1!', 31='2@', 32='3#', 33='4$')
% rating             (200) = integer; ranges from 1 to 4 (1 == definitely new, 2 == maybe new, 3 == maybe old, 4 = definitely old)
% accs               (200) = boolean; 1 is correct, 0 is incorrect (degree of confidence doesn't matter)
% actualTrialOnsets  (200) = double; timestamps, measured from screen flips
% actualTrialOffsets (200) = double; timestamps, measured from screen flips


%% Validation Criteria

% attnData:
% 
% trial IDs: there are 500 trials, and the list is of well-ordered integers (i.e. monotonically increasing)
% flip timing accuracy: difference between planned and actual offsets is never above a tolerance threshold (if it is, we have a low-level confounding variable)
% reaction times: values must be between 0.0 and 1.0
% accuracy: values must be 0 or 1
% response inputs: values must be 11 or 13
% correct response inputs: values must be 11 or 13

% memData Variables:
% 
% 


%% Validation Methods

% attnData Variables:
% 
% trial: length(trial) == 500; min(trial) == 1; max(trial) == 500

% memData Variables:
% 
% 


%% Validation Conditions

% i.e. what should the graphs look like


%% Implementation
% Validations are assumed to be the **default** case, meaning that a successful validation is *silent*.
% Only failure to validate will output a message.
% NOTE: Validation procedures do not factor in NaN values - please keep that in mind. Visualisation can help mitigate that pitfall.

% get the full path name of this script
projectDirectory = mfilename('fullpath');
% however, we want the path of the folder containing this script, not the path of the script itself
projectDirectory_decomposed = split(projectDirectory, '/'); % split the path...
projectDirectory_scriptFilename = projectDirectory_decomposed{end}; % ...and extract this script's filename...
projectDirectory_scriptFilename_length = length(projectDirectory_scriptFilename); % ...and get its length, so that we can...
projectDirectory = projectDirectory(1 : end - projectDirectory_scriptFilename_length); % ...extract the filename via indexing
% bam, now we have our project directory


% specify our default expectations
defaultTrials_attention = 500;
attentionOnsetThreshold = 0.100; % in seconds
keycode_h = 11;
keycode_j = 13;
defaultTrials_memory = 200;
proportionOld = round(0.50 * defaultTrials_memory);
proportionNew = round(0.50 * defaultTrials_memory);
keycode_1 = 30;
keycode_2 = 31;
keycode_3 = 32;
keycode_4 = 33;

% probably should refactor this terrible visualisation flag code
% viz = visualisation, a = attention, m = memory, # = section number of the validation code
viza1=0;
viza2=0;
viza3=0;
viza4=0;
viza5=0;
viza6=0;
viza7=0;
vizm1=0;
vizm2=0;
vizm3=0;
vizm4=0;
vizm5=0;
vizm6=0;
vizm7=0;
vizm8=0;


subjectID_start = 1;
subjectID_end = numberOfSubjects;

% loop through all the subject IDs, and validate, 1-by-1; this assumes that *no subject IDs were removed*
for subjectID = subjectID_start:subjectID_end
    
    % Subject-Specific Data Setup
    subjectID_displayString = [ 'Subject ID: ', num2str(subjectID), ', ' ];
    filename_attnData = dir([projectDirectory, 'data/', num2str(subjectID), '/', 'attnData_*']);
    filename_memData = dir([projectDirectory, 'data/', num2str(subjectID), '/', 'memData_*']);
    % attnData = ...
    % memData = ...
    load([filename_attnData.folder, '/',  filename_attnData.name]); % loads a variable with the name 'attnData'
    load([filename_memData.folder,  '/',  filename_memData.name]); % loads a variable with the name 'memData'
    
    
    
    % Validate Attention Data (i.e. part 1 of the study)
    
    % Verify trial IDs
    trials = rmmissing( attnData.trial );
    if (length(trials) == defaultTrials_attention) && (min(trials) == 1) && (max(trials) == defaultTrials_attention)
        % (default case)
    else
        disp( [subjectID_displayString, 'TRIAL ID VALIDATION FAILED (attention)'] );
    end
    if viza1==1; figure; hist(trials); title('trialIDs'); end
    
    % Verify trial onsets
    onsetDifferences = rmmissing( abs(attnData.plannedTrialOnsets - attnData.actualTrialOnsets) );
    if max(onsetDifferences) > attentionOnsetThreshold
        disp( [subjectID_displayString, 'TRIAL ONSET VALIDATION FAILED (attention)'] );
    else
        % (default case)
    end
    if viza2==1; figure; hist(onsetDifferences); title('onsetDifferences'); end
    
    % Verify response times
    responseTimes = rmmissing( attnData.rts );
    if (min(responseTimes) < 0.0) || (max(responseTimes) > 1.0)
        disp( [subjectID_displayString, 'RESPONSE TIME VALIDATION FAILED (attention)'] );
    else
        % (default case)
    end
    if viza3==1; figure; hist(responseTimes); title('responseTimes'); end
    
    % Verify accuracies
    accuracies = rmmissing( attnData.accs );
    if (length(unique(accuracies)) == 2) && (ismember(0, accuracies)) && (ismember(1, accuracies))
        % (default case)
    else
        disp( [subjectID_displayString, 'ACCURACY ENCODING VALIDATION FAILED (attention)'] );
    end
    if viza4==1; figure; hist(accuracies); title('accuracyOfResponses'); end
    
    % Verify subject's actual input responses
    subjectResponses = rmmissing( attnData.resps );
    if (length(unique(subjectResponses)) == 2) && (ismember(keycode_h, subjectResponses)) && (ismember(keycode_j, subjectResponses))
        % (default case)
    else
        disp( [subjectID_displayString, 'SUBJECT KEYCODE INPUT VALIDATION FAILED (attention)'] );
    end
    if viza5==1; figure; hist(subjectResponses); title('responses_keycodeInput'); end
    
    % Verify hypothetically-correct responses
    correctResponses = rmmissing( attnData.corrresps );
    if (length(unique(correctResponses)) == 2) && (ismember(keycode_h, correctResponses)) && (ismember(keycode_j, correctResponses))
        % (default case)
    else
        disp( [subjectID_displayString, 'KEYCODE ANSWER VALIDATION FAILED (attention)'] );
    end
    if viza6==1; figure; hist(correctResponses); title('correctAnswer_keycodeInput'); end
    
    % Verify category encoding of trials/pictures
    imageCategories = rmmissing( attnData.categs );
    if (length(unique(imageCategories)) == 2) && (ismember(1, imageCategories)) && (ismember(2, imageCategories))
        % (default case)
    else
        disp( [subjectID_displayString, 'TRIAL CATEGORY ENCODING VALIDATION FAILED (attention)'] );
    end
    if viza7==1; figure; hist(imageCategories); title('imageCategories'); end
    
    
    
    % Validate Memory Data (i.e. part 2 of the study)
    
    % Verify trial IDs
    trials = rmmissing( memData.trial );
    if (length(trials) == defaultTrials_memory) && (min(trials) == 1) && (max(trials) == defaultTrials_memory)
        % (default case)
    else
        disp( [subjectID_displayString, 'TRIAL ID VALIDATION FAILED (memory)'] );
    end
    if vizm1==1; figure; hist(trials); title('trialIDs'); end
    
    % Verify rarity-type category encoding of trials/pictures
    imageCategories = rmmissing( memData.attCategs );
    if (length(unique(imageCategories)) == 2) && (ismember(1, imageCategories)) && (ismember(2, imageCategories))
        % (default case)
    else
        disp( [subjectID_displayString, 'TRIAL RARITY ENCODING VALIDATION FAILED (memory)'] );
    end
    if vizm2==1; figure; hist(imageCategories); title('imageCategories_rarity'); end
    
    % Verify memory-task category encoding of trials/pictures
    imageCategories = rmmissing( memData.memCategs );
    if (length(unique(imageCategories)) == 4) && (ismember(1, imageCategories)) && (ismember(2, imageCategories)) ...
                                              && (ismember(3, imageCategories)) && (ismember(4, imageCategories))
        % (default case)
    else
        disp( [subjectID_displayString, 'TRIAL CATEGORY ENCODING VALIDATION FAILED (memory)'] );
    end
    if vizm3==1; figure; hist(imageCategories); title('imageCategories_''punnett'''); end
    
    % Verify trial IDs of trials containing old pictures, i.e. from previous task (if not from previous, then it is nan)
    trials = rmmissing( memData.attOrder );
    if (length(trials) == proportionOld) && (min(trials) >= 1) && (max(trials) <= defaultTrials_attention)
        % (default case)
    else
        disp( [subjectID_displayString, 'TRIAL ATTENTION-ID VALIDATION FAILED (memory)'] );
    end
    if vizm4==1; figure; hist(trials); title('trialIDs_old'); end
    
    % Verify response times
    responseTimes = rmmissing( memData.rts );
    if (min(responseTimes) < 0.0) %|| (max(responseTimes) > inf)
        disp( [subjectID_displayString, 'RESPONSE TIME VALIDATION FAILED (memory)'] );
    else
        % (default case)
    end
    if vizm5==1; figure; hist(responseTimes); title('responseTimes'); end
    
    % Verify subject's actual input responses
    subjectResponses = rmmissing( memData.resps );
    if (length(unique(subjectResponses)) == 4) && (ismember(keycode_1, subjectResponses)) ...
                                               && (ismember(keycode_2, subjectResponses)) ...
                                               && (ismember(keycode_3, subjectResponses)) ...
                                               && (ismember(keycode_4, subjectResponses))
        % (default case)
    else
        disp( [subjectID_displayString, 'SUBJECT KEYCODE INPUT VALIDATION FAILED (memory)'] );
    end
    if vizm6==1; figure; hist(subjectResponses); title('responses_keycodeInput'); end
    
    % Verify subject's rating response encoding
    subjectRatings = rmmissing( memData.rating );
    if (length(unique(subjectResponses)) == 4) && (ismember(1, subjectRatings)) ...
                                               && (ismember(2, subjectRatings)) ...
                                               && (ismember(3, subjectRatings)) ...
                                               && (ismember(4, subjectRatings))
        % (default case)
    else
        disp( [subjectID_displayString, 'SUBJECT RATING ENCODING VALIDATION FAILED (memory)'] );
    end
    if vizm7==1; figure; hist(subjectRatings); title('responses_ratingInput'); end
    
    % Verify accuracies
    accuracies = rmmissing( memData.accs );
    if (length(unique(accuracies)) == 2) && (ismember(0, accuracies)) && (ismember(1, accuracies))
        % (default case)
    else
        disp( [subjectID_displayString, 'ACCURACY ENCODING VALIDATION FAILED (memory)'] );
    end
    if vizm8==1; figure; hist(accuracies); title('accuracyOfResponses'); end
    
    
    
end


