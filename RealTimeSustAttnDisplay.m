function [attnData] = RealTimeSustAttnDisplay(varargin)
%
% Scene sustained attention experiment
%
%
% INPUTS:
% - subjectNum:  participant number [integer]
% - debug:       whether debugging [1 if debugging/0 if not]
%
% OUTPUTS
% - attnData:   structure of experiment design and behavioral responses
%
% Original by: Megan deBettencourt
% Modified by: Steven Cao
% Version: 1.1
% Last created: March 2014
% Last modified: January 2020

%% check inputs

% check that there is a sufficient number of inputs and that inputs are of the right format
if nargin > 2
    error('Too many inputs: only takes up to two, subjectNum and debug');
elseif nargin == 2
    if ~isnumeric(varargin{1}); error('subjectNum must be a number'); end
    if ((varargin{2}~=1) && (varargin{2}~=0)); error('debug must be either 1 (if debugging) or 0 (if not)'); end
    subjectNum = varargin{1};
    debug = varargin{2};
elseif nargin == 1
    if ~isnumeric(varargin{1}); error('subjectNum must be a number'); end
    subjectNum = varargin{1};
    debug = false;
else
    subjectNum = nan;
    debug = false;
end


%% Boilerplate

% if no arguments were provided, autofill the subjectNum so that it will be the highest existing ID
if isnan(subjectNum)
    folderInfo = dir('rtdata/');
    listSubjectEntries = [];
    for i=1:length(dir('rtdata/'))
        dataFolderName = folderInfo(i).name;
        dataFolderNumber = str2double(dataFolderName);
        if ~isempty(dataFolderNumber) && isnumeric(dataFolderNumber)
            listSubjectEntries = cat(2, listSubjectEntries, dataFolderNumber);
        end
    end
    maxSubjectID = max(listSubjectEntries);
    subjectNum = maxSubjectID + 1;
end

% make subject data folder, if it does not exist already
dataHeader = ['rtdata/' num2str(subjectNum)];
if ~exist(dataHeader, 'file'); mkdir(dataHeader); end

% random seed
seed = sum(100*clock);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));

% misc
Screen('Preference', 'SkipSyncTests', 1);
if ~debug; HideCursor; end
ListenChar(2); % suppress key inputs to matlab

% initialize system time calls
GetSecs;

%% Experimental Parameters

% main parameters
trialsPerRun = 500; % default is 500
proportionCounterCategory = 0.10;
nBufferTrials = 50; % "buffer trials" are the number of trials in the beginning, and the number of trials in the end of the run
                    % within such "buffer zones", counter trials are randomly-distributed, unlike the middle of the run which uses realtime triggering
nCounterTrials_GoodAttn = (trialsPerRun - nBufferTrials*2) * proportionCounterCategory*0.5; % default is 20
nCounterTrials_BadAttn  = (trialsPerRun - nBufferTrials*2) * proportionCounterCategory*0.5; % default is 20

% trial timing
stimDur = 1;        % secs
respWindow = .9;    % secs
IBI = 4;

% display parameters
textColor = 0;      % black text
textFont = 'Arial';
textSize = 20;
textSpacing = 25;   
fixColor = 0;       % black fixation dot before response
respColor = 255;    % white fixation dot after response
backColor = 127;    % gray background
imageSize = 256;    % assumed square 
fixationSize = 4;   % pixels
progWidth = 400;    % image loading progress bar width
progHeight = 20;    % image loading progress bar height
ScreenResX = 1600;  % TO BE CHANGED - check the resolution of the running room computers!!!
ScreenResY = 900;   % TO BE CHANGED - check the resolution of the running room computers!

nSubCategs=2;
OUTDOOR=1;
INDOOR=2;

%% Response Mapping and Counterbalancing

% skyra: use current design button box (keys 1,2,3,4)
LEFT = KbName('h');
RIGHT = KbName('j');

starttaskInstruct = ' ';


%% Initialize Screens

screenNumbers = Screen('Screens');

if debug
    %show on first screen
    screenNum = 0;
    
    %draw on only top left quarter of screen
    %[screenX screenY] = Screen('WindowSize',screenNum);
    screenX = 500;
    screenY = 500;
else
    %show on last screen (e.g. second monitor)
    screenNum = screenNumbers(end);
    
    %draw on full screen
    [screenX, screenY] = Screen('WindowSize',screenNum);
    
    %to ensure that the images are standardized (they take up the same degrees of the visual field) for all subjects
    if (screenX ~= ScreenResX) || (screenY ~= ScreenResY)
        fprintf('The screen dimensions may be incorrect. For screenNum = %d,screenX = %d (not 1152) and screenY = %d (not 864)', screenNum, screenX, screenY);
    end
end

% create main window
%mainWindow = Screen(screenNum,'OpenWindow',backColor,[0 0 screenX screenY]);
mainWindow = Screen(screenNum,'OpenWindow',backColor);

% details of main window
centerX = screenX/2; centerY = screenY/2;
Screen(mainWindow,'TextFont',textFont);
Screen(mainWindow,'TextSize',textSize);

% placeholder for images
imageRect = [0,0,imageSize,imageSize];

% position of images
centerRect = [centerX-imageSize/2,centerY-imageSize/2,centerX+imageSize/2,centerY+imageSize/2];

% position of fixation dot
fixDotRect = [centerX-fixationSize,centerY-fixationSize,centerX+fixationSize,centerY+fixationSize];

% image loading progress bar
progRect = [centerX-progWidth/2,centerY-progHeight/2,centerX+progWidth/2,centerY+progHeight/2];


%% Load Images

cd instructimages;
for categ=1:nSubCategs
    
    % move into the right folder
    if (categ == OUTDOOR)
        cd outdoor;
    elseif (categ == INDOOR)
        cd indoor;
    else
        error('Impossible category!');
    end
    
    % get filenames
    dirList{categ} = dir;  %#ok<AGROW>
    dirList{categ} = dirList{categ}(3:end); %#ok<AGROW> %skip . & ..
    if (~isempty(dirList{categ}))
        if (strcmp(dirList{categ}(1).name,'.DS_Store')==1)
            dirList{categ} = dirList{categ}(2:end);  %#ok<AGROW>
        end
        
        if (strcmp(dirList{categ}(end).name,'Thumbs.db')==1)
            dirList{categ} = dirList{categ}(1:(end-1));  %#ok<AGROW>
        end
        
        instructNumImages(categ) = length(dirList{categ}); %#ok<AGROW>
        
        if (instructNumImages(categ)>0)
            
            % get images
            for img=1:instructNumImages(categ)
                
                % update progress bar
                Screen('FrameRect',mainWindow,0,progRect,10);
                Screen('FillRect',mainWindow,0,progRect);
                Screen('FillRect',mainWindow,[255 0 0],progRect-[0 0 round((1-img/instructNumImages(categ))*progWidth) 0]);
                Screen('Flip',mainWindow);
                
                % read images
                instructImages{categ,img} = imread(dirList{categ}(img).name);  %#ok<AGROW>
            end
            
            % randomize order of images in each run
            instructImageShuffle{categ} = randperm(instructNumImages(categ)); %#ok<AGROW>
            cd ..;
        end
    else
        error('Need at least one image per directory!');
    end
end
cd ..;
Screen('Flip',mainWindow);


cd images;
for categ=1:nSubCategs
    
    % move into the right folder
    if (categ == OUTDOOR)
        cd sunoutdoor550;
    elseif (categ == INDOOR)
        cd sunindoor550;
    else
        error('Impossible category!');
    end
    
    % get filenames
    dirList{categ} = dir;  
    dirList{categ} = dirList{categ}(3:end); % skip . & ..
    if (~isempty(dirList{categ}))
        if (strcmp(dirList{categ}(1).name,'.DS_Store')==1)
            dirList{categ} = dirList{categ}(2:end);  
        end
        
        if (strcmp(dirList{categ}(end).name,'Thumbs.db')==1)
            dirList{categ} = dirList{categ}(1:(end-1));  
        end
        
        numImages(categ) = length(dirList{categ}); %#ok<AGROW>
        %numImages(categ) = 100; %FOR TESTING!!!!
        
        if (numImages(categ)>0)
            
            % get images
            for img=1:numImages(categ)
                
                % update progress bar
                Screen('FrameRect',mainWindow,0,progRect,10);
                Screen('FillRect',mainWindow,0,progRect);
                Screen('FillRect',mainWindow,[255 0 0],progRect-[0 0 round((1-img/numImages(categ))*progWidth) 0]);
                Screen('Flip',mainWindow);
                
                % read images
                imageNames{categ,img} = dirList{categ}(img).name;   %#ok<AGROW>
                images{categ,img} = imread(dirList{categ}(img).name);  %#ok<AGROW>
            end
            
            % randomize order of images in each run
            imageShuffle{categ} = randperm(numImages(categ)); %#ok<AGROW>
            cd ..;
        end
    else
        error('Need at least one image per directory!');
    end
end
cd ..;
Screen('Flip',mainWindow);


%% Output Files Setup

if subjectNum
    % open and set-up output file
    dataFile = fopen([dataHeader '/behavior.txt'],'a');
    fprintf(dataFile,'\n*********************************************\n');
    fprintf(dataFile,'* Sustained Attention and Recognition Memory Experiment - Sustained Attention Task v1.1\n');
    fprintf(dataFile,['* Date/Time: ' datestr(now,0) '\n']);
    fprintf(dataFile,['* Seed: ' num2str(seed) '\n']);
    fprintf(dataFile,['* Subject Number: ' num2str(subjectNum) '\n']);
    fprintf(dataFile,['* debug: ' num2str(debug) '\n']);
    fprintf(dataFile,'*********************************************\n\n');
end

% print header to command window
fprintf('\n*********************************************\n');
fprintf('* Sustained Attention and Recognition Memory Experiment - Sustained Attention Task v1.1\n');
fprintf(['* Date/Time: ' datestr(now,0) '\n']);
fprintf(['* Seed: ' num2str(seed) '\n']);
fprintf(['* Subject Number: ' num2str(subjectNum) '\n']);
fprintf(['* debug: ' num2str(debug) '\n']);
fprintf('*********************************************\n\n');


%% set up

if (mod(subjectNum,2)==1)
    counterCateg = INDOOR;
    prepotentCateg = OUTDOOR;
    OUTDOORRESP = LEFT; 
    INDOORRESP = RIGHT;
else
    counterCateg = OUTDOOR;
    prepotentCateg = INDOOR;
    OUTDOORRESP = RIGHT;
    INDOORRESP = LEFT;
end
counterCategResp = LEFT;
prepotentCategResp = RIGHT;

attnData.trialsPerRun       = trialsPerRun;
attnData.trial              = 1:trialsPerRun;
attnData.plannedTrialOnsets = nan(1,trialsPerRun);
attnData.actualTrialOnsets  = nan(1,trialsPerRun);
attnData.rts                = nan(1,trialsPerRun);
attnData.rtsEst             = nan(trialsPerRun,trialsPerRun);
attnData.rtsResid           = nan(trialsPerRun,trialsPerRun);
attnData.accs               = nan(1,trialsPerRun);
attnData.resps              = nan(1,trialsPerRun);
attnData.corrresps          = nan(1,trialsPerRun);
attnData.memPredict         = zeros(1,trialsPerRun);
attnData.forget             = zeros(1,trialsPerRun);
attnData.remember           = zeros(1,trialsPerRun);


nCounterTrials_StartBuffer = 0;
nCounterTrials_EndBuffer = 0;
targetNum_CounterTrials = ceil(nBufferTrials*proportionCounterCategory); % default is 5
attnData.categs = zeros(1,trialsPerRun);
% makes sure the number of counter trials is exactly 10% (and also that the very first trial isn't a counter trial)
while (attnData.categs(1)==counterCateg)  ||  nCounterTrials_StartBuffer ~= targetNum_CounterTrials  ||  nCounterTrials_EndBuffer ~= targetNum_CounterTrials
    
    % generate the appropriate number of prepotent and counter trials for our buffer zones (not for the middle zone since those are all just prepotent trials anyway)
    unrandomised_trial_order = cat(  2,   prepotentCateg*ones(1, (1.0-proportionCounterCategory)*nBufferTrials ) , counterCateg*ones(1, (proportionCounterCategory)*nBufferTrials )   );
    
    % the start buffer (trials 1:50) should contain 45 prepotent trials and 5 counter trials, with random distribution
    attnData.categs(1:nBufferTrials) = randsample( unrandomised_trial_order, nBufferTrials );
    % the experimental trials (middle of the experiment, i.e. trials 51:450) are all *initially* prepotent trials
    % (they can, however, be replaced by counter trials should the realtime triggering kick in)
    attnData.categs( (nBufferTrials+1):(trialsPerRun-nBufferTrials) ) = prepotentCateg;
    % the end buffer (trials 451:500) should contain 45 prepotent trials and 5 counter trials, with random distribution
    attnData.categs( (trialsPerRun-nBufferTrials+1):trialsPerRun ) = randsample( unrandomised_trial_order, nBufferTrials );
    
    % track how many counter trials are in our buffer zones so we can check if we've got the right amount
    nCounterTrials_StartBuffer = sum(attnData.categs( 1:nBufferTrials                             )==counterCateg);
    nCounterTrials_EndBuffer   = sum(attnData.categs( (trialsPerRun-nBufferTrials+1):trialsPerRun )==counterCateg);
    
end

% populate all of what the correct responses should be for each trial
attnData.corrresps((attnData.categs==prepotentCateg)) = counterCategResp;
attnData.corrresps((attnData.categs==counterCateg)) = prepotentCategResp;

categCounter = zeros(1,nSubCategs);
for iTrial = 1:trialsPerRun
    % update image counters
    categCounter(attnData.categs(iTrial)) = categCounter(attnData.categs(iTrial))+1;
    % reset counter and reshuffle images if list has been exhausted
    if (categCounter(attnData.categs(iTrial)) > numImages(attnData.categs(iTrial)))
        categCounter(attnData.categs(iTrial)) = 1; % start counter over, and reshuffle images
        imageShuffle{attnData.categs(iTrial)} = randperm(numImages(attnData.categs(iTrial)));
    end
    % get current images
    attnData.images(iTrial) = imageShuffle{attnData.categs(iTrial)}(categCounter(attnData.categs(iTrial))); 
end

%% Task Instructions

% clear screen
Screen(mainWindow,'FillRect',backColor);
Screen('Flip',mainWindow);
FlushEvents('keyDown');

% show instructions
welcomeInstruct = {'Welcome! In this task we are interested in how people categorize pictures of places';...
    'You will be shown images on the screen';...
    'Your task will be to decide whether each place is indoor or outdoor by pressing a button';...
    ['If the place is outdoor, respond with the ' KbName(OUTDOORRESP) ' key'];...
    ['If the place is indoor, respond with the ' KbName(INDOORRESP) ' key'];...
    'Use your index and middle finger of your right hand';...
    ' '};
for instruct=1:length(welcomeInstruct)
    DrawFormattedText(mainWindow, welcomeInstruct{1}, 'center', screenY*0.25);
    DrawFormattedText(mainWindow, welcomeInstruct{2}, 'center', screenY*0.30);
    DrawFormattedText(mainWindow, welcomeInstruct{3}, 'center', screenY*0.35);
    DrawFormattedText(mainWindow, welcomeInstruct{4}, 'center', screenY*0.50);
    DrawFormattedText(mainWindow, welcomeInstruct{5}, 'center', screenY*0.55);
    DrawFormattedText(mainWindow, welcomeInstruct{6}, 'center', screenY*0.60);
    DrawFormattedText(mainWindow, welcomeInstruct{7}, 'center', screenY*0.80);
    %tempBounds = Screen('TextBounds',mainWindow,welcomeInstruct{instruct});
    %Screen('drawtext',mainWindow,welcomeInstruct{instruct},centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing*(instruct-1),textColor);
    %clear tempBounds;
end
Screen('Flip',mainWindow);

% wait for experimenter to advance with space bar
FlushEvents('keyDown');
while KbCheck; end 
[~, ~, keyCode] = KbCheck(-1); % -1 checks all keyboards
while ~keyCode(KbName('space')) && ~keyCode(KbName('escape'))
    [~, ~, keyCode] = KbCheck(-1);
    if keyCode(KbName('escape')); sca; break; end
end
Screen(mainWindow,'FillRect',backColor);
Screen('Flip',mainWindow);

%% practice trials


iTrial = 1;
imageTex = Screen('MakeTexture',mainWindow,instructImages{OUTDOOR,1});
Screen('FillRect',mainWindow,backColor);
Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
Screen(mainWindow,'FillOval',fixColor,fixDotRect);
trialOnsets(iTrial) = Screen('Flip',mainWindow);
rts(iTrial) = NaN;
resps(iTrial) = NaN;
while KbCheck; end
while(resps(iTrial)~=OUTDOORRESP)
    % check for responses if none received yet
    if isnan(rts(iTrial))
        [keyIsDown, secs, keyCode] = KbCheck(-1); % -1 checks all keyboards
        if keyIsDown
            if (keyCode(OUTDOORRESP)) %|| keyCode(RIGHT))
                rts(iTrial) = secs-trialOnsets(iTrial);
                resps(iTrial) = find(keyCode,1);
                Screen('FillRect',mainWindow,backColor);
                Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
                Screen(mainWindow,'FillOval',respColor,fixDotRect);
                Screen('Flip',mainWindow);
            end
        end
    end
end

% show instructions
Screen('FillRect',mainWindow,backColor);
Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
Screen(mainWindow,'FillOval',respColor,fixDotRect);
topInstruct = {'Great!'};
tempBounds = Screen('TextBounds',mainWindow,topInstruct{1});
Screen('drawtext',mainWindow,topInstruct{1},centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing*-7,textColor);
clear tempBounds;
bottomInstruct = {'Now let''s try another place image',' '};
tempBounds = Screen('TextBounds',mainWindow,bottomInstruct{1});
Screen('drawtext',mainWindow,bottomInstruct{1},centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing*7,textColor);
clear tempBounds;
tempBounds = Screen('TextBounds',mainWindow,bottomInstruct{2});
Screen('drawtext',mainWindow,bottomInstruct{2},centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing*8,textColor);
clear tempBounds;
Screen('Flip',mainWindow);

% wait for experimenter to advance with space bar
FlushEvents('keyDown');
while KbCheck; end
[~, ~, keyCode] = KbCheck(-1); % -1 checks all keyboards
while ~keyCode(KbName('space'))
    [~, ~, keyCode] = KbCheck(-1);
end
Screen(mainWindow,'FillRect',backColor);
Screen('Flip',mainWindow);

iTrial = 1;
imageTex = Screen('MakeTexture',mainWindow,instructImages{INDOOR,1});
Screen('FillRect',mainWindow,backColor);
Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
Screen(mainWindow,'FillOval',fixColor,fixDotRect);
trialOnsets(iTrial) = Screen('Flip',mainWindow);
rts(iTrial) = NaN;
resps(iTrial) = NaN;
while KbCheck; end
while(resps(iTrial)~=INDOORRESP)
    % check for responses if none received yet
    if isnan(rts(iTrial))
        [keyIsDown, secs, keyCode] = KbCheck(-1); % -1 checks all keyboards
        if keyIsDown
            if (keyCode(INDOORRESP)) %keyCode(LEFT) ||
                rts(iTrial) = secs-trialOnsets(iTrial);
                resps(iTrial) = find(keyCode,1);
                Screen('FillRect',mainWindow,backColor);
                Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
                Screen(mainWindow,'FillOval',respColor,fixDotRect);
                Screen('Flip',mainWindow);
            end
        end
    end
end

% show instructions
Screen('FillRect',mainWindow,backColor);
Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
Screen(mainWindow,'FillOval',respColor,fixDotRect);
topInstruct = {'Great!'};
tempBounds = Screen('TextBounds',mainWindow,topInstruct{1});
Screen('drawtext',mainWindow,topInstruct{1},centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing*-7,textColor);
clear tempBounds;
bottomInstruct = {'Now let''s try a few images in a row. The images will be presented rapidly, so get ready!',' '};
tempBounds = Screen('TextBounds',mainWindow,bottomInstruct{1});
Screen('drawtext',mainWindow,bottomInstruct{1},centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing*7,textColor);
clear tempBounds;
tempBounds = Screen('TextBounds',mainWindow,bottomInstruct{2});
Screen('drawtext',mainWindow,bottomInstruct{2},centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing*8,textColor);
clear tempBounds;
Screen('Flip',mainWindow);

% wait for experimenter to advance with space bar
FlushEvents('keyDown');
while KbCheck; end
[~, ~, keyCode] = KbCheck(-1); % -1 checks all keyboards
while ~keyCode(KbName('space'))
    [~, ~, keyCode] = KbCheck(-1);
end
Screen(mainWindow,'FillRect',backColor);
Screen('Flip',mainWindow);

nInstructTrials = 10;
accs = nan(1,nInstructTrials);
instructBlockRep = 0;
while nanmean(accs)<.9 || isnan(nanmean(accs))
    
    instructBlockRep = 1;%instructBlockRep+1;
    categorder = [ones(1,5) 2*ones(1,5)];
    categorder = categorder(randperm(nInstructTrials));
    corrresps((categorder==prepotentCateg)) = prepotentCategResp;
    corrresps((categorder==counterCateg)) = counterCategResp;
    for iTrial = 1:10
        imageTex = Screen('MakeTexture',mainWindow,instructImages{categorder(iTrial),((instructBlockRep-1)*10)+1+iTrial});
        %
        Screen('FillRect',mainWindow,backColor);
        Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
        Screen(mainWindow,'FillOval',fixColor,fixDotRect);
        instructTrialOnsets(iTrial) = Screen('Flip',mainWindow);
        instructTrialTimeout(iTrial) = instructTrialOnsets(iTrial)+respWindow;
        rts(iTrial) = NaN;
        resps(iTrial) = NaN;
        
        while (GetSecs < instructTrialTimeout(iTrial)) && KbCheck; end
        while (GetSecs < instructTrialTimeout(iTrial))
                % check for responses if none received yet
                if isnan(rts(iTrial))
                    [keyIsDown, secs, keyCode] = KbCheck(-1); % -1 checks all keyboards
                    if keyIsDown
                        if (keyCode(LEFT) || keyCode(RIGHT))
                            rts(iTrial) = secs-instructTrialOnsets(iTrial);
                            resps(iTrial) = find(keyCode,1);
                            Screen('FillRect',mainWindow,backColor);
                            %if (stimOn) % leave image up if response before image duration
                            Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
                            %end
                            Screen(mainWindow,'FillOval',respColor,fixDotRect);
                            Screen('Flip',mainWindow);
                        end
                    end
                end
        end
        
        %accuracy
        if (resps(iTrial)==corrresps(iTrial)) %made correct response
            accs(iTrial) = 1;
        else
            accs(iTrial) = 0;
        end
        
    end
    
    Screen('FillRect',mainWindow,backColor);
    
    % instruction section summary
    if isnan(nanmean(accs))
        summaryText = {'We didn''t get any responses from you, so let''s try it again.',...
            ['This time, be sure to hit the ' KbName(OUTDOORRESP) ' key for outdoor scenes and the ' KbName(INDOORRESP) ' key for indoor scenes'],...
            ' ',' '};
    elseif nanmean(accs)>=.9
        summaryText = {sprintf('Great job! Your accuracy was %d%%',nanmean(accs)*100),...
            ' ',' '};
    elseif nanmean(accs)<.9
        summaryText = {sprintf('Your accuracy was only %d%%',nanmean(accs)*100),...
            'In order to continue to the main part of the experiment, you need to get 90% or better',...
            ' ','Let''s try a bit more practice before we start the task',...
            ['Remember to hit the ' KbName(OUTDOORRESP) ' key for outdoor scenes and the ' KbName(INDOORRESP) ' key for indoor scenes'],...
            'And try to respond quickly and accurately before the next image appears',...
            ' ',' '};
    end
    for instruct=1:length(summaryText)
        tempBounds = Screen('TextBounds',mainWindow,summaryText{instruct});
        Screen('drawtext',mainWindow,summaryText{instruct},centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing*(instruct-1),textColor);
        clear tempBounds;
    end
    Screen('Flip',mainWindow);
    
    % wait for experimenter to advance with space bar
    FlushEvents('keyDown');
    while KbCheck; end
    [~, ~, keyCode] = KbCheck(-1); % -1 checks all keyboards
    while ~keyCode(KbName('space'))
        [~, ~, keyCode] = KbCheck(-1);
    end
    Screen(mainWindow,'FillRect',backColor);
    Screen('Flip',mainWindow);
end

%% Experiment Instructions

% clear screen
Screen(mainWindow,'FillRect',backColor);
Screen('Flip',mainWindow);
FlushEvents('keyDown');

% show instructions
tempBounds = Screen('TextBounds',mainWindow,starttaskInstruct);
Screen('drawtext',mainWindow,starttaskInstruct,centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing,textColor);
clear tempBounds;
Screen('Flip',mainWindow);

% wait for experimenter to advance with space bar
FlushEvents('keyDown');
while KbCheck; end 
[~, ~, keyCode] = KbCheck(-1); % -1 checks all keyboards
while ~keyCode(KbName('space'))
    [~, ~, keyCode] = KbCheck(-1);
end
Screen(mainWindow,'FillRect',backColor);
Screen('Flip',mainWindow);


%% Start Experiment

Priority(MaxPriority(screenNum));
Screen(mainWindow,'FillRect',backColor);
Screen(mainWindow,'FillOval',fixColor,fixDotRect);
runStart = GetSecs;
Screen('Flip',mainWindow);
Priority(0);

% prepare for trial sequence
if subjectNum
    % trial number, time difference (between actual onset and planned onset timestamps), category ID, image ID, correct response ID, input response ID, accuracy of response, response time
    fprintf(dataFile,'trl\ttdf\tsct\tsimg\tcrsp\trsp\tacc\trt\n');
end
fprintf('trl\ttdf\tsct\tsimg\tcrsp\trsp\tacc\trt\n');
    
% timing
attnData.plannedTrialOnsets = runStart + 2 + (stimDur:stimDur:(stimDur*trialsPerRun));


%% start trial sequence

% assign new fields to attnData which are specific to realtime triggering
attnData.RT_Median = nan;
attnData.RT_Deviation = nan(1,trialsPerRun); % even though the size is 500, it's only applicable after the initial buffer
                            % because we initialise the values to nan, the initial buffer trials will default to an RT_Deviation value of nan

for iTrial=1:trialsPerRun
    
    % REALTIME TRIGGERING:
    % RT_Median is established in the first 50 trials, but only for accurate responses to prepotent trials
    % RT_Deviation of a single trial is (RT - RT_Median)
    % If RT_Deviation is high, attention state = bad = queue counter trial
    % If RT_Deviation is low, attention state = good = queue counter trial
    
    
    if (iTrial>nBufferTrials) && (iTrial<(trialsPerRun-nBufferTrials+1))
        
        %linearly regress rt (kept in for legacy data analysis)
        indFit = find(~isnan(attnData.rts(1:(iTrial-1))));
        attnData.linfit(:,iTrial) = polyfit(indFit,attnData.rts(indFit),1);
        attnData.rtsEst(iTrial,1:(iTrial-1)) = attnData.linfit(1,iTrial).*attnData.trial(1:(iTrial-1))+attnData.linfit(2,iTrial);
        attnData.rtsResid(iTrial,1:(iTrial-1)) = attnData.rts(1:(iTrial-1))-attnData.rtsEst(iTrial,1:(iTrial-1));
        attnData.meanRTresid(iTrial) = nanmean(attnData.rtsResid(iTrial,1:(iTrial-1)));
        attnData.stdRTresid(iTrial) = nanstd(attnData.rtsResid(iTrial,1:(iTrial-1)));
        
        % determine if we meet the criteria to trigger a counter trial
        % conditions for triggering one:
        % - the 3 trials before this one were all frequent
        % - not all responses/RTs for the past 3 trials were nan (i.e. lack of response) ***probably contentious (maybe none of them should be nan?)
        %   - if we are triggering a goodAttn or badAttn trial, that their respective alotted number of trials is not exceeded
        %   - some custom-defined threshold is broken (currently pending)
        if ~any(attnData.categs((iTrial-3):(iTrial-1))==counterCateg) && ~all(isnan(attnData.rts((iTrial-3):(iTrial-1))))
            if (sum(attnData.forget(1:iTrial))<nCounterTrials_BadAttn) && nanmean(attnData.rtsResid(iTrial,(iTrial-3):(iTrial-1)))<(attnData.meanRTresid(iTrial)-attnData.stdRTresid(iTrial))
                %forgettrial
                attnData.categs(iTrial) = counterCateg;
                attnData.memPredict(iTrial) = 1;
                attnData.forget(iTrial) = 1;
                attnData.remember(iTrial) = 0;
                attnData.corrresps(iTrial) = prepotentCategResp;
                categCounter(attnData.categs(iTrial)) = categCounter(attnData.categs(iTrial))+1;
                attnData.images(iTrial) = imageShuffle{attnData.categs(iTrial)}(categCounter(attnData.categs(iTrial)));
            elseif sum(attnData.remember(1:iTrial))<nCounterTrials_GoodAttn && nanmean(attnData.rtsResid(iTrial,(iTrial-3):(iTrial-1)))>(attnData.meanRTresid(iTrial)+attnData.stdRTresid(iTrial))
                %remember trial
                attnData.categs(iTrial) = counterCateg;
                attnData.memPredict(iTrial) = 1;
                attnData.forget(iTrial) = 0;
                attnData.remember(iTrial) = 1;
                attnData.corrresps(iTrial) = prepotentCategResp;
                categCounter(attnData.categs(iTrial)) = categCounter(attnData.categs(iTrial))+1;
                attnData.images(iTrial) = imageShuffle{attnData.categs(iTrial)}(categCounter(attnData.categs(iTrial)));
            end 
        end
    end
    
    % make textures
    attnData.imageName{iTrial} = imageNames{attnData.categs(iTrial),attnData.images(iTrial)};
    imageTex = Screen('MakeTexture',mainWindow,images{attnData.categs(iTrial),attnData.images(iTrial)});
    Screen('PreloadTextures',mainWindow,imageTex);
    
    % wait for trigger and show image
    FlushEvents('keyDown');
    Priority(MaxPriority(screenNum));
    Screen('FillRect',mainWindow,backColor);
    Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
    Screen(mainWindow,'FillOval',fixColor,fixDotRect);
    tRespTimeout = attnData.plannedTrialOnsets(iTrial)+respWindow; %response timeout
    
    attnData.actualTrialOnsets(iTrial) = Screen('Flip',mainWindow,attnData.plannedTrialOnsets(iTrial));
    while (GetSecs < tRespTimeout) && KbCheck; end  
    while(GetSecs < tRespTimeout)
        
        % check for responses if none received yet
        if isnan(attnData.rts(iTrial))
            [keyIsDown, secs, keyCode] = KbCheck(-1); % -1 checks all keyboards
            if keyIsDown
                if (keyCode(LEFT) || keyCode(RIGHT))
                    attnData.rts(iTrial) = secs-attnData.actualTrialOnsets(iTrial); 
                    attnData.resps(iTrial) = find(keyCode,1); 
                    Screen('FillRect',mainWindow,backColor);
                    %if (stimOn) % leave image up if response before image duration
                    Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
                    %end
                    Screen(mainWindow,'FillOval',respColor,fixDotRect);
                    Screen('Flip',mainWindow);
                end
            end
        end
    end
    
    %accuracy
    if (attnData.resps(iTrial)==attnData.corrresps(iTrial)) %made correct response
        attnData.accs(iTrial) = 1;
    else
        attnData.accs(iTrial) = 0;
    end
    
    
    % print trial results
    if subjectNum
        fprintf(dataFile,'%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%.3f\n',iTrial,attnData.actualTrialOnsets(iTrial)-attnData.plannedTrialOnsets(iTrial),attnData.categs(iTrial),attnData.images(iTrial),attnData.corrresps(iTrial),attnData.resps(iTrial),attnData.accs(iTrial),attnData.rts(iTrial));
    end
    fprintf('%d\t%.3f\t%d\t%d\t%d\t%d\t%d\t%.3f\n',iTrial,attnData.actualTrialOnsets(iTrial)-attnData.plannedTrialOnsets(iTrial),attnData.categs(iTrial),attnData.images(iTrial),attnData.corrresps(iTrial),attnData.resps(iTrial),attnData.accs(iTrial),attnData.rts(iTrial));
    
    % if we have just completed all of the initial buffer trials, then calculate RT_Median
    if iTrial==nBufferTrials
        valid_RTs = attnData.rts( attnData.categs==prepotentCateg & attnData.accs==1 );
        attnData.RT_Median = nanmedian(valid_RTs);
    end
    % if we are past the initial buffer, RT_Median is defined and we can now therefore calculate RT_Deviation for all subsequent trials
    if iTrial>nBufferTrials
        attnData.RT_Deviation(iTrial) = attnData.rts(iTrial) - attnData.RT_Median;
    end
    
end % trial loop

Screen('FillRect',mainWindow,backColor);
Screen(mainWindow,'FillOval',fixColor,fixDotRect);
endOfRun=Screen('Flip',mainWindow);

endText={'Great job! You are now done with part 1 of the experiment',' ',' '};
for instruct=1:length(endText)
    tempBounds = Screen('TextBounds',mainWindow,endText{instruct});
    Screen('drawtext',mainWindow,endText{instruct},centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing*(instruct-1),textColor);
    clear tempBounds;
end
Screen('Flip',mainWindow);


Screen('FillRect',mainWindow,backColor);
Screen(mainWindow,'FillOval',fixColor,fixDotRect);
Screen('Flip',mainWindow,endOfRun+IBI);
pause(IBI);

%% save

if subjectNum
    save([dataHeader '/attndata_' datestr(now,30)],'attnData','runStart');
end

% clean up and go home
sca;
ListenChar(1);
fclose('all');
end
