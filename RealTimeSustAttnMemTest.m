function [memData] = RealTimeSustAttnMemTest(varargin)
%
% Recognition memory test for scenes 
%
%
% REQUIRED INPUTS:
% - subjectNum:  participant number [any integer]
% - debug:       whether debugging [1 if debugging/0 if not]
%
% OUTPUTS
% - memData:   structure of experiment design and behavioral responses
%
% Original by: Megan deBettencourt
% Modified by: Steven Cao
% Version: 1.1
% Last created: April 2014
% Last modified: August 2019

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
% notice that we're only looking for the highest *already-existing* subject ID, since trial generation depends on the data generated by the previous (attention) task!
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
    subjectNum = max(listSubjectEntries);
end

% retrieve data folder reference; a new folder should not have to be created because the memory test is dependent on data generated from a previous task
dataHeader = ['rtdata/' num2str(subjectNum)];

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
% nTrials cannot be statically preallocated

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

ITI = .5;           % secs
IBI = 2;            % secs

nSubCategs=2;
OUTDOOR=1;
INDOOR=2;


%% Response Mapping and Counterbalancing

% skyra: use current design button box (keys 1,2,3,4)
ONE = KbName('1!');
TWO = KbName('2@');
THREE = KbName('3#');
FOUR = KbName('4$');

DEVICE = -1;

%starttaskinstruct = 'Please hit spacebar to start task';

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
        fprintf('The screen dimensions may be incorrect. For screenNum = %d,screenX = %d (not 1152) and screenY = %d (not 864)',screenNum, screenX, screenY);
    end
end

%create main window
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
    dirList{categ} = dirList{categ}(3:end); %skip . & ..
    if (~isempty(dirList{categ}))
        if (strcmp(dirList{categ}(1).name,'.DS_Store')==1)
            dirList{categ} = dirList{categ}(2:end);  
        end
        
        if (strcmp(dirList{categ}(end).name,'Thumbs.db')==1)
            dirList{categ} = dirList{categ}(1:(end-1));  
        end
        
        numImages(categ) = length(dirList{categ}); %#ok<AGROW>
        
        %numImages(categ) = 200;
        
        if (numImages(categ)>0)
            
            % get images
            for img=1:numImages(categ)
                
                % update progress bar
                Screen('FrameRect',mainWindow,0,progRect,10);
                Screen('FillRect',mainWindow,0,progRect);
                Screen('FillRect',mainWindow,[255 0 0],progRect-[0 0 round((1-img/numImages(categ))*progWidth) 0]);
                Screen('Flip',mainWindow);
                
                % read images
                imageNames{categ,img} = dirList{categ}(img).name;
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

% the confidence scale shown at the bottom during the experiment is actually a static image (so the four possible options you see are also separate static images)
cd images
for level = 1:5
    confRate{level} = imread(['conf/conf_rating_' num2str(level) '.jpg']); %#ok<AGROW>
    conf_textures{level} = Screen('MakeTexture',mainWindow,confRate{level}); %#ok<AGROW>
    
    if level == 1
        confSize = size(confRate{1});
        confRect = [0 0 confSize(2) confSize(1)];
        
        distConfRect = CenterRect(confRect,[0 0 screenX screenY]) + [0 200 0 200];%lower part
    end%if level
end%for level
cd ..

%%

% load attnData file
% note to self, directories should be handled in such a way that they are agnostic to the current directory
% right now this has to be run from the folder which contains the scripts
% due to the way macs work, they like to add a newline character to the end
% of the filename (yielded via ls()) for no reason, so we have to truncate that
fileDirectory = ls([dataHeader '/attndata' '_*']);
if fileDirectory(end) == newline
    fileDirectory = fileDirectory(1:end-1);
end
load( fileDirectory , 'attnData' );

%%

% to keep consistent with indoor/outdoor-frequent/rare designations from the previous (sustained attention) task
% for reference: (outdoor=1, indoor=2), (go=prepotent=frequent, nogo=counter=rare)
if mod(subjectNum,2)==1
    prepotentCateg=1;
    counterCateg=2;
    sustAttnPrepotentLocs = find(attnData.categs==prepotentCateg);
    sustAttnPrepotentImages = attnData.images(sustAttnPrepotentLocs);
    sustAttnCounterLocs = find(attnData.categs==counterCateg);
    sustAttnCounterImages = attnData.images(sustAttnCounterLocs);
else
    prepotentCateg=2;
    counterCateg=1;
    sustAttnPrepotentLocs = find(attnData.categs==prepotentCateg);
    sustAttnPrepotentImages = attnData.images(sustAttnPrepotentLocs);
    sustAttnCounterLocs = find(attnData.categs==counterCateg);
    sustAttnCounterImages = attnData.images(sustAttnCounterLocs);
end

nMemCategs = 4;
prepotentFamiliarTrials = 1; % prepotent means "frequent" category; familiar means it was in the previous task
counterFamiliarTrials = 2; % counter means "rare" category; familiar means it was in the previous task
prepotentForeignTrials = 3; % prepotent means "frequent" category; foreign means it was NOT in the previous task
counterForeignTrials = 4; % counter means "rare" category; foreign means it was NOT in the previous task

% not sure what this is
% by the way, the 'Shuffle()' function is capitalised because it is not a native MATLAB function
for i = 1:ceil(attnData.trialsPerRun/100)
    prepotentLocsChunk{i} = (find((sustAttnPrepotentLocs>((i-1)*100))&(sustAttnPrepotentLocs<((i-1)*100+101))));
    prepotentLocsShuffle{i} = Shuffle(prepotentLocsChunk{i});
    prepotentLocsSelect(i,:) = prepotentLocsShuffle{i}(1:10);
end

% memImageShuffle is a 4-slot cell array, with each slot containing images (or image IDs?) pertaining to one of the four categories
memImageShuffle{prepotentFamiliarTrials} = sustAttnPrepotentImages(Shuffle(reshape(prepotentLocsSelect,1,[])));
memImageShuffle{prepotentForeignTrials} = imageShuffle{prepotentCateg}(find(~ismember(imageShuffle{prepotentCateg},sustAttnPrepotentImages)));
memImageShuffle{counterFamiliarTrials} = sustAttnCounterImages(randperm(numel(sustAttnCounterImages)));
memImageShuffle{counterForeignTrials} = imageShuffle{counterCateg}(find(~ismember(imageShuffle{counterCateg},sustAttnCounterImages)));

% create a bunch of trials that will either pertain to the prepotent category or to the counter category, and randomise their order
memData.attCategs = Shuffle( [ prepotentCateg*ones( 1,100 ) counterCateg*ones( 1,2*numel(sustAttnCounterLocs) ) ] );
nTrials = numel(memData.attCategs);
memData.trial = 1:nTrials;

% the number of trials within each category (prepotent or counter) is matched with the number of such trials in the attention task
memData.memCategs(memData.attCategs==prepotentCateg) = Shuffle([prepotentFamiliarTrials*ones(1,50) prepotentForeignTrials*ones(1,50)]);
memData.memCategs(memData.attCategs==counterCateg)   = Shuffle([counterFamiliarTrials*ones(1,numel(sustAttnCounterLocs)) counterForeignTrials*ones(1,numel(sustAttnCounterLocs))]);

memCategCounter = zeros(1,nMemCategs);
for iTrial = 1:nTrials
    % update image counters
    memCategCounter(memData.memCategs(iTrial)) = memCategCounter(memData.memCategs(iTrial))+1;
    
    % get current images
    memData.images(iTrial) = memImageShuffle{memData.memCategs(iTrial)}(memCategCounter(memData.memCategs(iTrial))); 
%    memData.imageName{iTrial} = imageNames{memData.memCategs(iTrial),memData.images(iTrial)};
        
    if sum(memData.images(iTrial)==attnData.images)==0
        memData.attOrder(iTrial) = NaN;
    elseif (sum(memData.images(iTrial)==attnData.images)>=1) && (memData.memCategs(iTrial) <=2)
        tempInds = find(attnData.images==memData.images(iTrial));
        tempCategs = memData.attCategs(iTrial);
        memData.attOrder(iTrial) = tempInds(memData.attCategs(iTrial)==attnData.categs(tempInds));
    else
        memData.attOrder(iTrial) = NaN;
    end
end
memData.rts = nan(1,nTrials);
memData.resps = nan(1,nTrials);
memData.rating = nan(1,nTrials);
memData.accs = nan(1,nTrials);


%% Task Instructions

% clear screen
Screen(mainWindow,'FillRect',backColor);
Screen('Flip',mainWindow);
FlushEvents('keyDown');

% truncate some otherwise awkward-looking keycode names
truncated_hotkey_one   = KbName(ONE);   truncated_hotkey_one   = truncated_hotkey_one(1);
truncated_hotkey_two   = KbName(TWO);   truncated_hotkey_two   = truncated_hotkey_two(1);
truncated_hotkey_three = KbName(THREE); truncated_hotkey_three = truncated_hotkey_three(1);
truncated_hotkey_four  = KbName(FOUR);  truncated_hotkey_four  = truncated_hotkey_four(1);
% show instructions
welcomeInstruct = {'In part 2, we will be testing your memory for all of the photos that you just saw',...
    'This is the most critical part of the experiment, so please try your absolute best',...
    'You will be shown images on the screen. Your task will be to report whether you remember seeing the image',...
    ['If you remember it very strongly, press the ' truncated_hotkey_four ' key'],...
    ['If you remember it only a little bit, press the ' truncated_hotkey_three ' key'],...
    ['If you don''t think you remember seeing it, press the ' truncated_hotkey_two ' key'],...
    ['If you are sure you didn''t see it in part 1, press the ' truncated_hotkey_one ' key'],...
    'Press all keys with your index finger of your right hand',...
    'And please try to use all the keys for your responses',...
    'You will have as much time as needed to respond to each image.',...
    ' '};
% for instruct=1:length(welcomeInstruct)
%     tempBounds = Screen('TextBounds',mainWindow,welcomeInstruct{instruct});
%     Screen('drawtext',mainWindow,welcomeInstruct{instruct},centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing*(instruct-4),textColor);
%     clear tempBounds;
% end
DrawFormattedText(mainWindow, welcomeInstruct{1} , 'center', screenY*0.20);
DrawFormattedText(mainWindow, welcomeInstruct{2} , 'center', screenY*0.25);
DrawFormattedText(mainWindow, welcomeInstruct{3} , 'center', screenY*0.35);
DrawFormattedText(mainWindow, welcomeInstruct{4} , 'center', screenY*0.39);
DrawFormattedText(mainWindow, welcomeInstruct{5} , 'center', screenY*0.43);
DrawFormattedText(mainWindow, welcomeInstruct{6} , 'center', screenY*0.47);
DrawFormattedText(mainWindow, welcomeInstruct{7} , 'center', screenY*0.51);
DrawFormattedText(mainWindow, welcomeInstruct{8} , 'center', screenY*0.60);
DrawFormattedText(mainWindow, welcomeInstruct{9} , 'center', screenY*0.65);
DrawFormattedText(mainWindow, welcomeInstruct{10}, 'center', screenY*0.70);
DrawFormattedText(mainWindow, welcomeInstruct{11}, 'center', screenY*0.85);
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


% %% practice trials
% 
% iTrial = 1;
% imageTex = Screen('MakeTexture',mainWindow,instructImages{OUTDOOR,1});
% Screen('FillRect',mainWindow,backColor);
% Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
% Screen(mainWindow,'FillOval',fixColor,fixDotRect);
% trialOnsets(iTrial) = Screen('Flip',mainWindow);
% rts(iTrial) = NaN;
% resps(iTrial) = NaN;
% while KbCheck; end
% while isnan(resps(iTrial))
%     % check for responses if none received yet
%     if isnan(attnData.rts(iTrial))
%         [keyIsDown, secs, keyCode] = KbCheck(DEVICE); % -1 checks all keyboards
%         if keyIsDown
%             if (keyCode(LEFT)) %|| keyCode(RIGHT))
%                 rts(iTrial) = secs-trialOnsets(iTrial);
%                 resps(iTrial) = find(keyCode,1);
%                 Screen('FillRect',mainWindow,backColor);
%                 Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
%                 Screen(mainWindow,'FillOval',respColor,fixDotRect);
%                 Screen('Flip',mainWindow);
%             end
%         end
%     end
% end


%% Output Files Setup

if subjectNum
    % open and set-up output file
    dataFile = fopen([dataHeader '/behavior.txt'],'a');
    fprintf(dataFile,'\n*********************************************\n');
    fprintf(dataFile,'* Sustained Attention and Recognition Memory Experiment - Memory Test v1.1\n');
    fprintf(dataFile,['* Date/Time: ' datestr(now,0) '\n']);
    fprintf(dataFile,['* Seed: ' num2str(seed) '\n']);
    fprintf(dataFile,['* Subject Number: ' num2str(subjectNum) '\n']);
    fprintf(dataFile,['* debug: ' num2str(debug) '\n']);
    fprintf(dataFile,'*********************************************\n\n');
end

% print header to command window
fprintf('\n*********************************************\n');
fprintf('* Sustained Attention and Recognition Memory Experiment - Memory Test v1.1\n');
fprintf(['* Date/Time: ' datestr(now,0) '\n']);
fprintf(['* Seed: ' num2str(seed) '\n']);
fprintf(['* Subject Number: ' num2str(subjectNum) '\n']);
fprintf(['* debug: ' num2str(debug) '\n']);
fprintf('*********************************************\n\n');


%%

memRunStart = GetSecs;

if subjectNum
    % trial number, attention category ID (prepotent/counter = 1/2), memory category ID (see line 330), image ID, input response ID (keycode input), accuracy of response, response time
    fprintf(dataFile,'trl\tact\tmct\tsimg\trsp\tacc\trt\n');
end
fprintf(dataFile,'trl\tact\tmct\tsimg\trsp\tacc\trt\n');



%%

for iTrial = 1:nTrials
    % make textures
    imageTex = Screen('MakeTexture',mainWindow,images{memData.attCategs(iTrial),memData.images(iTrial)});
    Screen('PreloadTextures',mainWindow,imageTex);
    
    % wait for trigger and show image
    FlushEvents('keyDown');
    Priority(MaxPriority(screenNum));
    Screen('FillRect',mainWindow,backColor);
    Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
    Screen('DrawTexture',mainWindow,conf_textures{1},[],distConfRect);
    Screen(mainWindow,'FillOval',fixColor,fixDotRect);
    
    if iTrial==1
        memData.actualTrialOnsets(iTrial) = Screen('Flip',mainWindow,memRunStart+IBI);
    else
        memData.actualTrialOnsets(iTrial) = Screen('Flip',mainWindow,memData.actualTrialOffsets(iTrial-1)+ITI);
    end
    
    while KbCheck; end  
    while 1
        % check for responses if none received yet
        if isnan(memData.rts(iTrial))
            [keyIsDown, secs, keyCode] = KbCheck(DEVICE); % -1 checks all keyboards
            if keyIsDown
                if (keyCode(ONE)) || (keyCode(TWO)) || (keyCode(THREE)) || (keyCode(FOUR))
                    memData.rts(iTrial) = secs-memData.actualTrialOnsets(iTrial);
                    memData.resps(iTrial) = find(keyCode,1);
                    switch memData.resps(iTrial)
                        case ONE
                            memData.rating(iTrial) = 1;
                        case TWO
                            memData.rating(iTrial) = 2;
                        case THREE
                            memData.rating(iTrial) = 3;
                        case FOUR
                            memData.rating(iTrial) = 4;
                    end
                    
                    
                    if ((memData.memCategs(iTrial)==counterFamiliarTrials) || (memData.memCategs(iTrial)==prepotentFamiliarTrials)) && (memData.rating(iTrial)==3 || memData.rating(iTrial)==4)
                        memData.accs(iTrial) = 1;
                    elseif ((memData.memCategs(iTrial)== counterForeignTrials) || (memData.memCategs(iTrial)==prepotentForeignTrials)) && (memData.rating(iTrial)==1 || memData.rating(iTrial)==2)
                        memData.accs(iTrial) = 1;
                    else
                        memData.accs(iTrial) = 0;
                    end
                    break;
                end
            end
        end
    end
    
    while KbCheck; end
    
    Screen('FillRect',mainWindow,backColor);
    Screen('DrawTexture',mainWindow,imageTex,imageRect,centerRect);
    %Screen('DrawTexture',mainWindow,conf_textures{1},[],distConfRect);
    if ~isnan(memData.rating(iTrial))
        Screen('DrawTexture',mainWindow,conf_textures{memData.rating(iTrial)+1},[],distConfRect);
    else
        Screen('DrawTexture',mainWindow,conf_textures{1},[],distConfRect);
    end
    Screen(mainWindow,'FillOval',respColor,fixDotRect);
    Screen('Flip',mainWindow);
    WaitSecs(ITI); 
    
    Screen('FillRect',mainWindow,backColor);
    Screen(mainWindow,'FillOval',fixColor,fixDotRect);
    memData.actualTrialOffsets(iTrial) = Screen('Flip',mainWindow);
    
    % print trial results
    if subjectNum
        fprintf(dataFile,'%d\t%d\t%d\t%d\t%d\t%d\t%.3f\n',iTrial,memData.attCategs(iTrial),memData.memCategs(iTrial),memData.images(iTrial),memData.resps(iTrial),memData.accs(iTrial),memData.rts(iTrial));
    end
    fprintf('%d\t%d\t%d\t%d\t%d\t%d\t%.3f\n',iTrial,memData.attCategs(iTrial),memData.memCategs(iTrial),memData.images(iTrial),memData.resps(iTrial),memData.accs(iTrial),memData.rts(iTrial));
    
    
end


% show instructions
endtaskinstruct = 'Great Job! You are now done with the experiment.';
tempBounds = Screen('TextBounds',mainWindow,endtaskinstruct);
Screen('drawtext',mainWindow,endtaskinstruct,centerX-tempBounds(3)/2,centerY-tempBounds(4)/5+textSpacing,textColor);
clear tempBounds;
Screen('Flip',mainWindow);

WaitSecs(4);

%% save

if subjectNum
    save([dataHeader '/memdata_' datestr(now,30)],'memData','memRunStart');
end

% clean up and go home
sca;
ListenChar(1);
fclose('all');
end

