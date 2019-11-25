function JJF_17_6_v2
% Experiment: attend bullseye vs. attend fixation with isoeccentric band. 
%
% minor fixed made in v2
% now call wait_for_spacebar_v2 which will launch eye tracker calibration
% if "c" is pressed
%
% modified from JJF_17_2.the main tweaks were..
% 1. moved to an isoeccentric band of stimulus positions. 
% 2. as a result, there are not separate mapping and test runs. Instead, we
%    just have regular runs, and staircase runs. 
%
% Joshua Foster
% joshua.james.foster@gmail.com
% University of Chicago    
%--------------------------------------------------------------------------
% A few preliminary settings
%--------------------------------------------------------------------------
close all;
clc;
warning('off','MATLAB:dispatcher:InexactMatch'); % turn off case-mismatch manager (it's annoying)
KbName('UnifyKeyNames');

% seed the random generator
rng default % sets the generator to twister and sets the seed to 0 (as in matlab restarted)
rng shuffle % generates a new seed based on the clock
rngSettings = rng;

% save current directory
root = pwd;
% add folder with called functions to path
CalledFunctionsDir = [root,'\CalledFunctions\'];
addpath(CalledFunctionsDir);

% load gamma correction
load('EEG_A_GammaCorrection.mat');
p.calib = calib;
gammaTable = calib.gammaTable2;
clear p

% Turn of keyboard echoing
ListenChar(2);

%--------------------------------------------------------------------------
% General Experimental Parameters
%--------------------------------------------------------------------------
p.screenNum = 1;    % Grab a  for the display ID (should return 0 for single-monitor setups)
p.windowed = 0;     % 0 for full screen
p.port = 1;         % 0 = off, 1 = on
p.stimtrak = 1;     % 0 = off, 1 = on
p.eyeTrack = 1;     % 0 = off, 1 = on
p.eyeMode = 0;      % 1 = remote, 0 = chin rest
p.debug = 0;
p.isRefresh = 0;

% specific size of stimuli
p.vDist = 76; % viewing distance (cm)
p.px = 0.0275; % pixel size (cm)
p.fixSize = deg2pix(0.1,p.vDist,p.px);      % diameter of fixation point
p.stimtrakRect = [1870 4 1910 44];          % location of the stimtrak dot(in pixels) - top right corner
p.answerDist = deg2pix(0.8,p.vDist,p.px);   % position of reponse above fixation
p.questionDist = deg2pix(1.6,p.vDist,p.px);
p.cueSize = deg2pix(0.25,p.vDist,p.px);     % size of cue
p.bullseyeSize = deg2pix(1.6,p.vDist,p.px); % size of bullseye
p.bullseyeCycles = 13;                      % number of cycles in bullseye
p.stimEccentricity = deg2pix(4,p.vDist,p.px); % eccentricity of bullseye

p.nTrials = 40;
p.angPositions = [0:45:315];
p.nPositions = length(p.angPositions);

% Trial structure
p.nStim = 4; % number of stim per trial  
p.nTargEvents = 0:1; % # of target events/trial

% Timing
p.refreshRate = 120;
p.refreshCycle = 1/p.refreshRate;
p.startBufferDur = 0.4;
p.cueDur = 0.3;
p.isi = [0.5:p.refreshCycle:0.8]; % increased longest ISI from 0.7 in JJF_17_2 to 0.8 here
p.stimDur = 0.1;
p.endBufferDur = p.isi(1); % the shortest possible ISI duration       

% Color information
p.txtCol = 220;
p.white = 255;
p.foreCol = 128;
p.fix_ntLum = 200;
p.cueLum = 170; % brightness of cue (dim to minimize after-image)
p.bullseye_ntContrast = 0.9; % contrast of standard stim

% Responses
p.acceptKeys = [KbName('SPACE'),KbName('0'),KbName('1')];

%--------------------------------------------------------------------------
% Collect information
%--------------------------------------------------------------------------

% Subject number
prompt = {'Subject Number'};
defAns = {''}; 
box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);           
if length(box) == length(defAns)                                           
    p.subNum = num2str(box{1}); 
else
    return;
end

% load in subject info based on subject number
fname = [root,'\Subject Data\Sub',num2str(p.subNum),'\',num2str(p.subNum),'_SubInfo.mat'];
exists = exist(fname,'file'); 
if exists == 2
    load(fname)
    p.sub_info = sub_info; % save to p struct
    p.taskOrder = sub_info.taskOrder; % get task order from sub info
else
    error('subject info file does not exist');
end

% Run type
prompt = {'Run type? Staircase = 1, Regular = 2'};
defAns = {''};
box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);
if length(box) == length(defAns)
    p.runType = str2num(box{1});
else
    return;
end

% Run number
prompt = {'Run number?'};
defAns = {''};
box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);
if length(box) == length(defAns)
    p.runNum = str2num(box{1});
else
    return;
end

% Set bullseye contrast difference
prompt = {'Bullseye contrast difference?'};
defAns = {''};
box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);
if length(box) == length(defAns)
    p.bullseye_contrastDiff = str2num(box{1});
else
    return;
end

% Set fixation luminance difference
prompt = {'Luminance difference'};
defAns = {''};
box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);
if length(box) == length(defAns)
    p.fix_lumDiff = str2num(box{1});
else
    return;
end
     
%--------------------------------------------------------------------------
% Setup block sequence
%--------------------------------------------------------------------------
% staircase run
if p.runType == 1
    if p.taskOrder == 1
        p.attnCond = [1 2 1 2 1 2];
    else
        p.attnCond = [2 1 2 1 2 1];
    end
    staircase.attnCond = repmat(p.attnCond',1,p.nTrials); 
    staircase.fix = nan(length(p.attnCond),p.nTrials);
    staircase.fix_change = zeros(length(p.attnCond),p.nTrials);
    staircase.bullseye = nan(length(p.attnCond),p.nTrials);
    staircase.bullseye_change = zeros(length(p.attnCond),p.nTrials);
end

% regular run
if p.runType == 2
    if p.taskOrder == 1
        p.attnCond = [1 2 1 2];
    else
        p.attnCond = [2 1 2 1];
    end
end
   
% number of blocks
p.nBlocks = length(p.attnCond);

%--------------------------------------------------------------------------
% Build the stimulus display
%--------------------------------------------------------------------------
AssertOpenGL;                                                              % Make sure we're on an OpenGL-compatible machine

if p.windowed                                                              % if we're debugging, open up a small window
    [w,p.sRect] = Screen('OpenWindow',p.screenNum,p.foreCol,[120 120 1024,768]);
else
    [w,p.sRect] = Screen('OpenWindow',p.screenNum,p.foreCol);              % otherwise, go full-screen
    HideCursor;
    Priority(MaxPriority(w));                                              % set priority to max to discourage interruptions
end

% apply gamma correction
old = Screen('LoadNormalizedGammaTable',w,calib.gammaTable2*[1 1 1]);

% Turn on alpha blending. this makes drawing the stims much easier.
% Type Screen('BlendingFunction?') for more info.
Screen('BlendFunction',w,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

% Compute and store the center of the screen
p.xCenter = (p.sRect(3)-p.sRect(1))/2;
p.yCenter = (p.sRect(4)-p.sRect(2))/2;

% Foreground rectangle
p.foreRect = p.sRect; % set to the full screen

% Fixation rectangle
p.fixRect = CenterRect([0 0 p.fixSize p.fixSize],p.foreRect);

% Bullseye rectangle
p.bullseyeRect = CenterRect([0 0 p.bullseyeSize p.bullseyeSize],p.foreRect);

% cue rectangle
p.cueRect = CenterRect([0 0 p.cueSize p.cueSize],p.foreRect);

% make a dummy call to GetSecs to load the .dll before we need it
dummy = GetSecs; clear dummy;

%----------------------------------------------------
% Setup portcodes
%----------------------------------------------------
if p.port
    config_io;
    portCode = hex2dec('D050');
end

%----------------------------------------------------
% Initiate the Eyetracker.
%----------------------------------------------------
edfFile = []; % create so that wait_for_spacebar works, even when no eye tracking

if p.eyeTrack
    
    % initialization of the connection with the Eyelink Gazetracker
    if EyelinkInit()~= 1;return;end

    % load in default settings
    EyeLinkDefaults=EyelinkInitDefaults(w); 
    
    % update defaults
    EyeLinkDefaults.backgroundcolour = p.foreCol(1);    % match calibration background color to p.foreCol
    EyeLinkDefaults.calibrationtargetcolour = p.txtCol; % match calibration target color to p.txtCol
    EyeLinkDefaults.msgfontcolour = p.txtCol(1);        % change font color
    EyeLinkDefaults.imgtitlecolour = p.txtCol(1);       % change font color
    EyelinkUpdateDefaults(EyeLinkDefaults);             % update the eyetrack settings
   
    % setup eye tracker
    if p.eyeMode
        % remote mode
        Eyelink('command', 'elcl_select_configuration = RTABLER'); % remote mode
        Eyelink('command', 'calibration_type = HV5'); % 5-pt calibration
        Eyelink('command', 'sample_rate = 500'); % sample at 500 Hz for remote mode
    else
        % chin rest
        Eyelink('command', 'elcl_select_configuration = MTABLER'); % chin rest
        Eyelink('command', 'calibration_type = HV9'); % 9-pt calibration
        Eyelink('command', 'sample_rate = 1000'); % sample at 1000 Hz for chin rest
    end
    
    % stamp header in EDF file
    Eyelink('command', 'add_file_preamble_text','JJF_17_6'); % 
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    [width, height]=Screen('WindowSize',w);
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);
    % make sure that we get gaze data from the Eyelink
    % set EDF file contents using the file_sample_data and
    % file-event_filter commands
    Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,HTARGET,GAZERES,STATUS,INPUT');
    % set link data thtough link_sample_data and link_event_filter
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
    Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');
    % proportion commands to adjust size of calibrated area
    Eyelink('command', 'calibration_area_proportion 0.5 0.5')
    Eyelink('command', 'validation_area_proportion 0.5 0.5')
    
    % get host tracker version
    [v,vs]=Eyelink('GetTrackerVersion');
    
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    fprintf('Running experiment on version ''%d''.\n', v );
    
    % open file to record data to (filename must be fewer than 8 characters)
    if p.runType == 1 % staircase run
        edfFile = [num2str(p.subNum),'JF6S',num2str(p.runNum),'.edf']; 
    end
    if p.runType == 2 % regular run
        edfFile = [num2str(p.subNum),'JF6R',num2str(p.runNum),'.edf']; 
    end
    
    Eyelink('Openfile', edfFile);
    
end

%---------------------------------------------------
% Begin block loop
%---------------------------------------------------
for b = 1:p.nBlocks
        
    % Setup eyetracker
    if p.eyeTrack == 1
        EyelinkDoTrackerSetup(EyeLinkDefaults);
    end
    
    %----------------------------------------------------------------
    % Setup place to save files
    %----------------------------------------------------------------
    
    % Build an output directory if it doesn't already exist
    if ~exist([root, '\Subject Data\'],'dir')
        mkdir([root, '\Subject Data\']);
    end
    
    % subfolder for individual subject's data
    if ~exist([root, '\Subject Data\Sub',num2str(p.subNum)],'dir')
        mkdir([root, '\Subject Data\Sub',num2str(p.subNum)]);
    end

    % Build an output file and check to make sure that it doesn't exist yet either
    
    % staircase run
    if p.runType == 1
        fName = [root,'\Subject Data\Sub',num2str(p.subNum),'\',num2str(p.subNum),'_JJF_17_6_S',num2str(p.runNum),'_B',num2str(b),'.mat'];
    end
    % regular run
    if p.runType == 2
        fName = [root,'\Subject Data\Sub',num2str(p.subNum),'\',num2str(p.subNum),'_JJF_17_6_R',num2str(p.runNum),'_B',num2str(b),'.mat'];
    end

    
    % check that no such data file already exists
    if exist(fName)
        Screen('CloseAll');
        msgbox('File already exists!', 'modal');
        return;
    end
    
    %--------------------------------------------------------
    % Timing vectors
    %--------------------------------------------------------
    time.startBuffer.vblstamp = nan(1,p.nTrials);
    time.startBuffer.onset = nan(1,p.nTrials);
    time.startBuffer.flipstamp = nan(1,p.nTrials);
    time.startBuffer.missed = nan(1,p.nTrials);
    
    time.cue.vblstamp = nan(1,p.nTrials);
    time.cue.onset = nan(1,p.nTrials);
    time.cue.flipstamp = nan(1,p.nTrials);
    time.cue.missed = nan(1,p.nTrials);
        
    time.isi.vblstamp = nan(p.nTrials,p.nStim);
    time.isi.onset = nan(p.nTrials,p.nStim);
    time.isi.flipstamp = nan(p.nTrials,p.nStim);
    time.isi.missed = nan(p.nTrials,p.nStim);
    
    time.stim.vblstamp = nan(p.nTrials,p.nStim);
    time.stim.onset = nan(p.nTrials,p.nStim);
    time.stim.flipstamp = nan(p.nTrials,p.nStim);
    time.stim.missed = nan(p.nTrials,p.nStim);
        
    time.endBuffer.vblstamp = nan(1,p.nTrials);
    time.endBuffer.onset = nan(1,p.nTrials);
    time.endBuffer.flipstamp = nan(1,p.nTrials);
    time.endBuffer.missed = nan(1,p.nTrials);  
       
    %--------------------------------------------------------
    % Control parameters   
    %--------------------------------------------------------   
    % get the current condition
    stim.attnCond = ones(1,p.nTrials)*p.attnCond(b); % create trial-by-trial index of attn condition
    
    % note: the contrast/luminance variables will be overwritten in the
    % staircasing blocks...
    
    % contrast of bullseye
    stim.bullseye_ntContrast = ones(1,p.nTrials)*p.bullseye_ntContrast;
    stim.bullseye_contrastDiff = ones(1,p.nTrials)*p.bullseye_contrastDiff;
    if p.runType == 1 && b > 1 % carry forward staircasing variables
        stim.bullseye_contrastDiff = ones(1,p.nTrials)*staircase.bullseye(b-1,end);
    end
    stim.bullseye_tContrast = stim.bullseye_ntContrast - stim.bullseye_contrastDiff;
    
    % luminance of fixation
    stim.fix_ntLum = ones(1,p.nTrials)*p.fix_ntLum;
    stim.fix_lumDiff = ones(1,p.nTrials)*p.fix_lumDiff;
    if p.runType == 1 && b > 1 % carry forward staircasing variables
        stim.fix_lumDiff = ones(1,p.nTrials)*staircase.fix(b-1,end);
    end
    stim.fix_tLum = stim.fix_ntLum - stim.fix_lumDiff;
    
    % response vectors
    stim.resp = nan(1,p.nTrials);
    stim.acc = nan(1,p.nTrials);
    
    % randomly select the ISI durations: 
    tmpisi = datasample(p.isi,p.nTrials*p.nStim); % randomly sample ISIs
    stim.isi = reshape(tmpisi,[p.nTrials, p.nStim]); % reshape into matrix
        
    % determine order of stimulus positions
    repPerPos = p.nTrials/p.nPositions;
    tmpStimPos = repmat(p.angPositions,1,repPerPos);
    stim.stimPos = Shuffle(tmpStimPos);
    
    % determine order of target events
    stim.nTargs = setEventsPerTrial(p.nTrials,p.nTargEvents);
    
    % determine order of non-target events
    stim.nNonTargs = setEventsPerTrial(p.nTrials,p.nTargEvents);
    
    % preallocate matrices that specify target and non-target events
    stim.tMat = zeros(p.nTrials,p.nStim);
    stim.ntMat = zeros(p.nTrials,p.nStim);   
              
    %---------------------------------------------------------------------%
    % determine details of each trial
    %---------------------------------------------------------------------%
    for t = 1:p.nTrials
               
        % set stim coordinates
        stim.tycoord(t) = p.yCenter+sind(stim.stimPos(t))*p.stimEccentricity;
        stim.txcoord(t) = p.xCenter+cosd(stim.stimPos(t))*p.stimEccentricity;
        stim.bullseyeCoords(t,:) = CenterRectOnPointd(p.bullseyeRect,stim.txcoord(t),stim.tycoord(t));
        stim.cueCoords(t,:) = CenterRectOnPointd(p.cueRect,stim.txcoord(t),stim.tycoord(t));
        
        % determine event sequence for each trial
        stim.tMat(t,:) = setEventSeq(p.nStim,stim.nTargs(t));
        stim.ntMat(t,:) = setEventSeq(p.nStim,stim.nNonTargs(t));
        
    end 
    %---------------------------------------------------------------------%
    
    % specify whether stim of fix is the target
    if p.attnCond(b) == 1
        stim.bullseyeEventMat = stim.tMat;
        stim.fixEventMat = stim.ntMat;
    elseif p.attnCond(b) == 2
        stim.bullseyeEventMat = stim.ntMat;
        stim.fixEventMat = stim.tMat;
    end

    %-------------------------------------------------------
    % Start of block screen
    %-------------------------------------------------------
    HideCursor;
    % specify the attention condition
    if stim.attnCond == 1
        text1 = 'In this block you will ATTEND THE BULLSEYE';
        tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-140];
    elseif stim.attnCond == 2
        text1 = 'In this block you will ATTEND THE FIXATION POINT';
        tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-140];
    end
    Screen('FillRect',w,p.foreCol,p.foreRect);    % background screen
    Screen('DrawText', w, text1, tCenter1(1), tCenter1(2), [255 0 0]);
    Screen('FillOval',w,p.fix_ntLum,p.fixRect) % fixation point
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    WaitSecs(3);
    
    % Wait for a spacebar press to continue (or esc to abort)
    wait_for_spacebar_v2(p,stim,w,edfFile,EyeLinkDefaults)
    
    %-------------------------------------------------------
    % Trial loop
    %-------------------------------------------------------
    
    % create this here - constant across trials
    nt = makeBullseye(w,p.foreCol,p.bullseyeCycles,stim.bullseye_ntContrast(t)); % full-contrast standard
    
    for t = 1:p.nTrials
        
        % create bullseye textures for current trial. note: done on each trial for staircase procedure
        tg = makeBullseye(w,p.foreCol,p.bullseyeCycles,stim.bullseye_tContrast(t)); % low-contrast target
        fix_tLum = stim.fix_tLum(t); % luminance of fixation target
        
        %------------------------------------------------------------------
        % Initiate trial screen
        %------------------------------------------------------------------
        HideCursor;
        % specify the attention condition
        if stim.attnCond == 1
            text1 = 'Attend bullseye';
            tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-140];
        else
            text1 = 'Attend fixation';
            tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-140];
        end
        Screen('FillRect',w,p.foreCol,p.foreRect);    % background screen
        Screen('DrawText', w, text1, tCenter1(1), tCenter1(2), p.txtCol);
        Screen('FillOval',w,p.fix_ntLum,p.fixRect) % fixation point
%         Screen('FillOval',w,p.cueLum,squeeze(stim.cueCoords(t,:))); % Draw marker that shows bullseye position
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        WaitSecs(.1);
        
        % Wait for a spacebar press to continue (or esc to abort)
        wait_for_spacebar_v2(p,stim,w,edfFile,EyeLinkDefaults)
                        
        %------------------------------------------------------------------
        % Start-of-trial buffer display
        %------------------------------------------------------------------
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fix_ntLum,p.fixRect);
        Screen('DrawingFinished',w);
        
        % flip buffer display
        [time.startBuffer.vblstamp(t) time.startBuffer.onset(t) time.startBuffer.flipstamp(t) time.startBuffer.missed(t)] = Screen('Flip',w);
        if p.port == 1
            outp(portCode,b)
        end
        
        %------------------------------------------------------------------
        % tell eye tracker that the trial is beginning.
        %------------------------------------------------------------------
        if p.eyeTrack
            % set idle mode before start of recording
            Eyelink('Command', 'set_idle_mode');
            % Put block, trial number at the bottom of operater display
            Eyelink('command', 'record_status_message ''BLOCK %d TRIAL %d''', b, t)
            % start recording
            Eyelink('StartRecording')
            Eyelink('message', 'BLOCK %d ', b);
            Eyelink('message', 'TRIAL %d ', t);
            % mark zero-plot time in EDF file
            Eyelink('Message', 'TrialStart');  % may be a lag in recording (not a problem because the baseline screen is much longer than needed.
        end
        
        %------------------------------------------------------------------
        % cue display
        %------------------------------------------------------------------
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fix_ntLum,p.fixRect);
        Screen('FillOval',w,p.cueLum,squeeze(stim.cueCoords(t,:))); % Draw marker that shows bullseye position
        Screen('DrawingFinished',w);
        
        % flip cue display
        WaitSecs('UntilTime',time.startBuffer.onset(t)+p.startBufferDur-0.5*p.refreshCycle);
        [time.cue.vblstamp(t) time.cue.onset(t) time.cue.flipstamp(t) time.cue.missed(t)] = Screen('Flip',w);
        
        if p.port == 1
           outp(portCode,49) 
        end
        if p.eyeTrack
            Eyelink('message','CueOnset');
        end
        
        
        % loop through repeated ISI-stim displays
        for s = 1:p.nStim
            
            %----------------------------------------------
            % ISI display
            %----------------------------------------------
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fix_ntLum,p.fixRect);
            Screen('DrawingFinished',w);
            
            % flip ISI display
            if s == 1
                WaitSecs('UntilTime',time.cue.onset(t)+p.cueDur-0.5*p.refreshCycle);
            else
                WaitSecs('UntilTime',time.stim.onset(t,s-1)+p.stimDur-0.5*p.refreshCycle);
            end
        
            [time.isi.vblstamp(t,s) time.isi.onset(t,s) time.isi.flipstamp(t,s) time.isi.missed(t,s)] = Screen('Flip',w);
            if p.port == 1               
                outp(portCode,50)
            end
            if p.eyeTrack
                Eyelink('message','ISIOnset');
            end
            
            %----------------------------------------------
            % Stim display
            %----------------------------------------------
            Screen('FillRect',w,p.foreCol,p.foreRect);
            
            % determine chekerboard to present
            if stim.fixEventMat(t,s) == 1
                Screen('FillOval',w,fix_tLum,p.fixRect);
            else
                Screen('FillOval',w,p.fix_ntLum,p.fixRect);
            end
            
            % determine fixation point to present
            if stim.bullseyeEventMat(t,s) == 1
                Screen('DrawTexture',w,tg,[],stim.bullseyeCoords(t,:),[],1,[],[]);
            else
                Screen('DrawTexture',w,nt,[],stim.bullseyeCoords(t,:),[],1,[],[]);
            end
            
            % stimtrack
            if p.stimtrak
                Screen('FillOval',w,p.white,p.stimtrakRect); % draw circle for stimtrak
            end
            
            Screen('DrawingFinished',w);
            WaitSecs('UntilTime',time.isi.onset(t,s)+stim.isi(t,s)-0.5*p.refreshCycle);
                        
            % flip stim display
            [time.stim.vblstamp(t,s) time.stim.onset(t,s) time.stim.flipstamp(t,s) time.stim.missed(t,s)] = Screen('Flip',w);
            if p.port == 1               
                outp(portCode,51)
            end
            if p.eyeTrack
                Eyelink('message','StimOnset');
            end
            
        end % end of stim loop
               
        %----------------------------------------------
        % End-of-trial buffer display
        %----------------------------------------------
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fix_ntLum,p.fixRect);
        Screen('DrawingFinished',w);
        
        WaitSecs('UntilTime',time.stim.onset(t,p.nStim)+p.stimDur-0.5*p.refreshCycle);

        % flip buffer display
        [time.endBuffer.vblstamp(t) time.endBuffer.onset(t) time.endBuffer.flipstamp(t) time.endBuffer.missed(t)] = Screen('Flip',w);
        
         WaitSecs('UntilTime',time.endBuffer.onset(t)+p.endBufferDur-0.5*p.refreshCycle);
              
        % send end-of-trial trigger
        if p.port == 1
            outp(portCode,100+t)
        end
        % stop recording eye data for current trial
        if p.eyeTrack == 1
            Eyelink('message', 'TrialEnd');
            Eyelink('StopRecording');
        end;
         
        %----------------------------------------------
        % Get response
        %----------------------------------------------       
        % Draw screen without fix
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fix_ntLum,p.fixRect);
        text = 'target present?';
        tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter-p.questionDist];
        Screen('DrawText', w, text, tCenter1(1), tCenter1(2), p.txtCol);
        Screen('DrawingFinished',w);
        Screen('Flip',w);
        
        resp = getResponse(p,w,t,stim);
        resp = str2num(KbName(resp));
        stim.resp(t) = resp;
        
        % calculate accuracy
        if stim.resp(t) == stim.nTargs(t)
            stim.acc(t) = 1;
        else
            stim.acc(t) = 0;
        end
              
        %----------------------------------------------
        % Display feedback
        %----------------------------------------------
        if p.runType == 1 % only display for staircase blocks
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fix_ntLum,p.fixRect);
            if stim.acc(t) == 1
                feedback = 'Correct';
            else
                feedback = 'Incorrect';
            end
            tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, feedback))/2 p.yCenter-p.answerDist];
            Screen('DrawText', w, feedback, tCenter1(1), tCenter1(2), p.txtCol);
            Screen('DrawingFinished',w);
            Screen('Flip',w);
            WaitSecs(0.5);
        end
        
        %------------------------------------------------------------------
        % Staircase procedure
        %------------------------------------------------------------------   
        if p.runType == 1 % if staircase block
            
            % get accuracy from current trial
            currTrialAcc = stim.acc(t);
            % get accuracy from last three trials (only starting at trial 3)
            recentTrialAcc = nan(1,3);
            if t > 2
                recentTrialAcc = [stim.acc(t-2:t)];
            end

            % staircase bullseye contrast
            if p.attnCond(b) == 1
                % get current value
                currVal = stim.bullseye_contrastDiff(t);
                % get recent values
                if t > 2
                    recentVals = [stim.bullseye_contrastDiff(t-2:t)];
                else
                    recentVals = nan(3);
                end
                % do staircasing
                [newVal changeDir] = runStaircase(b,t,currVal,currTrialAcc,recentTrialAcc,recentVals );
                % save variables to staircase struct
                staircase.bullseye(b,t) = newVal;
                staircase.bullseye_change(b,t) = changeDir; 
                staircase.fix(b,t) = stim.fix_lumDiff(t);
                % save new value
                if t < p.nTrials
                    stim.bullseye_contrastDiff(t+1) =  newVal;
                    stim.bullseye_tContrast(t+1) = p.bullseye_ntContrast - stim.bullseye_contrastDiff(t+1);
                end
            end
            
            % staircase fixation luminance
            if p.attnCond(b) == 2
                % get current value
                currVal = stim.fix_lumDiff(t);
                % get recent values
                if t > 2
                    recentVals = [stim.fix_lumDiff(t-2:t)];
                else
                    recentVals = nan(3);
                end
                % do staircasing
                [newVal changeDir] = runStaircase(b,t,currVal,currTrialAcc,recentTrialAcc,recentVals);
                % save variables to staircase struct
                staircase.fix(b,t) = newVal;
                staircase.fix_change(b,t) = changeDir;
                staircase.bullseye(b,t) = stim.bullseye_contrastDiff(t);
                % save new value
                if t < p.nTrials
                    stim.fix_lumDiff(t+1)=  newVal;
                    stim.fix_tLum(t+1) = p.fix_ntLum - stim.fix_lumDiff(t+1);
                end
            end
        end
       
        %------------------------------------------------------------------

        % save data file at the end of each trial
        if p.runType == 1
            save(fName,'p','stim','staircase'); % only save staircase struct for staircase blocks
        else
            save(fName,'p','stim');
        end
        
    end % end of trial loop
    
    %----------------------------------------------
    % End of block feedback
    %---------------------------------------------- 
    if p.runType ~= 1 % don't do for staircase blocks
    accuracy =  round(100*(mean(stim.acc))); % calculate block accuracy
    
    % present feedback
    text1 = ['Accuracy: ' num2str(accuracy) '%'];
    tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-120];
    text2 = 'Press the spacebar to continue.';
    tCenter2 = [p.xCenter-RectWidth(Screen('TextBounds', w, text2))/2 p.yCenter-70];
    Screen('FillRect',w,p.foreCol,p.foreRect);            % Draw the foreground window
    Screen('FillOval',w,p.fix_ntLum,p.fixRect);           % Draw the fixation point
    Screen('DrawText', w, text1, tCenter1(1), tCenter1(2), p.txtCol);
    Screen('DrawText', w, text2, tCenter2(1), tCenter2(2), p.txtCol);
    Screen('DrawingFinished', w);
    Screen('Flip', w);
    WaitSecs(.5);
    
    % Wait for a spacebar press to continue
    wait_for_spacebar(p,w)
    end
    %-----------------------------------------------

end % end of block loop

ListenChar(0);
% pack up and go home
if p.eyeTrack == 1
    grabEDF(edfFile);
end
LoadIdentityClut(w); % restore identity CLUT
Screen('CloseAll');
ShowCursor('Arrow');
