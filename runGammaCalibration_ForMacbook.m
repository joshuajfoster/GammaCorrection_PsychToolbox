    function runGammaCalibration_v2
% This is a wrapper function to run the built-in PsychToolbox monitor
% calibration. Produces a file that includes the gamma correction to be
% used.
%
% Instructions for use:
%
% 1. The BenQ monitor should be set to the standard color setting. To do
% this, menu-system-reset all, then menu-picture-set to standard.
% 2. Put some electrical tape over the monitor controls so people don't mess 
% with the settings!
% 3. Place the photometer at the correct viewing distance (e.g., where the
% chin rest is). Keep the photometer as steady as possble because
% variations in viewing distance will affect the readings.
% 4. The script will present a number of different  shades of grey (not
% 50,currently set to 33). Take a reading for each and report on the
% numpad.
%
% notes: 
% you can chose thenumber of readings to take (p.numMeaures). Steve Shevell recommends at least 16.
% you should try this a couple of times to make sure you're getting consistent results.
% The CalibirateMonitorPhotometer_JJFfix(function below) is the standard function from PTB
% with 'UnifyKeys' added at the top so it runs properly.
%
% Joshua J. Foster
% January 22, 2017
%
% v2 update (7.9.2018)
% script now displays input, we can edit if we need to (JJF made changes)
% levels presented in random order (Janna Wennberg made changes)

% ask for experimenter details
prompt = {'Experimenter Initials'};
defAns = {''}; 
box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);             % actually make GUI
if length(box) == length(defAns)                                            % simple check for enough input, otherwise bail
    p.expName = num2str(box{1}); 
else
    return;
end

% ask for room details
prompt = {'Room details  (e.g., EEG_A or BEH_C)'};
defAns = {''}; 
box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);             % actually make GUI
if length(box) == length(defAns)                                            % simple check for enough input, otherwise bail
    p.roomName = num2str(box{1}); 
else
    return;
end

% ask for date
prompt = {'Date mm.dd.yyyy'};
defAns = {''}; 
box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);             % actually make GUI
if length(box) == length(defAns)                                            % simple check for enough input, otherwise bail
    p.date = num2str(box{1}); 
else
    return;
end

Screen('Preference', 'SkipSyncTests', 1);

p.numMeasures = 17 ; % Steve Shevell recommends at least 16 levels
p.nSamps = 2; % number of measuremens for each level

s = max(Screen('Screens'));     s = 1; 
[win,p.sRect] = Screen('OpenWindow',s,0,[]);

% load in identity gamma table
LoadIdentityClut(win);

% Screen('FillRect',win,[255 255 255],p.sRect);
% text1 = ['Adjust monitor brightness until white is ~ 120 cd/m2'];
% text2 = ['then press spacebar'];
% Screen('DrawText', win, text1,[0],[0],[200 200 200]);
% Screen('DrawText', win, text2,[0],[30],[200 200 200]);
% Screen('DrawingFinished', win);
%  Screen('Flip',win);
% % Wait for a space bar press to con tinue
% while 1
%     [keyIsDown,secs,keyCode]=KbCheck;
%     if keyIsDown
%         kp = find(keyCode);
%         if kp == 32
%             break;
%         end
%     end
% end
% 
% % close screen
% Screen('CloseAll');
% 
% % ask for monitor brightness
% prompt = {'What was the brightness value on the monitor needed to achieve desired brightness?'};
% defAns = {''}; 
% box = inputdlg(prompt,  'Enter Subject Information', 1, defAns);             % actually make GUI
% if length(box) == length(defAns)                                            % simple check for enough input, otherwise bail
%     calib.monitorBrightness  = num2str(box{1}); 
% else
%     return;
% end

% Run calaibration script
[calib.gammaTable1, calib.gammaTable2, calib.displayBaseline, calib.displayRange, calib.displayGamma, calib.maxLevel calib.measurements] = CalibrateMonitorPhotometer_JoshFix(p.numMeasures,p.nSamps,[]);

% save file to be used later
fName = [p.roomName,'_GammaCorrection.mat']
save(fName,'p','calib');
saveas(gcf,[p.roomName,'_GammaCorrection_Fits'])

%%-------------------------------------------------------------------------
function [ gammaTable1, gammaTable2, displayBaseline, displayRange, displayGamma, maxLevel measurements] = CalibrateMonitorPhotometer_JoshFix(numMeasures,nSamps,screenid)
% [gammaTable1, gammaTable2, displayBaseline, displayRange. displayGamma, maxLevel] = CalibrateMonitorPhotometer([numMeasures=9][, screenid=max])
%
% A simple calibration script for analog photometers.
%
% Use CalibrateMonSpd() if you want to do more fancy calibration with
% different types of photometers or special devices like Bits+ or DataPixx,
% assuming you know how to operate CalibrateMonSpd() that is...
%
% numMeasures (default: 9) readings are taken manually, and the readings
% are fit with a gamma function and piecewise cubic splines. numMeasures -
% 1 should be a power of 2, ideally (9, 17, 33, etc.). The corresponding
% linearized gamma tables (1 -> gamma, 2 -> splines) are returned, as well
% as the display baseline, display range in cd/m^2 and display gamma. Plots
% of the two fits are created as well. Requires fit tools.
%
% If the normalized gamma table is not loaded, then the cd/m^2 value of a
% screen value can be figured out by the formula: cdm2 =
% displayRange*(screenval/maxLevel).^(1/displayGamma) + displayBaseline.
%
% Generally, you will want to load the normalized gamma tables and use them
% in Screen('LoadNormalizedGammaTable'). For example:
%
% [gammaTable1, gammaTable2] = CalibrateMonitorPhotometer;
% %Look at the outputted graphs to see which one gives a better fit
% %Then save the corresponding gamma table for later use
% gammaTable = gammaTable1;
% save MyGammaTable gammaTable
% 
% %Then when you're ready to use the gamma table:
% load MyGammaTable
% Screen('LoadNormalizedGammaTable', win, gammaTable*[1 1 1]);
%
%
% History:
% Version 1.0: Patrick Mineault (patrick.mineault@gmail.com)
% 22.10.2010 mk Switch numeric input from use of input() to use of
%               GetNumber(). Restore gamma table after measurement. Make
%               more robust.
% 19.08.2012 mk Some cleanup.
%  4.09.2012 mk Use Screen('ColorRange') to adapt number/max of intensity
%               level to given range of framebuffer.

KbName('UnifyKeyNames');

% global vals;
% global inputV;

    if (nargin < 1) || isempty(numMeasures)
        numMeasures = 9;
    end

    input(sprintf(['When black screen appears, point photometer, \n' ...
           'get reading in cd/m^2, input reading using numpad and press enter. \n' ...
           'A screen of higher luminance will be shown. Repeat %d times. ' ...
           'Press enter to start'], numMeasures));
       
    psychlasterror('reset');    
    try
        if nargin < 2 || isempty(screenid)
            % Open black window on default screen:
            screenid = max(Screen('Screens'));
              screenid  = 1;
        end
        
        % Open black window:
        win = Screen('OpenWindow', screenid,0);
        maxLevel = Screen('ColorRange', win);

        % Load identity gamma table for calibration:
        LoadIdentityClut(win);

        measurements  = nan(nSamps,numMeasures); % observed luminance
        inputV = [0:(maxLevel+1)/(numMeasures - 1):(maxLevel+1)];
        inputV(end) = maxLevel;
 
        % repeat each level nSamp m any times
        for s = 1:nSamps
            shuffledInputV = Shuffle(inputV);
            for i = 1:numMeasures
                level = shuffledInputV(i);
                  
                %----------------------------------------------------------
                % get a response
                %----------------------------------------------------------
                while 1
                
                % get a response (type in number, can backspace if needed)
                Screen('FillRect',win,level);
                text = ['RGB = ',num2str(level),',  Luminance = '];
                resp = GetEchoNumber(win,text ,10,10,[200 0 0],level); 
                
                % check the response
                Screen('FillRect',win,level);
                text = ['Is ',num2str(resp),' correct?  (1 = yes, 0 = no)'];
                corrCheck = GetEchoNumber(win,text ,10,10,[200 0 0],level);
                
                % exit loop if response is correct
                if corrCheck == 1 
                    break 
                end
                
                end
                %----------------------------------------------------------

                % save response to the measurements matrix
                k = inputV == level;
                measurements(s,k) = resp;
            end
        end
        
        % Restore normal gamma table and close down:
        Screen('CloseAll');
    catch 
        RestoreCluts;
        Screen('CloseAll');
        psychrethrow(psychlasterror);
    end

    % save before averaging
    vals = mean(measurements); % average across repeated measurements
    
    displayRange = range(vals);
    displayBaseline = min(vals);   
    
    %Normalize values
    vals = (vals - displayBaseline) / displayRange;
    inputV = inputV/maxLevel;
    
    if ~exist('fittype');
        fprintf('This function needs fittype() for automatic fitting. This function is missing on your setup.\n');
        fprintf('Therefore i can''t proceed, but the input values for a curve fit are available to you by\n');
        fprintf('defining "global vals;" and "global inputV" on the command prompt, with "vals" being the displayed\n');
        fprintf('values and "inputV" being user input from the measurement. Both are normalized to 0-1 range.\n\n');
        error('Required function fittype() unsupported. You need the curve-fitting toolbox for this to work.\n');
    end
    
    %Gamma function fitting
    g = fittype('x^g');
    fittedmodel = fit(inputV',vals',g);
    displayGamma = fittedmodel.g;
    gammaTable1 = ((([0:maxLevel]'/maxLevel))).^(1/fittedmodel.g); 
    firstFit = fittedmodel([0:maxLevel]/maxLevel);
    
    %Spline interp fitting
    fittedmodel = fit(inputV',vals','splineinterp');
    secondFit = fittedmodel([0:maxLevel]/maxLevel); 
    
    figure;
    plot(inputV, vals, '.', [0:maxLevel]/maxLevel, firstFit, '--', [0:maxLevel]/maxLevel, secondFit, '-.');
    legend('Measures', 'Gamma model', 'Spline interpolation');
    title(sprintf('Gamma model x^{%.2f} vs. Spline interpolation', displayGamma));
    
    %Invert interpolation
    fittedmodel = fit(vals',inputV','splineinterp');
     gammaTable2 = fittedmodel([0:maxLevel]/maxLevel);
     
     
     
     
     
    
return;
