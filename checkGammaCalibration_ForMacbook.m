function checkGammaCalibration_v2(numMeasures, screenid)
% Check that the gamma correction was successful. numMeasures specifies
% the number of observations to take (default = 9). Must set path to find
% calibration file (line 12).

% Adapted from PTB's CalibrateMonitorPhotometer function.
%
% Joshua J. Foster
% December, 29, 2016
%
% v2 update (7.9.2018)
% script now displays input, we can edit if we need to (JJF made changes)
% levels presented in random order (Janna Wennberg made changes)

Screen('Preference', 'SkipSyncTests', 1);

% specify location of the gamma table
load('MyMacbookPro_GammaCorrection.mat')
gammaTable = calib.gammaTable2;

KbName('UnifyKeyNames');

if (nargin < 1) || isempty(numMeasures)
    numMeasures = 17;
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
        screenid = 1;
    end
    
    % Open black window:
    win = Screen('OpenWindow', screenid, 0);
    maxLevel = Screen('ColorRange', win);
        
    % load in gamma table
    %     LoadIdentityClut(win);
    old = Screen('LoadNormalizedGammaTable',win,gammaTable*[1 1 1])
    
    vals = nan(1,numMeasures);
    inputV = [0:(maxLevel+1)/(numMeasures - 1):(maxLevel+1)]; 
    inputV(end) = maxLevel;
    shuffledInputV = Shuffle(inputV);
    for i = 1:numMeasures
        
        %----------------------------------------------------------
        % get a response
        %----------------------------------------------------------
        level = shuffledInputV(i);
        
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
        k = inputV == level;       
        vals(k) = resp; % save to vals
    end
        
    % Restore normal gamma table and close down:
    LoadIdentityClut(win)
    Screen('CloseAll');
catch %#ok<*CTCH>
%     LoadIdentityClut(win)
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
    
displayRange = range(vals);
displayBaseline = min(vals);

%Normalize values
norm_vals = (vals - displayBaseline) / displayRange;
norm_inputV = inputV/maxLevel;

%calculate residuals
resids = norm_vals-norm_inputV;

% save file
fName = ['GammaCalTest.mat']
save(fName,'vals','inputV','norm_vals','norm_inputV','resids');

% show values scattered around y = x.
figure;
plot(norm_inputV,norm_inputV,'--k'); hold on;
scatter(norm_inputV,norm_vals,'r','filled');
ylabel('Luminance (normalized)')
xlabel('Specified intensity (normalized)')

% plot residuals
figure;
plot(norm_inputV,zeros(1,length(norm_inputV)),'--k'); hold on;
scatter(norm_inputV,resids,'r','filled');
ylabel('Residual')
xlabel('Specified intensity (normalized)')
ylim([-0.1 0.1])


return;