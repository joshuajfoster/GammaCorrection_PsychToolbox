% use PTB's screen function to access the current color lookup table (CLUT)
% JJF 
% 1.22.2017

% open window
s = max(Screen('Screens'));[win,p.sRect] = Screen('OpenWindow',s,0,[120 120 1024,768]);
% access current clut: gammatable is the CLUT
[gammatable] = Screen('ReadNormalizedGammaTable',win); 
% close window
Screen('CloseAll');