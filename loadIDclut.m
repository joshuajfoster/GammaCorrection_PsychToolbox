function loadIDclut(s)
% restore the identity CLUT
% 6.14.2017

% open window
[win,p.sRect] = Screen('OpenWindow',s,0,[120 120 1024,768]);
% access current clut: gammatable is the CLUT
LoadIdentityClut(win); 
% close window
Screen('CloseAll');