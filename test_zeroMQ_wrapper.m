%little utility script to see if the zeroMQwrapper mex file is working
%properly. If you get errors, you might need to check whether the right mex
%file for your platform is on the path. You could also try executing this
%from within the same directory as that mex file. If that doesn't work you
%might need to compile the cpp file for your platform (see installation
%instructions).

%url = 'localhost:5556'; % or, e.g., //'tcp://10.71.212.19:5556 if GUI runs on another machine...
url = 'tcp://127.0.0.1:5556';
%url = 'tcp://184.171.84.30:5556';    % if GUI runs on another machine...
fprintf('\n\nopening thread to %s', url)

handle = zeroMQwrapper('StartConnectThread',url);
pause(.5)
fprintf('\nsuccess, handle is %g', handle)
fprintf('\ntrying to send some messages...')
zeroMQwrapper('Send',handle ,'ChangeDirectory');
pause(.1)
zeroMQwrapper('Send',handle ,'NewDesign nGo_Left_Right');
pause(.1)
zeroMQwrapper('Send',handle ,'AddCondition Name GoRight TrialTypes 1 2 3');
pause(.1)
zeroMQwrapper('Send',handle ,'AddCondition Name GoLeft TrialTypes 4 5 6');
pause(.1)


% indicate trial type number #2 has started
zeroMQwrapper('Send',handle ,'TrialStart 2');
pause(.1)

% indicate that trial has ended with outcome "1"
zeroMQwrapper('Send',handle ,'TrialEnd 1');
pause(.1)

zeroMQwrapper('CloseThread',handle);

fprintf('\nIf you didn''t see any errors, zeroMQwrapper is probably working fine. ')
fprintf('\nNote that this does not check to see if open-ephys actually received ')
fprintf('\nany of those messages. For that, check the open-ephys message console.')
fprintf('\nIf you saw "Connection was closed. Cannot send", it means that')
fprintf('\nthe zeroMQwrapper file probably works, but there was a problem')
fprintf('\nopening the connection thread. In this case, try it again.\n')