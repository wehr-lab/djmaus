function djPrefs
%Local preferences
%specifics for the local machine are at top
%preferences for individual users at the bottom
global pref SP

% machine-wide default prefs:
pref.users={'lab','mike-home','mike-kombi','mike-linux','test','test2','kip','kat','Hanna','apw','Emily','Nick','Tai'};
%try
    pref.soundcarddeviceID=GetAsioLynxDevice; %note: using the asio device for both input and output
%catch
 %   pref.soundcarddeviceID=40;
%end
pref.num_soundcard_outputchannels=4;
pref.Soundchannel1=1;
pref.Soundchannel2=nan; %we're not doing binaural on this rig at the moment
pref.Soundcardtriggerchannel=3;
pref.Laserchannel=4;
pref.Shockchannel=2;
pref.CampulseOn=nan; %we're not using camera pulses on this rig
pref.CampulseOff=nan; %we're not using camera pulses on this rig
pref.camera_rec=nan; %we're not using camera pulses on this rig
pref.reqlatencyclass=0; %Rig1: saw dropouts with 0, fewer 1, fewer with 2, still some with 3
pref.suggestedLatency=.1;
pref.SoundFs=192000;
pref.maxSPL=80;
pref.allmouseIDs='';
pref.root=fileparts(which(mfilename));
pref.windowpos=[940    62   420   660];  %djmaus GUI position [88    62   363   660]

pref.local =0; %1 for local communication (djmaus and open-ephys on same
%computer), 0 for remote (djmaus and open-ephys on different computers)
if pref.local %same computer
    switch computer
        case 'MACI64'
            pref.zmqurl='tcp://localhost:5556'; %seems to work for mac
            pref.mexpath='mac';
        case 'GLNXA64'
            pref.zmqurl='tcp://127.0.0.1:5556'; %seems to work for linux
            pref.mexpath='unix';
        case 'PCWIN64'
            pref.zmqurl='tcp://127.0.0.1:5556'; %seems to work for windows
            pref.mexpath='windows';
    end
    pref.datapath=fullfile(pref.root, 'Data\lab'); %open-ephys data is acquired on this computer
    pref.remotedatapath=pref.datapath; %these are the same thing for local acquisition
        pref.datahost=''; %changed from '\', mw 080217
else %different computer
    switch computer
        case 'MACI64'
            pref.mexpath='mac';
        case 'GLNXA64'
            pref.mexpath='unix';
        case 'PCWIN64'
            pref.mexpath='windows';
    end
    %specific zmq url for the open-ephys computer
%    pref.zmqurl='tcp://184.171.85.38:5556';
    pref.zmqurl='tcp://184.171.84.72:5556';
%    pref.datahost='\\wehrrig1b';
    pref.datahost='o:'; %mw08-02-17
    pref.remotedatapath='d:\lab\djmaus\Data\lab'; %what the open-ephys datapath looks like on the other computer
   % pref.datapath='\\wehrrig1b\d\lab\djmaus\Data\lab'; %what the open-ephys datapath looks like to this computer
   %switching to o: mw08-02-17 
   pref.datapath='o:\lab\djmaus\Data\lab'; %what the open-ephys datapath looks like to this computer
end
pref.stimuli=fullfile(pref.root, 'stimuli');
%not currently used since "launch OE" button still buggy:
pref.OEpath='D:\lab\plugin-GUI\Builds\VisualStudio2013\x64\Release64\bin\open-ephys.exe &'; %make sure to put & at the end

% individual user prefs:
if ~isfield(SP, 'user')
    SP.user='lab'; %default user
end

try
    switch SP.user
case 'Tai'
pref.datapath='\\wehrrig1b\d\lab\djmaus\Data\Tai';
pref.remotedatapath='d:\lab\djmaus\Data\Tai';
case 'Nick'
%pref.datahost='n:';
 pref.datahost='\\wehrrig1b\';
 pref.datapath='\\wehrrig1b\e'; %mw 10122018
% pref.datapath='n:\Nick';
pref.remotedatapath='e:\NickNick';
case 'Emily'
pref.datapath='\\wehrrig1b\d\lab\djmaus\Data\Emily';
pref.remotedatapath='d:\lab\djmaus\Data\Emily';
case 'apw'
pref.datapath='\\wehrrig1b\d\lab\djmaus\Data\apw';
pref.remotedatapath='d:\lab\djmaus\Data\apw';
case 'Hanna'
pref.datapath='\\wehrrig1b\d\lab\djmaus\Data\Hanna';
pref.remotedatapath='d:\lab\djmaus\Data\Hanna';
case 'kat'
pref.datapath='\\wehrrig1b\d\lab\djmaus\Data\kat';
pref.remotedatapath='d:\lab\djmaus\Data\kat';
case 'kip'
pref.datapath='\\wehrrig1b\d\lab\djmaus\Data\kip';
pref.remotedatapath='d:\lab\djmaus\Data\kip';


        case 'mike-home'
            %pref.stimuli='/home/mike/djmaus/stimuli/mike';
            pref.datapath='/Users/mikewehr/Documents/Data';
        case 'mike-linux'
            pref.stimuli='/home/mike/djmaus/stimuli/mike';
            pref.datapath='/home/Data';
        case 'mike-kombi'
            pref.stimuli='/Users/mikewehr/Documents/Analysis/djmaus-master/stimuli/mike';
            pref.datapath='/Users/mikewehr/Documents/Data2';
            pref.soundcarddeviceID=3;
    end
end



%the following sections are generated automatically, please don't edit below here








%saved mouseIDs for lab
switch SP.user
case  'lab'
	pref.allmouseIDs={'','7092','7093','7094','7165','7305','7306','7315','7400','7400-light','7400-noise','7401','7401-light','7401-noise','7437','7515','7516','7681','7682','7684','7700','7730','7731','7733','7744','7747','7748','7750','7774','7812','7815','7817','7818','7856','7887','7890','7984','7test','8211','8282','8283','8287','8288','8568','A','A1IC1','A1IC2','AC','ACIC1','ACIC2','ACIC3','Fred','LGN','MGLG','N003','N010','N011','NTSR','NTSRa','NTSRb','PV3','pv3','test','test64'};
end


































%saved mouseIDs for kip
switch SP.user
case  'kip'
	pref.allmouseIDs={'','1111','3B13A','3B13B','7165','7305','7306','7315','7436','7437','8092','8093','8094','8096','8281','8282','8283','N00','N002','N003','N004','N006','N007','N008','N009','N010','N011','N012','N013','N014','N015','PV3','XXXX','test'};
end



















%saved mouseIDs for kat
switch SP.user
case  'kat'
	pref.allmouseIDs={'','3B13A','3B13B','3b13b','7306','7436'};
end







%saved mouseIDs for Hanna
switch SP.user
case  'Hanna'
	pref.allmouseIDs={'','3B13A','7515','7516','HH1','HH2','HH3'};
end





























































































































































































%saved mouseIDs for apw
switch SP.user
case  'apw'
	pref.allmouseIDs={'','54','55','56','57','7437','7515','7516','7681','7682','7684','7744','7747','7851','7853','7856','7857','7860','7861','7863','7865','NTSRa','NTSRb','ntsr2','ntsra','ntsrb','pv5111lt','pv511ln','pv511lt','pv511rn','pv511rt','pv511tt','pv514lt','pv514nt','pv514rt','pv514tt','pv62rt','test'};
end













































































































































































%saved mouseIDs for Emily
switch SP.user
case  'Emily'
	pref.allmouseIDs={'','7684','7744','7745','7747','7748'};
end



%saved mouseIDs for Nick
switch SP.user
case  'Nick'
	pref.allmouseIDs={'','6666','7092','7093','7094','7165','7305','7306','7315','7400','7400-light','7400-noise','7401','7401-light','7401-noise','7437','7515','7516','7681','7682','7684','7700','7730','7731','7733','7744','7747','7748','7750','7774','7812','7815','7817','7818','7856','7887','7890','7984','7test','8211','8282','8283','8287','8288','8568','8888','A','A1IC1','A1IC2','AC','ACIC1','ACIC2','ACIC3','Fred','LGN','MGLG','N003','N010','N011','NTSR','NTSRa','NTSRb','PV3','pv3','test','test64'};
end



































































%saved mouseIDs for Tai
switch SP.user
case  'Tai'
	pref.allmouseIDs={'','7697','7698','7766'};
end




































































































































































































































































































































































































































































































































































































































































































