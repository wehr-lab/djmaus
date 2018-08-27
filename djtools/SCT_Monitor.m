function SCT_Monitor(datadir, StartAcquisitionSec, Events, all_channels_data, all_channels_timestamps, all_channels_info)

%plots soundcard triggers, timestamps, and network messages as a sanity
%check and way to monitor the data processing for a given data file
%this used to be a commented-out stanza in each processing function but I'm
%pulling it out for modularity
%mw 1-30-2017

fprintf('\nrunning SCT_Monitor to examine soundcard triggers, timestamps, and network messages...')

%   I'm running the soundcard trigger (SCT) into ai1 as another sanity check.
SCTfname=getSCTfile(datadir);
if isempty(SCTfname)
    warning('could not find soundcard trigger file')
else
    [SCTtrace, SCTtimestamps, SCTinfo] =load_open_ephys_data(SCTfname);
end

%     %here I'm loading a data channel to get good timestamps - the ADC timestamps are screwed up
%     datafname=getContinuousFilename( pwd, 1 );
%     [scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(datafname);

%SCTtimestamps=SCTtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition


%sanity check
fid=fopen('temp.txt', 'w');
for i=1:length(all_channels_timestamps)
    fprintf(fid, '%f, eventType %d, Id %d, channel %d\n', ...
        all_channels_timestamps(i), all_channels_info.eventType(i), ...
        all_channels_info.eventId(i), all_channels_data(i));
end
fclose(fid);



%messages is a list of all network event, which includes the stimuli
%messages sent by djmaus, as well as the "ChangeDirectory" and
%"GetRecordingPath" messages sent by djmaus, as well as 2 initial system
%messages. I strip out the stimulus (sound) event and put them in "Events."
%Events is a list of sound event, which were sent by djmaus with the
%'TrialType' flag.

%sanity check, continued
% fid=fopen('temp.txt', 'a');
% for i=1:length(Events)
%     fprintf(fid, '%f, %s, freq %f\n', ...
%         Events(i).message_timestamp_sec, Events(i).type, Events(i).frequency);
% end
% fclose(fid);





figure
hold on
%set(gcf, 'pos', [-1853 555 1818 420]);
SCTtrace=SCTtrace./max(abs(SCTtrace));
%     scaledtrace=scaledtrace./max(abs(scaledtrace));
plot(SCTtimestamps, SCTtrace, 'b') %plotting with data timestamps to work around wierd bug in ADC timestamps
% plot(datatimestamps, scaledtrace, 'm') %plotting with data timestamps to work around wierd bug in ADC timestamps


%let's plot some other signals too
stimfile=getStimfile(datadir); %mw 08.30.2107 old: sprintf('%s_ADC2.continuous', node);
laserfile=getLaserfile(datadir); %mw 08.30.2017 old: sprintf('%s_ADC2.continuous', node);
[stimtrace, stimtimestamps, stiminfo] =load_open_ephys_data(stimfile);
stimtimestamps=stimtimestamps-StartAcquisitionSec;
[lasertrace, lasertimestamps, laserinfo] =load_open_ephys_data(laserfile);
lasertimestamps=lasertimestamps-StartAcquisitionSec;
plot(stimtimestamps, .1*stimtrace, 'm') 
plot(SCTtimestamps, .1*lasertrace, 'c') 
title('SCT Monitor sanity check')


hold on
fprintf('plotting %d Events... ', length(Events))
%plot "software trigs" i.e. network messages in red o's
for i=1:length(Events)
    plot(Events(i).message_timestamp_sec, .25, 'ro');
    plot(Events(i).soundcard_trigger_timestamp_sec, 1, 'g*');
    text(Events(i).message_timestamp_sec, .5, sprintf('network message #%d', i))
    text(Events(i).message_timestamp_sec, .75, sprintf('%s', Events(i).type))
    text(Events(i).soundcard_trigger_timestamp_sec, 1.25, sprintf('SCT #%d', i))
end

%all_channels_info.eventType(i) = 3 for digital line in (TTL), 5 for network Events
%all_channels_info.eventId = 1 for rising edge, 0 for falling edge
%all_channels_data(i) is the digital input line channel


% % plot TTL SCTs in green ^=on, v=off
% for i=1:length(all_channels_timestamps)
%     if all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==1 & all_channels_data(i)==2
%         plot(all_channels_timestamps(i), 1, 'g^')
%         text(all_channels_timestamps(i), 1, 'TTL on/off')
%     elseif all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==0 & all_channels_data(i)==2
%         plot(all_channels_timestamps(i), 1, 'gv')
%     end
% end
% plot TTL SCTs in green ^=on, v=off
for i=1:length(all_channels_timestamps)
    if all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==1 & all_channels_data(i)==2
        plot(all_channels_timestamps(i)-StartAcquisitionSec, 1, 'g^')
        text(all_channels_timestamps(i)-StartAcquisitionSec, 1, 'TTL on/off')
    elseif all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==0 & all_channels_data(i)==2
        plot(all_channels_timestamps(i)-StartAcquisitionSec, 1, 'gv')
    end
end


for i=1:length(Events)
    xlim([Events(i).message_timestamp_sec-2 Events(i).message_timestamp_sec+3])
    ylim([-5 2])
    pause(.002)
   % pause
end

%Plot successive timestamp values in order to see if they skip ahead or
%double back, which would cause a problem in plotting any data. TH
%2017-11-21
figure
hold off
plot(SCTtimestamps, 'b')
xlabel("Index")
ylabel("SCT time stamp values")
title("Time Stamps. -Look for deviations from straight line.")

