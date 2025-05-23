function [bad_channels, magnitude, channel, shank]=ReadImpedanceTestXML
%reads impedance test xml file generated by open ephys, for use with silicon
%probes to check for bad channels.%mw 3.17.2025
%at the moment it just reads the first '*.xml' file found in the current
%directory. We could upgrade it to find the
%most recent relevant impedance file, either automatically or manually.
%Note that we have to correct the channel numbers in the xml file accroding to the
% channel map file chanMap128.mat which should be in the DataRoot
% directory, and we use DataRoot from dirs.mat to find the channel map
% file.
%
%usage:
%  [bad_channels, magnitude, channel, shank]=ReadImpedanceTestXML
%  looks in current directory and loads the first xml file it finds
bad_channels=[];

d=dir('*.xml');
if isempty(d) fprintf('\n%s: could not find impedance test file', mfilename), return, end
try
    load dirs
    load(fullfile(DataRoot, 'chanMap128.mat'))
catch
    error('cannot find chanMap128.mat file, which is required to interpret channel numbers')
end

%here you could find the most recent file relative to youre recording session using d.date
filename=d(1).name;
s=xml2struct(filename);

chan=0;
for sh=1:2
    for ch=1:64
        chan=chan+1;
        shank(chan)=sh;
        channel(chan)=chan;
        magnitude(chan)=str2num(s.IMPEDANCES.HEADSTAGE{sh}.CHANNEL{ch}.Attributes.magnitude);
        name{chan}=s.IMPEDANCES.HEADSTAGE{sh}.CHANNEL{ch}.Attributes.name;
        number(chan)=str2num(s.IMPEDANCES.HEADSTAGE{sh}.CHANNEL{ch}.Attributes.number);
        phase(chan)=str2num(s.IMPEDANCES.HEADSTAGE{sh}.CHANNEL{ch}.Attributes.phase);
    end

end

bad_chan_thresh=2*std(magnitude);
bad_chan_thresh=1e6; %1 MOhm
bad_channels=sort(chanMap(find(magnitude>bad_chan_thresh)));

figure
subplot(211)
stem(chanMap, magnitude)
set(gca, 'yscale', 'log')
ylabel('magnitude (Ohms)')
line(xlim, [1 1]*bad_chan_thresh)
title(['impedance from file ', filename])
text(chanMap, magnitude, int2str(chanMap'))

subplot(212)
stem(chanMap, phase)
ylabel('phase')
xlabel('channel')
text(chanMap, phase, int2str(chanMap'))
set(gcf, 'pos',[192 818 2140 420]);

fprintf('\n%s: ', mfilename)
fprintf('\nreading file %s from %s', filename, d(1).date)
fprintf('\nusing channel map %s', fullfile(DataRoot, 'chanMap128.mat'))
fprintf('\nusing magnitude threshold of %.1f MOhm', bad_chan_thresh/1e6)
fprintf('\ndetected bad channels: ')
fprintf('%d ', bad_channels)



