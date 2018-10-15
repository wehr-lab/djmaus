function PlotGPIASflashtrain_PSTH_single(varargin)

%plots a single file of clustered spiking GPIAS flashtrain data from djmaus
%this has a train of laser pulses with varying pulewidth
%
% usage: PlotflashtrainGPIAS_PSTH(datapath, t_filename, [xlimits],[ylimits], [binwidth])
% (xlimits, ylimits, binwidth are optional)
%
%Processes data if outfile is not found;

rasters=1;
force_reprocess=1;

if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};

t_filename=varargin{2};

try
    xlimits=varargin{3};
catch
    xlimits=[];
end
try
    ylimits=varargin{4};
catch
    ylimits=[];
end
try
    binwidth=varargin{5};
catch
    binwidth=5;
end

if force_reprocess
    fprintf('\nForce re-process\n')
    ProcessGPIASflashtrain_PSTH_single(datadir,  t_filename, xlimits, ylimits, binwidth);
end

[p,f,ext]=fileparts(t_filename);
split=strsplit(f, '_');
ch=strsplit(split{1}, 'ch');
channel=str2num(ch{2});
clust=str2num(split{end});

outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
fprintf('\nchannel %d, cluster %d', channel, clust)
fprintf('\n%s', t_filename)
fprintf('\n%s', outfilename)

cd(datadir)

if exist(outfilename,'file')
    load(outfilename)
    fprintf('\nloaded outfile')
else
    fprintf('\ncould not find outfile, calling ProcessGPIASflashtrain_PSTH_single...')
    ProcessGPIASflashtrain_PSTH_single(datadir,  t_filename, xlimits, ylimits, binwidth);
    load(outfilename);
end

%if xlimits don't match, force preprocess
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        fprintf('\nPlot called with xlimits [%d %d] but xlimits in outfile are [%d %d], calling ProcessAsymGPIAS_PSTH_single...', xlimits(1), xlimits(2), out.xlimits(1), out.xlimits(2))
        ProcessGPIASflashtrain_PSTH_single(datadir,  t_filename, xlimits, ylimits, binwidth);
        load(outfilename);
    end
end

%we have the following dimensions:
% M1(numgapdurs, numpulsewidths, nreps).spiketimes
% mM1(numgapdurs, numpulsewidths).spiketimes
% mM1spikecount(numgapdurs, numpulsewidths)
%note that numpulsewidths includes 0, the laser-off condition
%thus there is no M1 OFF

numpulsewidths=out.numpulsewidths;
numgapdurs=out.numgapdurs;
pulsewidths=out.pulsewidths;
gapdurs=out.gapdurs;
gapdelay=out.gapdelay;
samprate=out.samprate; %in Hz
mM1=out.mM1;
M1=out.M1;
nreps=out.nreps;
soa=out.soa;
if isfield(out, 'LaserRecorded')
    LaserRecorded=out.LaserRecorded;
    M1Laser=out.M1Laser; % a crash here means this is an obsolete outfile. Set force_reprocess=1 up at the top of this mfile. (Don't forget to reset it to 0 when you're done)
    mM1Laser=out.mM1Laser;
else
    LaserRecorded=0;
end
if isfield(out, 'StimRecorded')
    StimRecorded=out.StimRecorded;
    M1Stim=out.M1Stim;
    mM1Stim=out.mM1Stim;
else
    StimRecorded=0;
end

fs=10; %fontsize

% %find optimal axis limits
if isempty(xlimits)
    xlimits(1)=-1.5*max(gapdurs);
    xlimits(2)=2*soa;
end
fprintf('\nusing xlimits [%d-%d]', xlimits(1), xlimits(2))

if isempty(ylimits)
    ymax=0;
    for gdindex=1:numgapdurs
        for pwindex=1:numpulsewidths
            st=mM1(gdindex,pwindex).spiketimes;
            nr=nreps(gdindex,pwindex);
            X=xlimits(1):binwidth:xlimits(2); %specify bin centers
            [N, x]=hist(st, X);
            N=N./nr; %normalize to spike rate (averaged across trials)
            N=1000*N./binwidth; %normalize to spike rate in Hz
            ymax= max(ymax,max(N));
        end
    end
    ylimits=[-.3 ymax];
end

M1Lasertemp=M1Laser(find(M1Laser~=0));
minlaser=min(M1Lasertemp(:));
maxlaser=max(M1Lasertemp(:));
nM1Laser=M1Laser-minlaser;
nM1Laser=nM1Laser./abs(maxlaser);

% figure('position',[161 1669 600 750])
% p=0;
% subplot1(numpulsewidths, numgapdurs, 'Max', [.95 .96], 'Gap', [.05 .01 ])
% for pwindex=1:numpulsewidths
%     for gdindex=1:numgapdurs
%         p=p+1;
%         subplot1(p)
%         hold on
%         
%         for rep=1:nreps(gdindex,pwindex)
%             Lasertrace=squeeze(nM1Laser(gdindex,pwindex,rep, :));
%             %Lasertrace=Lasertrace -mean(Lasertrace(1:100));
%             % Lasertrace=Lasertrace-min(M1Laser(:));
%             % Lasertrace=Lasertrace./max(M1Laser(:));
%             Lasertrace=5*Lasertrace;
%             plot( t, Lasertrace, 'c')
%         end
%     end
% end
    
%plot the tuning curve
figure('position',[257   144 600 750])

p=0;
subplot1(numpulsewidths, numgapdurs, 'Max', [.95 .96], 'Gap', [.05 .01 ])
    for pwindex=1:numpulsewidths
for gdindex=1:numgapdurs
        
        p=p+1;
        subplot1(p)
        hold on
        spiketimes1=mM1(gdindex,pwindex).spiketimes; %spiketimes are in ms relative to gap termination
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(spiketimes1, X);
        N=N./nreps(gdindex,pwindex); %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        offset=0;
        yl=ylimits;
        inc=(yl(2))/max(nreps(:));
        if rasters==1
            for n=1:nreps(gdindex,pwindex)
                spiketimes2=M1(gdindex,pwindex, n).spiketimes;
                offset=offset+inc;
                h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
            end
        end
        bar(x, N,1,'facecolor','k','edgecolor','k');
        line([0 0],[ylim],'color','m')
        line(-[(gapdurs(gdindex)) (gapdurs(gdindex))],[ylim],'color','m')
        line(xlimits, [0 0], 'color', 'k')
        ylimits2(2)=ylimits(2)*2.2;
        ylimits2(1)=-2;
        ylim(ylimits2)
        
        xlim(xlimits)
        set(gca, 'fontsize', fs)
        %set(gca, 'xticklabel', '')
        %set(gca, 'yticklabel', '')
        if StimRecorded
            Stimtrace=squeeze(mM1Stim(gdindex,pwindex, :));
            Stimtrace=Stimtrace -mean(Stimtrace);
            height=.1*diff(ylimits2); %reasonable height
            %height=.5*diff(ylimits2); %magnified height for scrutinizing stimulus
            Stimtrace=height*Stimtrace;
            t=1:length(Stimtrace);
            t=1000*t/out.samprate; %convert to ms
            t=t+out.xlimits(1); %correct for xlim in original processing call
            offset=-range(Stimtrace);
            ylimits2(1)=ylimits2(1)+2*offset;
            %  offset=ylimits(1)+.1*diff(ylimits);
            plot(t, Stimtrace+offset, 'm')
            ylim(ylimits2)
        else
            %do nothing
        end
        if LaserRecorded
            height=.5*diff(ylimits2);
            offset=-height;
            for rep=1:nreps(gdindex,pwindex)
                Lasertrace=squeeze(nM1Laser(gdindex,pwindex,rep, :));
                % Lasertrace=Lasertrace -mean(Lasertrace(1:100));
               % Lasertrace=Lasertrace-min(M1Laser(:));
               % Lasertrace=Lasertrace./max(M1Laser(:));
                Lasertrace=height*Lasertrace;
                plot( t, Lasertrace, 'c')
            end
            ylimits2(1)=ylimits2(1)-2*range(Lasertrace);
            ylim(ylimits2)
        end
        
        text(xlimits(1)-10, ylimits(2), sprintf('%d\n%.1f', gapdurs(gdindex), pulsewidths(pwindex)), 'color', 'r')
        if gdindex<numgapdurs
            set(gca, 'yticklabel', '');
        end
    end
end

subplot1(1)
h=title(sprintf('%s: \ntetrode%d cell %d, nreps: %d-%d',datadir,channel,out.cluster,min(nreps(:)),max(nreps(:))));
set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')



%label amps and freqs
% p=0;
% for gdindex=1:numgapdurs
%     p=p+1;
%     subplot1(p)
%     vpos=ylimits(2);
%     text(xlimits(1), vpos, sprintf('%d', gapdurs(gdindex)), 'color', 'r')
%     if gdindex<numgapdurs
%         set(gca, 'yticklabel', '');
%     end
% end
%turn on ytick for bottom-most plot
set(gca, 'yticklabelmode', 'auto');








