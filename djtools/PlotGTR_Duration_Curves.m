function PlotGTR_Duration_Curves(varargin )
%UNTITLED SummaryTuning of this function goes here
%   Detailed explanation goes here

if nargin==0
    fprintf('\nno input');
    return;
end
newdatadir1=varargin{1};
filename=varargin{2};



        [p,f,ext]=fileparts(filename);
        split=strsplit(f, '_');
        ch=strsplit(split{1}, 'ch');
        channel=str2num(ch{2});
        clust=str2num(split{end});
        outfilename1=sprintf('outPSTH_ch%dc%d.mat',channel, clust) ;
        
        fprintf('\n%s', outfilename1')
        load(outfilename1)

for gd=1:out.numgapdurs
            gapdurs=1:out.gapdurs;
            for rep=1:out.nrepsOFF(gd)
                stOFF=out.M1OFF(gd, 1, rep).spiketimes;
                GTRspikecountOFF=length(find(stOFF>0 & stOFF<76));
                GTRsOFF(gd, rep)=GTRspikecountOFF;
                %                 OFFSETspikecountOFF=length(find(stOFF(-200) & stOFF<0));
                %                 OFFSETsOFF(gd, rep)=OFFSETspikecountOFF
            end
            
            %             test(gd)=sum(OFFSETsOFF(gd, :))
            GTRspikesOFF(gd)= sum(GTRsOFF(gd, :)); %use 'k' for black
            
            for rep=1:out.nrepsON(gd)
                stON=out.M1ON(gd, 1, rep).spiketimes;
                GTRspikecountON=length(find(stON>0 & stON<76));
                GTRsON(gd, rep)=GTRspikecountON;
            end
            
            GTRspikesON(gd)= sum(GTRsON(gd, :));  %use 'b' for blue
            
        end
end

