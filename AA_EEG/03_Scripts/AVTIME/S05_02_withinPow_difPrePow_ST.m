clear
restoredefaultpath();
restoredefaultpath;
addpath(genpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231'));
addpath('D:\Other\Programs\Matlab\CircStat2012a');
addpath('D:\Experiments\StartTimeExp\AA_EEG\03_Scripts\AVTIME\scripts');
loclog = 'D:\Experiments\StartTimeExp\AA_EEG\02_Data\Logs\';
root_dir = 'D:\Experiments\StartTimeExp\AA_EEG\';
cd('D:\Experiments\StartTimeExp\AA_EEG\03_Scripts\AVTIME\')
vall = dir('*AVTIME*own.txt');
ppnames = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'41';'42';'43'};
AS = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1 2 1];

load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\S01_02_outlierdef'], 'ol');
load('D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Eff_SSVEP', 'freq','int');
loc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\PreProcData\';
AN = [1 2; 2 1];
expname = 'AVTIME';
tpv = [50 175];

%% Last block, resp based on previous power
n = {'incor';'cor'};
%pra = {-20:-1; -15:-1; -10:-1; -7:-1; -5:-1; -3:-1; 0};
%infor.pr = {'PowerPreTrials';'PowerCurrentTrial'};
infor.asna = {'nonassociated';'associated'};
infor.resp = {'incorrect';'correct'};
for pp = 1:length(ppnames)
    load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_' ppnames{pp}], 'FF');
    frinx = nearest(FF.freq, freq(pp));
    it=1;
    for pr = 1:20
        for asna = 1:2
            for resp = 1:2
                trl = find(FF.trialinfo(:,4) == asna-1  & FF.trialinfo(:,7)+1 == resp & FF.trialinfo(:,9) >= 0);
                PPi = [];
                for trc = 1:length(trl)
                    in = find(ismember(FF.trialinfo(:,11)-FF.trialinfo(trl(trc),11),((pr)*-1)) & FF.trialinfo(:,9) >= 0);
                    if length(in) >= 1
                        PPi(end+1,:) = mean(log(FF.pow(in,:,frinx)),1);
                    end
                end
                PP(pp,pr,asna,resp,:) = mean(PPi);
                % particpipants, pre or current power, association tp?,
                % incor/cor
            end
        end
    end
end

%% based on topo and MC cor
% add pr as another dimension in the cluster
%for pr = 1:size(PP,2)
V1.individual = squeeze(PP(:,end:-1:1,1,2,:)-PP(:,end:-1:1,1,1,:));
V1.individual = permute(V1.individual, [1 3 2]);
V1.time = [-1*size(PP,2):1:-1];
V1.dimord = 'subj_chan_time';
load([loc '\eyeblinkcorrectionRMT_AVTIME1'], 'datAud_preproc_RM');
datAud_preproc_RM.label{16}= 'Fp1';
datAud_preproc_RM.label{17}= 'Fp2';
V1.label = datAud_preproc_RM.label;
V2 = V1;
V2.individual = squeeze(PP(:,end:-1:1,2,2,:)-PP(:,end:-1:1,2,1,:));
V2.individual = permute(V2.individual, [1 3 2]);

cfg = [];
cfg.latency = [-20 -1];
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.method = 'montecarlo';
cfg.correctm = 'cluster';
cfg.clusterthreshold = 'nonparametric_individual';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 1;

cfgn = [];
cfgn.method = 'distance';
cfgn.neighbourdist = 0.25;
cfg.neighbours = ft_prepare_neighbours(cfgn, datAud_preproc_RM);
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 1000;
cfg.design(1,1:2*length(ppnames))  = [ones(1,length(ppnames)) 2*ones(1,length(ppnames))];
cfg.design(2,1:2*length(ppnames))  = [1:length(ppnames) 1:length(ppnames)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

st = ft_timelockstatistics(cfg, V1,V2);
%end
save(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\S05_02_PowAccInterStat'],'st');

%% plot the whole shebang
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\S05_02_PowAccInterStat'],'st');
load('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231\template\layout\easycapM1');
cfg = [];
cfg.layout = lay;
cfg.parameter = 'stat';
cfg.highlightchannel = {st.label{sum(st.posclusterslabelmat==1,2) >= 1}};
cfg.highlight = 'on';
ft_topoplotER(cfg,st);
colorbar
set(gcf, 'position',[500 500 375 270]);
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_02_PowAccInterTopo'], 'Color', 'rgb');
print(['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_02_PowAccInterTopo'],'-dpng');

%%
close all
cfg = [];
cfg.layout = lay;
cfg.channel = {st.label{sum(st.posclusterslabelmat==1,2) >= 1}};
cfg.parameter = 'stat';
ft_singleplotER(cfg,st);
ylim = get(gca, 'ylim');
tw = find(sum(st.posclusterslabelmat==1,1));
x = area(st.time(tw), ylim(2).*ones(1,size(tw,2)).*ylim(2));
x.FaceAlpha = 1; x.EdgeAlpha = 1; x.FaceColor = [1 0 0];
set(gcf, 'position',[500 500 375 270]);
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_02_PowAccInter'], 'Color', 'rgb');
print(['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_02_PowAccInter'],'-dpng');

%%
colorbar
%exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_00_PowAccInterTopo_' infor.pr{pr}], 'Color', 'rgb');
%print(['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_00_PowAccInterTopo_' infor.pr{pr}],'-dpng');
