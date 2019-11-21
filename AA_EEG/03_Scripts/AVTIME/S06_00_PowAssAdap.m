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

load('D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Eff_SSVEP', 'freq','int');
loc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\PreProcData\';
AN = [1 2; 2 1];
expname = 'AVTIME';
tpv = [50 175];

%% extract behavioral response:
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_all'], 'AV');
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_1'], 'FF');
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Behdata.mat'], 'ppnames', 'vars');
ppto = [1:40];
clear proplp amtr
allfr = {[1:100], [102:201]};
allfrab17 = {[1:201], [203:400]};
tp = [50 175];
bl = [4];
for pp = 1:length(ppto)
    for t = 1:2
        for fr = 1:2
            inx = find(vars(ppto(pp)).SOUND2 == fr & vars(ppto(pp)).TIMEPOINT == tp(t)& ismember(vars(ppto(pp)).BLOCK, bl));
            amlp = find(vars(ppto(pp)).RESPONSE(inx) == AN(AS(ppto(pp)),fr));
            proplp(pp,fr,t) = length(amlp)./length(inx);
            amtr(pp,fr,t) = length(inx);
        end
    end
end
dif = proplp(:,:,2)-proplp(:,:,1);
difN = dif(:,2)-dif(:,1); % association effect

%% remove influential cases based on cooks distance
ch = [18];
afinx = arrayfun(@(x) nearest(FF.freq, freq(x)), 1:40);
cor2 = arrayfun(@(x,y) squeeze(AV.savepowvl(x,ch,y)), 1:40, afinx, 'uniformoutput', 1)';
es = difN;
tb = table(cor2, es);
restoredefaultpath();
for ch = 1:length(FF.label)
    tb{:,1} = arrayfun(@(x,y) squeeze(AV.savepowvl(x,ch,y)), ppto, afinx, 'uniformoutput', 1)';
    lm = fitlm(tb, 'cor2~1+es');
    IC(ch,:) = lm.Diagnostics.CooksDistance > 4/length(ppto);    
end
ppto = find(sum(IC) <= 5);

%% plot with and without removed influential cases
addpath(genpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231'));
figure
ch = [18];
col = double(sum(IC) <= 5)';
scatter(es, cor2, 200,col, '.');
t = zeros(3,10);
t(1:3,1:5) = 0.7;
colormap(t');
[r p] = corr(es(ppto), cor2(ppto));
b = polyfit(es(ppto),cor2(ppto),1);
hold on
plot([-0.5 0.6], b(1).*[-0.5 0.6]+b(2), 'r');
xlabel('Association effect');ylabel('log power');
text(0.2, 1, ['r=' num2str(round(r, 2))]);
text(0.2, 0.8, ['p=' num2str(round(p, 3))]);
set(gca, 'xlim', [-0.6 0.7]);
set(gcf, 'position', [300 300 280 300]);
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S06_00_PowAssAdapFz'], 'Color', 'rgb');
print(['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S06_00_PowAssAdapFz'],'-dpng')

%% based on topo and MC cor
restoredefaultpath();
addpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231');
ppto=1:40;
corA = arrayfun(@(x,y) squeeze(AV.savepowvl(x,:,y))', ppto, afinx(ppto), 'uniformoutput', 0)';
clear V1 V2 st
V1.individual(:,:,1) = cell2mat(corA')';
V1.time = [0];
V1.dimord = 'subj_chan_time';
V1.label = FF.label;

load([loc '\eyeblinkcorrectionRMT_AVTIME1'], 'datAud_preproc_RM');  
cfg = [];
cfg.latency = [-1 1];
cfg.statistic = 'ft_statfun_correlationT_SO';
cfg.method = 'montecarlo';
cfg.correctm = 'cluster';
cfg.clusterthreshold = 'nonparametric_individual'; 
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 1;

cfgn = [];
cfgn.method = 'distance';
datAud_preproc_RM.label{16}= 'Fp1';
datAud_preproc_RM.label{17}= 'Fp2';
V1.label = datAud_preproc_RM.label;

cfgn.neighbourdist = 0.25;
cfg.neighbours = ft_prepare_neighbours(cfgn, datAud_preproc_RM);
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 1000;
cfg.design = es(ppto);

st = ft_timelockstatistics(cfg, V1);
%ave(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\S06_00_PowAssAdapStat_v2'],'st');

%% do the topoplot
addpath(genpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231'));
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\S06_00_PowAssAdapStat_v2'],'st');
load('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231\template\layout\easycapM1');

V1.avg(:,1) = mean(V1.individual);
V1.avg(:,2) = mean(V1.individual);
V1.dimord = 'chan_time';
subplot(121)
cfg = [];
cfg.highlightchannel = find(st.negclusterslabelmat == 1);
cfg.highlightchannel = find(st.prob <= 0.05);
cfg.highlightsymbol = '*';
cfg.highlight = 'off';
cfg.highlightcolor = [1 0 1];
cfg.xlim = [-0.7 0.7];
cfg.comment = 'no';
cfg.layout = lay;
ft_topoplotER(cfg, V1);
colorbar

V1.avg(:,1) = st.rho;
V1.avg(:,2) = st.rho;
V1.dimord = 'chan_time';
subplot(122)
cfg = [];
cfg.highlightchannel = find(st.negclusterslabelmat == 1);
cfg.highlightchannel = find(st.prob <= 0.05);
cfg.highlightsymbol = '*';
cfg.highlight = 'on';
%cfg.highlightcolor = [1 0 1];
cfg.xlim = [-0.5 0.5];
cfg.zlim = [-0.4 0.4];
cfg.layout = lay;
cfg.comment = 'no';
ft_topoplotER(cfg, V1);
colorbar
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S06_00_PowAssAdapTopo'], 'Color', 'rgb');
print(['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S06_00_PowAssAdapTopo'],'-dpng')
