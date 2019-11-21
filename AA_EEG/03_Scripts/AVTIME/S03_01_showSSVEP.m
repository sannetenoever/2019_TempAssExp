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
saveloc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\';
load([saveloc 'S01_02_outlierdef'], 'ol');
varnames = {'accuracy';'sound';'ass';'block';'pp';'curval'; 'dev'; 'pow'};
bl = [0:5];
pr=1;
pra = {-7:-1, 0};
infor.pr = {'PowerPreTrials';'PowerCurrentTrial'};
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\S05_00_PowAccInterStat' infor.pr{pr}],'st');
ch = find(st.posclusterslabelmat==1);

for pp = 1:length(ppnames)
    load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_' ppnames{pp}], 'FF');
    load([loc '\eyeblinkcorrectionRMT_' expname ppnames{pp}], 'datAud_preproc_RM');
    load([loc '\trialinfo_' expname ppnames{pp}], 'trlinf');
    trlinf = trlinf(setdiff(1:size(trlinf,1), datAud_preproc_RM.badtr),:);
    
    % align to auditory    
    in = find(ismember(trlinf(:,9), bl));
    if str2double(ppnames{pp}) < 42
        cfg = [];
        cfg.trials = in;
        cfg.offset = round(-0.110.*datAud_preproc_RM.fsample);
        datAud_preproc_RM = ft_redefinetrial(cfg, datAud_preproc_RM);
    else % different due to computer update
        cfg = [];
        cfg.trials = in;
        cfg.offset = 0;
        datAud_preproc_RM = ft_redefinetrial(cfg, datAud_preproc_RM);
    end
    datAud_preproc_RM.trialinfo = trlinf;
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [1 20];
    x = ft_preprocessing(cfg, datAud_preproc_RM);    
    cfg = [];
    cfg.offset = zeros(size(trlinf,1),1);
    cfg.offset(trlinf(:,2) == 175) = 256*(0.5/freq(pp));
    x = ft_redefinetrial(cfg, x);
    cfg = [];
    cfg.toilim = [-2 0];
    x = ft_redefinetrial(cfg, x);
    cfg = [];
    cfg.keeptrials = 'yes';
    cfg.vartrllength = 2;
    x  = ft_timelockanalysis(cfg, x);
    
    % get pretrial power:
    frinx = nearest(FF.freq, freq(pp));
    PPi = [];
    clear PPi SSVEPtr TimeSSVEPtr
    for trc = 1:length(FF.trialinfo)
        in = find(ismember(FF.trialinfo(:,11)-FF.trialinfo(trc,11),pra{pr}) & FF.trialinfo(:,9) >= 0);
        PPi(trc) = squeeze(mean(mean(log(FF.pow(in,ch,frinx)),1),2));
        SSVEPtr(trc,:) = squeeze(mean(mean(log(FF.pow(in,ch,:)),1),2));
        TimeSSVEPtr(trc,:) = squeeze(mean(mean(x.trial(in,ch,:),1),2));
    end
    used = find(~isnan(PPi));
      
    SSVEPtrace(pp,1,:) = nanmean(SSVEPtr(PPi(used) < prctile(PPi(used),50),:));
    SSVEPtrace(pp,2,:) = nanmean(SSVEPtr(PPi(used) > prctile(PPi(used),50),:));
    TimeSSVEP{pp}(1,:) = nanmean(TimeSSVEPtr(PPi(used) < prctile(PPi(used),50),:));
    TimeSSVEP{pp}(2,:) = nanmean(TimeSSVEPtr(PPi(used) > prctile(PPi(used),50),:));      
end

%% show the SSVEP, FFT wise
close all
addpath('D:\Other\Programs\Matlab\extrausefullScripts');
foi = [2:0.5:15];
clear allpp
for pp = 1:length(ppnames)
    subplot(6,7,pp)
    plot(foi, squeeze(SSVEPtrace(pp,:,:)));
    title(num2str(freq(pp)));
    x = nearest(foi,freq(pp));
    if x <= 6
        allpp(pp,:,:)= NaN;
        allpp(pp,:,[x-x+1:x+6]+x-x+1) = SSVEPtrace(pp,:,[x-x+1:x+6]);
    else
        allpp(pp,:,:) = SSVEPtrace(pp,:,[x-6:x+6]);
    end
end
figure
subplot(2,1,1)
tosqrt = sum(~isnan(squeeze(mean(allpp,2))))';
h = shadedErrorBar([-6:6]./2, squeeze(nanmean(allpp(:,1,:),1)), squeeze(nanstd(allpp(:,1,:),1))./sqrt(tosqrt));
hold on
h2 = shadedErrorBar([-6:6]./2, squeeze(nanmean(allpp(:,2,:),1)), squeeze(nanstd(allpp(:,2,:),1))./sqrt(tosqrt), 'r');
legend([h.mainLine h2.mainLine], {'low';'high'});
subplot(2,1,2)
plot(foi, squeeze(SSVEPtrace(32,:,:)));
set(gcf, 'position',[500 500 375 270]);

%% show the SSVEP, trace itself
figure
clear SSVEP
for pp = 1:length(ppnames)
    subplot(6,7,pp)
    tp = find([-2:1/256:0] >= -2);
    SSVEP{pp}(1,:) = TimeSSVEP{pp}(1,tp)-mean(TimeSSVEP{pp}(1,tp));
    SSVEP{pp}(2,:) = TimeSSVEP{pp}(2,tp)-mean(TimeSSVEP{pp}(2,tp));
    plot([-2:1/256:0], squeeze(SSVEP{pp})');    
    set(gca,'xlim', [-1 0]);
end

legend({'low';'high'});

%% show the SSVEP align based on stimulus... kinda moving window
close all
clear tc
for pp = 1:length(ppnames)
    tp = find([-2:1/256:0] >= -2);
    SSVEP{pp}(1,:) = TimeSSVEP{pp}(1,tp)-mean(TimeSSVEP{pp}(1,tp));
    SSVEP{pp}(2,:) = TimeSSVEP{pp}(2,tp)-mean(TimeSSVEP{pp}(2,tp));
    tc(pp,:) = linspace(-freq(pp)*2,1, size(tp,2));
end
[m inMa] = max(freq);
[m inMi] = min(freq);

mvv = zeros(40,2,length(tc)-1);
for pp = 1:length(ppnames)       
    for tpc = 1:length(tc)-1
        in = find(tc(pp,:) >= tc(inMa,tpc) & tc(pp,:) < tc(inMa,tpc+1));
        mvv(pp,:,tpc) = nanmean(SSVEP{pp}(:,in),2);
    end
end
x = squeeze(mvv(:,1,:));

%% make the figure
p = 31;
figure
subplot(2,2,1)
plot(tc(inMa,1:length(tc)-1)+0.16, squeeze(mean(mvv,1)));
set(gca,'xlim',[-4 0]);
subplot(2,2,3)
plot([-2:1/256:0]+0.16, squeeze(SSVEP{p})');
set(gca,'xlim',[-0.75 0]);

subplot(2,4,3)
tosqrt = sum(~isnan(squeeze(mean(allpp,2))))';
h = shadedErrorBar([-6:6]./2, squeeze(nanmean(allpp(:,1,:),1)), squeeze(nanstd(allpp(:,1,:),1))./sqrt(tosqrt));
hold on
h2 = shadedErrorBar([-6:6]./2, squeeze(nanmean(allpp(:,2,:),1)), squeeze(nanstd(allpp(:,2,:),1))./sqrt(tosqrt), 'r');
set(gca,'xlim',[-3 3]);
%legend([h.mainLine h2.mainLine], {'low';'high'});
subplot(2,4,7)
plot(foi, squeeze(SSVEPtrace(p,:,:)));
set(gca, 'xlim', [2 15]);
set(gcf, 'position',[500 500 400 300]);
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S03_01_FFT_' infor.pr{pr}], 'Color', 'rgb');
print(['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S03_01_FFT_' infor.pr{pr}],'-dpng');

