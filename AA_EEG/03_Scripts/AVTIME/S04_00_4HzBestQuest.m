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

%% extract behavioral response:
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_all'], 'AV');
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_1'], 'FF');
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Behdata.mat'], 'ppnames', 'vars');
ppto = setdiff(1:40,ol);
clear proplp amtr
allfr = {[1:100], [102:201]};
allfrab17 = {[1:201], [203:400]};
tp = [50 175];
bl = [0 1 2 3 4];
for pp = 1:length(ppto)
    for t = 1:2
        for fr = 1:2
            inx = find(vars(ppto(pp)).SOUND2 == fr & vars(ppto(pp)).TIMEPOINT == tp(t)& ismember(vars(ppto(pp)).BLOCK, bl));
            amlp = find(vars(ppto(pp)).RESPONSE(inx) == AN(AS(ppto(pp)),fr));
            proplp(pp,fr,t) = length(amlp)./length(inx);
            amtr(pp,fr,t) = length(inx);
        end
        %proplp(pp,fr,t) = proplp(pp,fr,t)./mean(proplp(pp,:,t));
    end
end
dif = proplp(:,:,2)-proplp(:,:,1);
difN = dif(:,2)-dif(:,1); % association effect

%% plot association effect with frequency and frequency difference
addpath(genpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231'));

close all
t2{1} = abs(freq-4);
t2{2} = freq;
l{1} = 'absolute difference (Hz)';
l{2} = 'frequency (Hz)';
figure
for it = 1:2
    subplot(1,2,it);
    scatter(difN, t2{it}(ppto), 200,[0 0 0], '.');    
    [r p] = corr(difN, t2{it}(ppto)');
    b = polyfit(difN,t2{it}(ppto)',1);
    hold on
    plot([-0.5 0.6], b(1).*[-0.5 0.6]+b(2), 'r');
    xlabel('Association effect');ylabel(l{it});
    set(gca, 'xlim', [-0.2 0.2]); 
    y = get(gca,'ylim'); x = get(gca,'xlim');
    text(x(1)+0.01, y(2)-0.03*diff(y), ['r=' num2str(round(r, 2))]);
    text(x(1)+0.01, y(2)-0.1*diff(y), ['p=' num2str(round(p, 3))]);       
    set(gca, 'ylim', y);
end
set(gcf, 'position', [300 300 500 200]);
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S04_00_4HzBest'], 'Color', 'rgb');
print('D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S04_00_4HzBest','-dpng')
