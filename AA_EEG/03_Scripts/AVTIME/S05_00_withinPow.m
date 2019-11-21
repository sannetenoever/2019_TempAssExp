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
pra = {-7:-1, 0};
infor.pr = {'PowerPreTrials';'PowerCurrentTrial'};
infor.asna = {'nonassociated';'associated'};
infor.resp = {'incorrect';'correct'};
for pp = 1:length(ppnames)
    load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_' ppnames{pp}], 'FF');
    frinx = nearest(FF.freq, freq(pp));
    it=1;
    for pr = 1:2
        for asna = 1:2
            for resp = 1:2
                trl = find(FF.trialinfo(:,4) == asna-1  & FF.trialinfo(:,7)+1 == resp & FF.trialinfo(:,9) >= 0);
                PPi = [];
                for trc = 1:length(trl)
                    in = find(ismember(FF.trialinfo(:,11)-FF.trialinfo(trl(trc),11),pra{pr}) & FF.trialinfo(:,9) >= 0);
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
pr = 1;
V1.individual = squeeze(PP(:,pr,1,2,:)-PP(:,pr,1,1,:));
V1.time = [0];
V1.dimord = 'subj_chan_time';
load([loc '\eyeblinkcorrectionRMT_AVTIME1'], 'datAud_preproc_RM');
datAud_preproc_RM.label{16}= 'Fp1';
datAud_preproc_RM.label{17}= 'Fp2';
V1.label = datAud_preproc_RM.label;
V2 = V1;
V2.individual = squeeze(PP(:,pr,2,2,:)-PP(:,pr,2,1,:));

cfg = [];
cfg.latency = [-1 1];
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
%save(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\S05_00_PowAccInterStat' infor.pr{pr}],'st');

%% plot the interaction:
pr=1;
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\S05_00_PowAccInterStat' infor.pr{pr}],'st');
load('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231\template\layout\easycapM1');
cfg = [];
cfg.layout = lay;
cfg.parameter = 'stat';
cfg.highlightchannel = {st.label{st.posclusterslabelmat==1}};
cfg.highlight = 'on';
ft_topoplotER(cfg,st);
colorbar
set(gcf, 'position',[500 500 375 270]);
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_00_PowAccInterTopo_' infor.pr{pr}], 'Color', 'rgb');
print(['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_00_PowAccInterTopo_' infor.pr{pr}],'-dpng');

%% repeat separately for the 30% time point and 70% time point effect:
pr = 1;
V1.individual = squeeze(PP(:,pr,1,1,:)); % incor
V1.time = [0];
V1.dimord = 'subj_chan_time';
load([loc '\eyeblinkcorrectionRMT_AVTIME1'], 'datAud_preproc_RM');
datAud_preproc_RM.label{16}= 'Fp1';
datAud_preproc_RM.label{17}= 'Fp2';
V1.label = datAud_preproc_RM.label;
V2 = V1;
V2.individual = squeeze(PP(:,pr,1,2,:)); % cor

cfg = [];
cfg.latency = [-1 1];
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
cfg.alpha = 0.05;
cfg.numrandomization = 1000;
cfg.design(1,1:2*length(ppnames))  = [ones(1,length(ppnames)) 2*ones(1,length(ppnames))];
cfg.design(2,1:2*length(ppnames))  = [1:length(ppnames) 1:length(ppnames)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

st30 = ft_timelockstatistics(cfg, V1,V2);
V1.individual = squeeze(PP(:,pr,2,1,:)); % incor
V2.individual = squeeze(PP(:,pr,2,2,:)); % cor
st70 = ft_timelockstatistics(cfg, V1,V2);
save(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\S05_00_PowAccSimpleStat' infor.pr{pr}],'st30','st70');

%% look at simple effects
pr=1;
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\S05_00_PowAccInterStat' infor.pr{pr}],'st70', 'st30');
load('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231\template\layout\easycapM1');
cfg = [];
cfg.layout = lay;
cfg.parameter = 'stat';
cfg.highlightchannel = {st70.label{st70.posclusterslabelmat==1}};
cfg.highlight = 'on';
ft_topoplotER(cfg,st70);
colorbar
set(gcf, 'position',[500 500 375 270]);
%exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_00_PowAccInterTopo_' infor.pr{pr}], 'Color', 'rgb');
%print(['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_00_PowAccInterTopo_' infor.pr{pr}],'-dpng');

%% look at the power
figure
ppto = 1:40;
adv = 0.5; % just added for the plot
ch = find(st.posclusterslabelmat==1);
voi = squeeze(mean(PP(:,pr,:,:,ch),5))+adv;
bar(squeeze(mean(voi,1))); hold on
% make within subject error bar for nonass and ass seperately
for asna = 1:2
    t = squeeze(voi(:,asna,:));
    temp = squeeze(t)-repmat(mean(t,2),[1 2])-mean(t(:));
    mb(asna,:) = squeeze(std(temp(ppto,:))./sqrt(length(ppto)));
end
errorbar([0.86 1.14; 1.86 2.14],squeeze(mean(voi,1)),mb, 'linestyle', 'none');
set(gca, 'xtick', [1 2], 'xlim', [0.5 2.5], 'xticklabel', {'non-associated';'associated'});
x = get(gca, 'yticklabel');
set(gca,'yticklabel', cellfun(@(x) num2str(str2double(x)-adv),x,'uniformoutput',0));
legend({'incor','cor'},'location','southeast');
ylabel('log power');
%[h p ci t] = ttest(squeeze(voi(:,1,1)), squeeze(voi(:,1,2)))
%[h p ci t] = ttest(squeeze(voi(:,2,1)), squeeze(voi(:,2,2)))
set(gcf, 'position',[500 500 375 270]);
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_00_SimplePowAcc_' infor.pr{pr}], 'Color', 'rgb');
print(['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_00_SimplePowAcc_' infor.pr{pr}],'-dpng')

%% big matrix for R
saveloc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\';
load([saveloc 'S01_02_outlierdef'], 'ol');
varnames = {'accuracy';'sound';'ass';'block';'pp';'curval'; 'dev'; 'pow'};
bigmat = [];
for pp = 1:length(ppnames)
    load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_' ppnames{pp}], 'FF');
    
    % get pretrial power:
    frinx = nearest(FF.freq, freq(pp));
    PPi = [];
    for trc = 1:length(FF.trialinfo)
        in = find(ismember(FF.trialinfo(:,11)-FF.trialinfo(trc,11),pra{pr}) & FF.trialinfo(:,9) >= 0);        
        PPi(trc) = squeeze(mean(mean(log(FF.pow(in,ch,frinx)),1),2));        
    end
    used = find(~isnan(PPi));   
    
    inxT = size(bigmat,2) + 1:size(bigmat,2) + length(used);  
    bigmat(1,inxT) = FF.trialinfo(used,7); % accuracy
    bigmat(2,inxT) = FF.trialinfo(used,3); % sound
    bigmat(3,inxT) = FF.trialinfo(used,4); % association
    bigmat(4,inxT) = FF.trialinfo(used,9) >= 3; % block
    bigmat(5,inxT) = pp; % pp
    bigmat(6,inxT) = FF.trialinfo(used,5); % currentvalue
    bigmat(7,inxT) = mat2gray(FF.trialinfo(used,11)); % time development
    bigmat(8,inxT) = zscore(PPi(used));     
end
bigmat(:,ismember(bigmat(5,:),ol)) = [];
bigmat(6,:) = zscore(bigmat(6,:)); % zscore difficulty
bigmat = bigmat';
save([saveloc 'S05_00_BIGMAT'], 'bigmat');
save([saveloc 'S05_00_BIGMAT_varnames'], 'varnames');

%% now look at behavior only for the low power trials:
addpath(genpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231'));
tp = [50 175];
pr=1;
for pp = 1:length(ppnames)
    load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_' ppnames{pp}], 'FF');
    
    % only extract lowest 50% pretrial power:
    frinx = nearest(FF.freq, freq(pp));
    PPi = [];
    for trc = 1:length(FF.trialinfo)
        in = find(ismember(FF.trialinfo(:,11)-FF.trialinfo(trc,11),pra{pr}) & FF.trialinfo(:,9) >= 0);        
        PPi(trc) = squeeze(mean(mean(log(FF.pow(in,ch,frinx)),1),2));        
    end
    used = find(~isnan(PPi));
    trialinfo = FF.trialinfo(used,:);
    PPi = PPi(used);
    inx = find(PPi < prctile(PPi,50));
    
    for it = 1:2
        if it == 1
            trlinf = trialinfo(inx,:);
        else
            trlinf = trialinfo(setdiff(1:length(trialinfo),inx),:);
        end
        
        blc = {[0 1]; [2 3 4]};
        for bl = 1:size(blc,1)
            for fr = 1:2
                for t = 1:2
                    inx = find(trlinf(:,3) == fr & trlinf(:,4)== t-1 & ismember(trlinf(:,9), blc{bl}));            
                
                    %inx = find(trlinf(:,3)==fr & trlinf(:,2) == tp(t)& ismember(trlinf(:,9), blc{bl}));
                    amlp = sum(trlinf(inx,7));
                    allv.proplp(pp,it,bl,fr,t) = amlp./length(inx);
                    allv.pp(pp,it,bl,fr,t) = pp;
                    allv.bl(pp,it,bl,fr,t) = bl;
                    allv.t(pp,it,bl,fr,t) = t;
                    allv.s(pp,it,bl,fr,t) = fr;
                    allv.pow(pp,it,bl,fr,t) = it;
                    amtr(pp,it,bl,fr,t) = length(inx);
                end
            end
        end
    end
end

%%
% check for outliers
v = squeeze(allv.proplp(:,2,:,2)-allv.proplp(:,2,:,1));
d = mean(v,2);
d = v(:,2);
ol = [find(d > mean(d) + 2.5*std(d)) find(d < mean(d) - 2.5*std(d))];
d = v(:,1);
ol = [ol find(d > mean(d) + 2.5*std(d)) find(d < mean(d) - 2.5*std(d))]
ppto = setdiff(1:length(ppnames),ol);

figure
it=1;
leg = {'30%', '70%'};
tl = {'first part';'second part'};
for bl = 1:2
    subplot(1,2,bl)
    bar(squeeze(mean(allv.proplp(ppto,it,bl,:,:),1)));
    hold on
    errorbar([0.85 1.15; 1.85 2.15], squeeze(mean(allv.proplp(ppto,it,bl,:,:),1)), squeeze(std(allv.proplp(ppto,it,bl,:,:),[],1))./sqrt(length(ppto)), 'k','linestyle','none', 'linewidth',3)
    set(gca, 'ylim', [0.5 0.65], 'xticklabel', {'SOUND A';'SOUND B'})
    if bl == 1
        legend(leg, 'location', 'southeast');
    end
    ylabel('proportion correct');
    title(tl{bl})
    set(gcf, 'position', [500 200 600 250]);
end
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_00_BehLowPower_' infor.pr{pr}], 'Color', 'rgb');
print(['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S05_00_BehLowPower_' infor.pr{pr}],'-dpng')

%% do the RM anova
it = 1;
restoredefaultpath();
addpath('D:\Other\Programs\Matlab\extrausefullScripts\');
ft.propcor = allv.proplp(ppto,it,:,:,:);
ft.pp = allv.pp(ppto,it,:,:,:);
ft.block = allv.bl(ppto,it,:,:,:);
ft.associated = allv.t(ppto,it,:,:,:);
ft.sound = allv.s(ppto,it,:,:,:);

tb = table(ft.pp(:),ft.block(:),ft.associated(:),ft.sound(:), ft.propcor(:), 'VariableNames', {'pp', 'BL','ASS','S','pCOR'});
tb = uni2multiTB(tb, [2:4], 5, 1);
F1 = [1 1 1 1 2 2 2 2]';
F2 = [1 1 2 2 1 1 2 2]';
F3 = [1 2 1 2 1 2 1 2]';
withdes = table(F1,F2,F3, 'VariableNames', {'BLOCK';'ASS';'SOUND'});
withdes.BLOCK = categorical(withdes.BLOCK);
withdes.ASS = categorical(withdes.ASS);
withdes.SOUND = categorical(withdes.SOUND);

rm = fitrm(tb,'BL1ASS1S1-BL2ASS2S2~1','WithinDesign',withdes);
ranovatbl = ranova(rm, 'WithinModel','BLOCK*ASS*SOUND');

%% just to ensure that both seperate are significant (although no
% interaction so strickly not necessary)
F1 = [1 1 2 2]';
F2 = [1 2 1 2]';
withdes = table(F1, F2, 'VariableNames', {'BLOCK';'ASS'});
withdes.ASS = categorical(withdes.ASS);

