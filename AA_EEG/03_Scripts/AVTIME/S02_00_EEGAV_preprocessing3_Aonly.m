clear
ppn = 41

expname = 'AVTIME';

addpath(genpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231'));
addpath('D:\Experiments\StartTimeExp\AA_EEG\03_Scripts\AVTIME');
loclog = 'D:\Experiments\StartTimeExp\AA_EEG\02_Data\Logs\';
addpath('D:\Other\Programs\Matlab\extrausefullScripts');

ppnames = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'41';'42';'43';};
pp = ppnames{ppn};

namefile = [pp '_' expname '.vhdr'];
loc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\PreProcData\';

% preprocessing for one file
cd D:\Experiments\StartTimeExp\AA_EEG\

%%
% first load the layout
load(['D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231\template\layout\easycapM1']);
%load([loc '\eyeblinkcorrection_' expname pp], 'datAud_preproc');
load([loc '\firstdat_' expname pp], 'datA', 'trlAud', 'trlVis');
load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Behdata.mat'], 'vars');
load([loc '\eventdata_' expname pp], 'eventdata');

if pp == '34' | pp == '37'
    x = load([loc '\eventdata_' expname '_part2' pp], 'eventdata');
    y = load([loc '\firstdat_' expname '_part2' pp], 'datA', 'trlAud', 'trlVis');
end

cfg = [];
cfg.detrend = 'yes';
cfg.resamplefs = 256;
dat_rmbchan = ft_resampledata(cfg, datA);
dat_rmbchan.sampleinfo = round(datA.sampleinfo/(2500/256));
if diff(dat_rmbchan.sampleinfo) ~= length(dat_rmbchan.trial{1})-1
    dat_rmbchan.sampleinfo = [dat_rmbchan.sampleinfo(1)+1 dat_rmbchan.sampleinfo(2)];
end
trlAud(:,1:3) = round(trlAud(:,1:3)/(2500/256));
trlVis(:,1:3) = round(trlVis(:,1:3)/(2500/256));
cfg = [];
cfg.trl = trlAud(trlAud(:,4)<15,:);
datAud_preproc = ft_redefinetrial(cfg, dat_rmbchan);

if pp == '34' | pp == '37'
    cfg = [];
    cfg.detrend = 'yes';
    cfg.resamplefs = 256;
    dat_rmbchan2 = ft_resampledata(cfg, y.datA);
    dat_rmbchan2.sampleinfo = round(y.datA.sampleinfo/(2500/256));
    if diff(dat_rmbchan2.sampleinfo) ~= length(dat_rmbchan2.trial{1})-1
        dat_rmbchan2.sampleinfo = [dat_rmbchan2.sampleinfo(1)+1 dat_rmbchan2.sampleinfo(2)];
    end
    y.trlAud(:,1:3) = round(y.trlAud(:,1:3)/(2500/256));
    y.trlVis(:,1:3) = round(y.trlVis(:,1:3)/(2500/256));
    cfg = [];
    cfg.trl =y.trlAud(y.trlAud(:,4)<15,:);
    datAud_preproc2 = ft_redefinetrial(cfg, dat_rmbchan2);
    inx2 = find(cellfun(@(x) strcmp(x, 'S 15'), {x.eventdata.value}));
end

inx = find(cellfun(@(x) strcmp(x, 'S 15'), {eventdata.value}));
%{eventdata(inx+1).value}'
if pp == '12'
    cfg = [];
    cfg.trials = [2:501];
    datAud_preproc = ft_redefinetrial(cfg, datAud_preproc);
end
if pp == '7'
     datAud_preproc.trialinfo = [datAud_preproc.trialinfo'; vars(str2double(pp)).TIMEPOINT([1:455 457:500]); vars(str2double(pp)).SOUND2([1:455 457:500]); vars(str2double(pp)).ASS([1:455 457:500]); vars(str2double(pp)).currentvalue([1:455 457:500]); vars(str2double(pp)).RESPONSE([1:455 457:500]); vars(str2double(pp)).RESC([1:455 457:500]); vars(str2double(pp)).RT([1:455 457:500]); vars(str2double(pp)).BLOCK([1:455 457:500])]';
elseif pp == '16'
    datAud_preproc.trialinfo = [datAud_preproc.trialinfo'; vars(str2double(pp)).TIMEPOINT([1:283 285:500]); vars(str2double(pp)).SOUND2([1:283 285:500]); vars(str2double(pp)).ASS([1:283 285:500]); vars(str2double(pp)).currentvalue([1:283 285:500]); vars(str2double(pp)).RESPONSE([1:283 285:500]); vars(str2double(pp)).RESC([1:283 285:500]); vars(str2double(pp)).RT([1:283 285:500]); vars(str2double(pp)).BLOCK([1:283 285:500])]';
elseif strcmp(pp, '31')
    datAud_preproc.trialinfo = [datAud_preproc.trialinfo'; vars(str2double(pp)).TIMEPOINT([2:500]); vars(str2double(pp)).SOUND2([2:500]); vars(str2double(pp)).ASS([2:500]); vars(str2double(pp)).currentvalue([2:500]); vars(str2double(pp)).RESPONSE([2:500]); vars(str2double(pp)).RESC([2:500]); vars(str2double(pp)).RT([2:500]); vars(str2double(pp)).BLOCK([2:500])]';
elseif strcmp(pp, '34')
    datAud_preproc.trialinfo = [datAud_preproc.trialinfo'; vars(str2double(pp)).TIMEPOINT([1:200]); vars(str2double(pp)).SOUND2([1:200]); vars(str2double(pp)).ASS([1:200]); vars(str2double(pp)).currentvalue([1:200]); vars(str2double(pp)).RESPONSE([1:200]); vars(str2double(pp)).RESC([1:200]); vars(str2double(pp)).RT([1:200]); vars(str2double(pp)).BLOCK([1:200])]';
    datAud_preproc2.trialinfo = [datAud_preproc2.trialinfo'; vars(str2double(pp)).TIMEPOINT([201:500]); vars(str2double(pp)).SOUND2([201:500]); vars(str2double(pp)).ASS([201:500]); vars(str2double(pp)).currentvalue([201:500]); vars(str2double(pp)).RESPONSE([201:500]); vars(str2double(pp)).RESC([201:500]); vars(str2double(pp)).RT([201:500]); vars(str2double(pp)).BLOCK([201:500])]';
elseif strcmp(pp, '37')
    in1 = [1:413];
    in2 = [415:499];
    datAud_preproc.trialinfo = [datAud_preproc.trialinfo'; vars(str2double(pp)).TIMEPOINT([in1]); vars(str2double(pp)).SOUND2([in1]); vars(str2double(pp)).ASS([in1]); vars(str2double(pp)).currentvalue([in1]); vars(str2double(pp)).RESPONSE([in1]); vars(str2double(pp)).RESC([in1]); vars(str2double(pp)).RT([in1]); vars(str2double(pp)).BLOCK([in1])]';
    datAud_preproc2.trialinfo = [datAud_preproc2.trialinfo'; vars(str2double(pp)).TIMEPOINT([in2]); vars(str2double(pp)).SOUND2([in2]); vars(str2double(pp)).ASS([in2]); vars(str2double(pp)).currentvalue([in2]); vars(str2double(pp)).RESPONSE([in2]); vars(str2double(pp)).RESC([in2]); vars(str2double(pp)).RT([in2]); vars(str2double(pp)).BLOCK([in2])]';
elseif strcmp(pp, '43')
    in1 = [1:500];
    datAud_preproc.trialinfo = [datAud_preproc.trialinfo'; vars(ppn).TIMEPOINT([in1]); vars(ppn).SOUND2([in1]); vars(ppn).ASS([in1]); vars(ppn).currentvalue([in1]); vars(ppn).RESPONSE([in1]); vars(ppn).RESC([in1]); vars(ppn).RT([in1]); vars(ppn).BLOCK([in1])]';   
else    
    in1 = [1:500];
    datAud_preproc.trialinfo = [datAud_preproc.trialinfo'; vars(ppn).TIMEPOINT([in1]); vars(ppn).SOUND2([in1]); vars(ppn).ASS([in1]); vars(ppn).currentvalue([in1]); vars(ppn).RESPONSE([in1]); vars(ppn).RESC([in1]); vars(ppn).RT([in1]); vars(ppn).BLOCK([in1])]';   
end

% check trial length (due to resampling problem now)
[yz l] = min(cellfun(@(x) length(x), datAud_preproc.time));
datAud_preproc.trial = cellfun(@(x) x(:,1:yz), datAud_preproc.trial,'uniformoutput',0);
datAud_preproc.time = cellfun(@(x) datAud_preproc.time{l}, datAud_preproc.trial,'uniformoutput',0);
if pp == '34' | pp == '37'
    datAud_preproc2.trial = cellfun(@(x) x(:,1:yz), datAud_preproc2.trial,'uniformoutput',0);
    datAud_preproc2.time = cellfun(@(x) datAud_preproc.time{l}, datAud_preproc2.trial,'uniformoutput',0);
    datAud_preproc = ft_appenddata([],datAud_preproc, datAud_preproc2);
end

cfg = [];
cfg.keeptrials = 'yes';
tempd = ft_timelockanalysis(cfg, datAud_preproc);
trlinf = datAud_preproc.trialinfo;
if pp == '6' | pp == '10' | ppn == 15 | ppn == 24 | pp == 33
    offset = [17; (trlAud(2:500,1)-trlAud(2:500,3))-(trlVis(:,1)-trlVis(:,3))]; 
elseif pp == '7' 
    offset = [(trlAud(:,1)-trlAud(:,3))-(trlVis([1:455 457:500],1)-trlVis([1:455 457:500],3))];    
elseif ppn == 16
    offset = [(trlAud(:,1)-trlAud(:,3))-(trlVis(1:499,1)-trlVis(1:499,3))];
elseif ppn == 16
    offset = [(trlAud(:,1)-trlAud(:,3))-(trlVis([1:283 285:500],1)-trlVis([1:283 285:500],3))];
elseif ppn == 12   
     offset = [(trlAud(2:501,1)-trlAud(2:501,3))-(trlVis(2:501,1)-trlVis(2:501,3))]; 
elseif ppn == 31
    offset = [17; (trlAud(2:499,1)-trlAud(2:499,3))-(trlVis(:,1)-trlVis(:,3))]; 
elseif ppn == 34
    offset = [48; (trlAud(2:200,1)-trlAud(2:200,3))-(trlVis(:,1)-trlVis(:,3)); (y.trlAud(:,1)-y.trlAud(:,3))-(y.trlVis(:,1)-y.trlVis(:,3))];
elseif ppn == 37
     offset = [(trlAud(:,1)-trlAud(:,3))-(trlVis(1:end-1,1)-trlVis(1:end-1,3)); (y.trlAud(:,1)-y.trlAud(:,3))-(y.trlVis(1:end-1,1)-y.trlVis(1:end-1,3))];
else
    offset = (trlAud(:,1)-trlAud(:,3))-(trlVis(:,1)-trlVis(:,3));
end
trlinf(:,end+1) = offset;
trlinf(:,end+1) = 1:size(trlinf,1);

save([loc '\trialinfo_' expname pp], 'trlinf');

%% do via EEGlab because of easier vis
%eeglab
addpath(genpath('D:\Other\Programs\Matlab\eeglab14_1_1b'));
TMPREJ=[];
eegplot(permute(tempd.trial,[2 3 1]),'winlength',3,'command','pippo');

%% additional channels to remove?
badchn = {'100'};
tmp = 1;
while ~isempty(badchn)
    badch(tmp) = sscanf(badchn{1}, '%d');
    badchn = inputdlg('enter bad channels, when no more, press cancel');
    tmp = tmp+1;
end
badch(1) = [];
close all

%%
clear trialrej trialrejV
if ~isempty(TMPREJ)
    [trialrej elecrej]=eegplot2trial(TMPREJ,size(tempd.trial,3),size(tempd.trial,1));
else
    trialrej=[];
end

restoredefaultpath();
addpath(genpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231'));

cfg = [];
cfgn = [];
cfg.badchannel = datAud_preproc.label(badch);
cfgn.layout = lay;
cfgn.method = 'distance';
cfgn.neighbourdist = 0.25;
cfg.neighbours = ft_prepare_neighbours(cfgn, datAud_preproc);
cfg.layout = lay;
datAud_preproc = ft_channelrepair(cfg, datAud_preproc);

datAud_preproc.trialinfo(:,8) = 1:size(datAud_preproc.trialinfo,1);

cfg = [];
cfg.trials = find(~trialrej);
datAud_preproc_RM = ft_redefinetrial(cfg, datAud_preproc);

%% Before applying ICA to reduce artifacts, check which is the maximum number of independent components via singular value decomposition
cfg = [];
cfg.keeptrials = 'yes';
tempdz = ft_timelockanalysis(cfg, datAud_preproc_RM);

tempd = permute(tempdz.trial, [2 3 1]);
tempd = tempd(1:31,:,:);

tmpdata=reshape(tempd,[size(tempd,1) size(tempd,2)*size(tempd,3)]);
tmpdata=tmpdata-repmat(mean(tmpdata,2),[1,size(tempd,2)*size(tempd,3)]);
tempCpca=svd(tmpdata);
th=0.1;% this threshold is arbitrary but it works in most cases
figure,semilogy(tempCpca,'.-')% plot of the eigenvalues of the data matrix
a=ylim;
hold on
Cpca=max(find(tempCpca>th))
line([Cpca Cpca],[a(1) a(2)],'Color','k','LineStyle','--');
hold off
% if the vertical dotted line correctly separates the eigenvalues close to
% zero from the ones much higher than zero, Cpca correctly corresponds to
% the total number of independent components that can be estimated from the
% data matrix. if this is not the case, you have to change the threshold above or to manually correct Cpca

% move into EEGlab structu
addpath(genpath('D:\Other\Programs\Matlab\eeglab14_1_1b'));
PathName = 'D:\Experiments\StartTimeExp\AA_EEG\02_Data\EEG\';
namefile = ['1_' expname '.vhdr'];
[EEG s] = pop_loadbv(PathName,namefile);
EEGn.chanlocs = EEG.chanlocs(1:31);
l = readlocs('D:\Other\Programs\Matlab\eeglab14_1_1b\sample_locs\Standard-10-20-Cap81.locs');
[o inxa inxb] = intersect(upper({EEGn.chanlocs.labels}), upper({l.labels}));
EEGn.chanlocs = [l(inxb)];
amch = length(o);
EEGn.data=tempd(inxa,:,:);
EEGn.nbchan=amch;
EEGn.times = tempdz.time;
EEGn.nbchan = size(EEGn.data,1);
EEGn.pnts = size(EEGn.data,2);
EEGn.trials = size(EEGn.data,3);
EEGn.srate = 256;

%% ICA
addpath('D:\Other\Programs\Matlab\extrausefullScripts');
EEGn=ICA_analysis(EEGn,EEGn.data,Cpca);% run ICA using the runica method from EEGlab

EEGpart.ICA=EEGn;
EEGn.comp2remove=[];

%%
ICA_selectandremove(EEGn);

%%
TMPREJ=[];
eegplot(data_GUI.compproj,'winlength',3,'command','pippo');

%%
datAud_preproc_RM = datAud_preproc_RM;
datAud_preproc_RM.trial = arrayfun(@(x) data_GUI.compproj(:,:,x), [1:size(data_GUI.compproj,3)], 'uniformoutput', 0);
datAud_preproc_RM.badtr = find(trialrej);
%datAud_preproc_RM.time = datAud_preproc.time(~trialrej);
%datAud_preproc_RM.trialinfo = datAud_preproc.trialinfo(~trialrej);
datAud_preproc_RM.comprej = data_GUI.comp2remove;
datAud_preproc_RM.label = datAud_preproc.label(inxa);
preICA = tempdz;

save([loc '\eyeblinkcorrectionRMT_' expname pp], 'datAud_preproc_RM', 'preICA');




