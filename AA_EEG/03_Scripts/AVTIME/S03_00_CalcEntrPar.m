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

%%
foi = [2:0.5:15];
bl = [0:5];
for pp = 1:length(ppnames)
    load([loc '\trialinfo_' expname ppnames{pp}]);
    load([loc '\eyeblinkcorrectionRMT_' expname ppnames{pp}], 'datAud_preproc_RM');
    trlinf = trlinf(setdiff(trlinf(:,end),datAud_preproc_RM.badtr),:);
   
    % align to visual
    %offset = trlinf(:,10);
    %in = find(ismember(trlinf(:,9), bl));
    %cfg = [];
    %cfg.trials = in;
    %cfg.offset = offset(in);
    %datn = ft_redefinetrial(cfg, datAud_preproc_RM);
%     clear FF    
%     cfg = [];
%     cfg.toilim = [-2 2];
%     dat = ft_redefinetrial(cfg, datn);
%     for it = 1:length(dat.trial)
%         x = nearest(dat.time{it}, 0);
%         l = length(dat.time{it})-x;
%         dat.trial{it}(:,x:end) = repmat(dat.trial{it}(:,x), [1 l+1]);
%     end
%     dat.trialinfo = trlinf;

    % align to auditory    
    in = find(ismember(trlinf(:,9), bl));
    if str2double(ppnames{pp}) < 42
        cfg = [];
        cfg.trials = in;
        cfg.offset = round(-0.110.*datAud_preproc_RM.fsample);
        datn = ft_redefinetrial(cfg, datAud_preproc_RM);
    else % different due to computer update
        cfg = [];
        cfg.trials = in;
        cfg.offset = 0;
        datn = ft_redefinetrial(cfg, datAud_preproc_RM);
    end
    datn.trialinfo = trlinf;
    
    load('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231\template\layout\easycapM1.mat');
    for fr = 1:length(foi)
        cfg = [];
        cfg.toilim = [-4/foi(fr) 0];
        datt = ft_redefinetrial(cfg, datn);
        
        cfg =[];
        cfg.foi = foi(fr);
        cfg.method = 'mtmfft';
        cfg.taper = 'hanning';
        cfg.output = 'fourier';
        frt = ft_freqanalysis(cfg, datt);
        if fr == 1
           fra = frt;
        else
           fra.fourierspctrm(:,:,fr) = frt.fourierspctrm;
           fra.freq(fr) = frt.freq;
        end
    end
    FF.fourierspctrm = fra.fourierspctrm;
    FF.freq = fra.freq;
    FF.label = fra.label;
    FF.dimord = fra.dimord;
    FF.trialinfo = datn.trialinfo;
    FF.angle = angle(FF.fourierspctrm);
    FF.pow = abs(FF.fourierspctrm).^2;
    
    AVi.savemrvla = circ_r(FF.angle);
    AVi.saveangvl = circ_mean(FF.angle);
    AVi.savepowvl = mean(log(FF.pow));
    
    AV.savemrvla(pp,:,:) = circ_r(FF.angle);
    AV.saveangvl(pp,:,:) = circ_mean(FF.angle);
    AV.savepowvl(pp,:,:) = mean(log(FF.pow));

    save(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_' ppnames{pp}], 'FF','fra', 'AVi');
end

save(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\entrresp_all'], 'AV');


