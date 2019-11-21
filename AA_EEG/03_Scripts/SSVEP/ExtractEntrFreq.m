% Information section
clear
clc
root_dir = 'D:\Experiments\StartTimeExp\AA_EEG\02_Data\EEG\'; % the location of the EEG data
loclog = 'D:\Experiments\StartTimeExp\AA_EEG\02_Data\Logs\'; % the location of the EEG data
locfieldtrip = 'D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231\'; % the location of the fieldtrip files
figloc = 'D:\Experiments\StartTimeExp\AA_EEG\05_Figures\SSVEP\';

addpath('D:\Experiments\StartTimeExp\AA_EEG\03_Scripts\SSVEP');
addpath(genpath(locfieldtrip));
addpath('D:\Other\Programs\Matlab\CircStat2012a\');
addpath('D:\Other\Programs\Matlab\extrausefullScripts\');    

cd(root_dir)

%%
for pp = [1:37 42:43]
    close all   
    if pp == 7 
        partname = [num2str(pp) '_SSVEP'];
    else
        partname = [num2str(pp) '_SSVEP']; % the name of the participant 
        lf = dir(loclog);
        lf(find(cellfun(@(x) length(x) < 7, {lf.name}))) =[];
        if pp < 10
            dloc = find(cellfun(@(x) strcmp(x(1:7), [num2str(pp) '_SSVEP']) && contains(x, ['own']), {lf.name}));
        else
            dloc = find(cellfun(@(x) contains(x, [num2str(pp) '_SSVEP']) && contains(x, ['own']), {lf.name}));
        end
        logname = lf(dloc).name;
    end
    
    %% Preprocessing section
    cfg = [];
    cfg.dftfreq = 'yes';
    cfg.dftfreq = [50 100 150];
    cfg.demean = 'yes';
    cfg.dataset = [root_dir '\' partname '.vhdr'];
    dat = ft_preprocessing(cfg);
    
    cfg = [];
    cfg.reref = 'yes';
    cfg.refchannel = 'all';
    dat = ft_preprocessing(cfg, dat);
    
    cfg = [];
    cfg.trialdef.prestim = 2;
    cfg.trialdef.poststim = 0;
    cfg.trialdef.eventtype = 'Stimulus';
    cfg.trialdef.eventvalue = {'S107'; 'S108';'S109'; 'S110'; 'S111'; 'S112'; 'S113'; 'S115'; 'S117'; 'S120'};
    cfg.dataset= [root_dir '\' partname '.vhdr'];
    trialdef = ft_definetrial(cfg);
    
    cfg = [];
    cfg.trl = trialdef.trl;
    dat = ft_redefinetrial(cfg, dat);
    
    cfg = [];
    dat = ft_resampledata(cfg, dat);
    if pp == 29
        cfg = [];
        cfg.trials = 2:121;
        dat = ft_redefinetrial(cfg,dat);
    end
    
    % add the trial information of the logfile:
    dat = AddLFinfo([loclog logname], dat);
    
    %%
    vislabel = {'Cz';'FC1';'FC2';'CP1';'CP2'};
    %vislabel = {'Oz';'O1';'O2'};
    
    %cfg =[ ];
    %cfg.trials = 1:120;
    datn = dat;
    
    % automatic removal of bad trials. (only for the relevant channel)
    [val inla] = intersect(datn.label, vislabel);
    
    clear badtr
    for lv = 1:length(inla)
        d = cellfun(@(x) x(inla(lv),:)', datn.trial, 'uniformoutput', 0);
        d = cell2mat(d);
        allv = var(d);
        badtr(lv,:) = allv > mean(allv) + std(allv).*2.5 | allv < mean(allv) - std(allv).*2.5;
    end
    goodtr = find(sum(badtr) == 0);
    
    cfg =[ ];
    cfg.trials = goodtr;
    datn = ft_redefinetrial(cfg, datn);
    
    % now do the power and phase analysis:
    foi = 60./[7 8 9 10 11 12 13 15 17 20];
    int = [7 8 9 10 11 12 13 15 17 20];
    for it = 1:length(foi)
        cfg = [];
        cfg.trials = find(datn.trialinfo == int(it)+100);
        cfg.toilim = [-5/foi(it) 0];
        dattemp = ft_redefinetrial(cfg, datn);
        
        cfg = [];
        %cfg.pad = 'nextpow2';
        cfg.foilim = [1 25];
        cfg.method = 'fft';
        cfg.output = 'fourier';
        cfg.taper = 'hanning';
        %cfg.channel = inla;
        frout{it} = ft_freqanalysis(cfg, dattemp);
        frout{it}.angle = squeeze(circ_r(angle(frout{it}.fourierspctrm)));
        frout{it}.power = squeeze(mean(abs(frout{it}.fourierspctrm).^2));
        for ch = 1:length(vislabel)
            myfit = fittype('a + b*log(x)',...
                'dependent',{'y'},'independent',{'x'},...
                'coefficients',{'a','b'});
            [f] = fit(frout{it}.freq',frout{it}.power(ch,:)',myfit);
            coeffvalues(f);
            frout{it}.powmin1f(ch,:) = frout{it}.power(ch,:)-f(frout{it}.freq)';
        end
    end
    save(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\SSVEP\' partname], 'dat', 'frout');
    
    %% plot here
    ac = {'r';'m';'c';'g';'b';'r:';'m:';'c:';'g:';'b:'};
    for it = 1:length(foi)
        acind = find(frout{it}.freq > foi(it)-1 & frout{it}.freq < foi(it)+1);
        negind = find(frout{it}.freq > foi(it)-2 & frout{it}.freq < foi(it)-1);
        posind = find(frout{it}.freq > foi(it)+1 & frout{it}.freq < foi(it)+2);
        peakvITC(it,:) = (mean(frout{it}.angle(:,acind),2) - mean(frout{it}.angle(:,[negind posind]),2))./mean(frout{it}.angle(:,[negind posind]),2);
        peakvPOW(it) = mean((mean(frout{it}.power(inla,acind),2) - mean(frout{it}.power(inla,[negind posind]),2))./mean(frout{it}.power(inla,[negind posind]),2));
        v = nearest(frout{it}.freq, foi(it));
        vITC(it,:) = frout{it}.angle(:,v);
        %vPOW(it) = mean(frout{it}.power(v),1);
        vPOW(it,:) = frout{it}.power(:,v);
    end
    
    figure
    subplot(1,3,1)
    for it = 1:length(foi)
        plot(frout{it}.freq, mean(frout{it}.power(inla,:),1), ac{it}, 'linewidth', 2);
        hold on
    end
    set(gca, 'xlim', [2.5 9]);
    legend(num2str(round(foi'*10)/10), 'location', 'northeast');
    title('power');
    
    subplot(1,3,2)
    for it = 1:length(foi)
        plot(frout{it}.freq, mean(frout{it}.angle(inla,:),1), ac{it}, 'linewidth', 2);
        hold on
    end
    set(gca, 'xlim', [2.5 9]);
    legend(num2str(round(foi'*10)/10), 'location', 'southeast');
    title('ITC');
    
    % subplot(1,4,3)
    % plot(foi,mean(peakvITC(:,inla),2))
    % title('ITCdiffBas');
    subplot(1,3,3)
    plot(foi,mean(vITC(:,inla),2))
    set(gca, 'xlim', [2.5 9]);
    title('ITCabs');
    set(gcf, 'position', [300 200 1100 500]);
    print([figloc '\' partname 'avgs'],  '-dpng');
    exportfig(gcf, [figloc '\' partname 'avgs']);
    
    %% topography:
    load([locfieldtrip '\template\layout\easycapM1']);
    pldat.label = frout{1}.label;
    pldat.freq = foi(end:-1:1);
    pldat.dimord = 'chan_freq';
    pldat.powspctrm = vITC(end:-1:1,:)';
    %pldat.powspctrm = vPOW(end:-1:1,:)';
    
    figure
    for it = 1:length(foi)
        subplot(2,5,it);
        cfg = [];
        cfg.layout = lay;
        %cfg.zlim = [0.1 0.6];
        cfg.xlim = [pldat.freq(it)-0.1 pldat.freq(it)+0.1];
        ft_topoplotER(cfg, pldat);
    end
    set(gcf, 'position', [300 200 1200 650]);   
    print([figloc '\' partname 'topo'],  '-dpng');
    exportfig(gcf, [figloc '\' partname 'topo']);
    
    %% max ITC freq:
    clc
    [mv inx] = max(mean(vITC(:,inla),2));
    display(['maximum value at frequency ' num2str(foi(inx))])
    display(['corresponds to frame interval '  num2str(int(inx))])           
end

%% show histogram of used freq
close all
load('D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Eff_SSVEP');

figure
histogram(freq,10,'FaceAlpha', 1); xlabel('freq (Hz)'); ylabel('amount');
set(gcf, 'position', [300 200 280 200]);

print([figloc '\histSSVEP'],  '-dpng');
exportfig(gcf, [figloc '\histSSVEP']);   
    
%% show topo of used freq
ppus = [1:37 41:43];
load('D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Eff_SSVEP');
foi = 60./[7 8 9 10 11 12 13 15 17 20];

figure('units','normalized','outerposition',[0 0 1 1])
for pp = 1:length(ppus)    
    partname = [num2str(ppus(pp)) '_SSVEP'];   
    load(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\SSVEP\' partname], 'dat', 'frout');
    
    % nearest freq
    frin = nearest(foi, freq(pp));
    v = nearest(frout{frin}.freq, freq(pp));
    
    % topography:
    load([locfieldtrip '\template\layout\easycapM1']);
    clear pldat
    pldat.label = frout{1}.label;
    pldat.freq = 1;
    pldat.dimord = 'chan_freq';
    pldat.itc = frout{frin}.angle(:,v);
    pldat.pow = squeeze(mean(log(abs(frout{v}.fourierspctrm(:,:,v)).^2)))';
    alld.itc(pp,:) = frout{frin}.angle(:,v);
    alld.pow(pp,:) = squeeze(mean(log(abs(frout{v}.fourierspctrm(:,:,v)).^2)));
    pldat.label{1} = 'Fp1';
    pldat.label{2} = 'Fp2';  
    
    subplot(7,6,pp);
    cfg = [];
    cfg.parameter ='itc';
    cfg.layout = lay;
    cfg.xlim = [0.9 1.1];
    ft_topoplotER(cfg, pldat);    
end
print([figloc '\topobestsInd' cfg.parameter],  '-dpng');
exportfig(gcf, [figloc '\topobestsInd' cfg.parameter], 'color', 'rgb')

%%
figure
pl = {'pow';'itc'};
pl
for it = 1:2
    subplot(1,2,it)
    pldat.powspctrm = squeeze(mean(alld.(pl{it}),1))';
    cfg = [];
    cfg.layout = lay;
    cfg.xlim = [0.9 1.1];
    ft_topoplotER(cfg, pldat);
    colorbar
end
set(gcf, 'position', [300 200 600 300]);
print([figloc '\topobestsMean'],  '-dpng');
exportfig(gcf, [figloc '\topobestsMean'], 'Color', 'rgb')
