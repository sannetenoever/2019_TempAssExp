ppnames = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'41';'42';'43'};
expname = 'AVTIME';
loc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\PreProcData\';
load('D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Eff_SSVEP', 'freq','int');

for ppc = 1:length(ppnames)    
    pp = ppnames{ppc};
    load([loc '\firstdat_' expname pp], 'datA', 'trlAud', 'trlVis');
    load([loc '\eventdata_' expname pp], 'eventdata');
    load([loc '\eyeblinkcorrectionRMT_' expname pp], 'datAud_preproc_RM', 'preICA');
    
    toi = find(datAud_preproc_RM.time{1} < 0 & datAud_preproc_RM.time{1} > -5/freq(ppc));
    clear eb
    for tc = 1:length(datAud_preproc_RM.trial)
        e = squeeze(preICA.trial(tc,32,:));
        datAud_preproc_RM.trial{tc}(32,:) = e
        datAud_preproc_RM.trial{tc}(33,:) = squeeze(preICA.trial(tc,1,:));
        eb(tc) = isempty(find(e(toi) < -100)); 
    end
    datAud_preproc_RM.label{32} = 'eye';
    datAud_preproc_RM.label{33} = 'preFP1';   
            
%     cfg = [];
%     cfg.viewmode = 'vertical';
%     ft_databrowser(cfg, datAud_preproc_RM);
    alltr{ppc} = eb;
end
save('D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\S02_01_blinktrl', 'alltr');


