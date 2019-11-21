ppnames = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'41';'42';'43';};
loc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\PreProcData\';
expname = 'AVTIME';

for it = 1:length(ppnames)
    load([loc '\eventdata_' expname ppnames{it}], 'eventdata');
    ind = find(cellfun(@(x) strcmp(x, 'S  5'), {eventdata.value}));
    int(it) = round(diff([eventdata(ind(end-1:end)).sample]./2500)./(1/60));
end

freq = 1./(int./60);

save('D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Eff_SSVEP', 'freq','int');