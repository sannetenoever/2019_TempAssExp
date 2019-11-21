% read the AV data first
clear vars
restoredefaultpath();
addpath(genpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231'));
addpath('D:\Experiments\StartTimeExp\AA_EEG\03_Scripts\AVTIME');
loclog = 'D:\Experiments\StartTimeExp\AA_EEG\02_Data\Logs\';
root_dir = 'D:\Experiments\StartTimeExp\AA_EEG\';
cd(loclog)
vall = dir('*AVTIME*own.txt');
ppnames = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'38';'39';'41';'42';'43'};
%int = [12 15 11 12 10 13 11 12 17 8 15 15 15 11 17 17 11 17 11 11 13 17 10 11 13 11 10 9 9 17 10 9 12 15 17 11 10];
%freq = 60./int;

%% load data AVTIME
for pp = 39:length(ppnames)
    S01_01_EEGAV_preprocessing(ppnames{pp}, 'AVTIME', 1);
    if pp == 34
        S01_01_EEGAV_preprocessing(ppnames{pp}, 'AVTIME_part2', 1);
    elseif pp == 37
        S01_01_EEGAV_preprocessing(ppnames{pp}, 'AVTIME_part2', 1);
    end
end

%% behavioral analyses
S01_02_readbehav_py_group_propcor;
S01_03_readbehav_py_group_thres;
S01_04_effectiveSSVEP;

% then do S02_00_EEGAV_preprocessing3_Aonly
