%% read the files
clear vars
addpath('D:\Other\Programs\Matlab\extrausefullScripts');
loclog = 'D:\Experiments\StartTimeExp\AA_EEG\02_Data\Logs\';
cd(loclog)
vall = dir('*AVTIME*own.txt');
ppnames = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'41';'42';'43'};
AS = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1 2 1];
allfr = {[1:100], [102:201]};
allfrab17 = {[1:201], [203:400]};
AN = [1 2; 2 1];
for ppn = 1:length(ppnames)
    x = find(cellfun(@(x) strcmp([ppnames{ppn} '_'], x(1:length(ppnames{ppn})+1)), {vall.name}));
    fid = fopen(vall(x).name);
    linesinf = fgetl(fid);
    s = textscan(linesinf, '%s');
    nums = fscanf(fid, '%f');
    nums = reshape(nums, [length(s{1}), length(nums)/length(s{1})]);
    for i = 1:length(s{1})
        vars(ppn).(s{1}{i}) = nums(i,:);
    end
    fclose(fid);
    vars(ppn).SOUND2 = ones(size(vars(ppn).SOUND));
    vars(ppn).RESC = zeros(size(vars(ppn).SOUND));   
    if ppn >= 17
        af = allfrab17;
    else
        af = allfr;
    end
    inx = find(ismember(vars(ppn).SOUND, af{AN(AS(ppn),1)}));    
    vars(ppn).SOUND2(inx) = 1;
    inx = find(ismember(vars(ppn).SOUND, af{AN(AS(ppn),2)}));
    vars(ppn).SOUND2(inx) = 2;
    inx = find(ismember(vars(ppn).SOUND, af{1}) & vars(ppn).RESPONSE == 1);    
    vars(ppn).RESC(inx) = 1;
    inx = find(ismember(vars(ppn).SOUND, af{2})& vars(ppn).RESPONSE == 2);   
    vars(ppn).RESC(inx) = 1;
    vars(ppn).ASS = zeros(size(vars(ppn).SOUND));
    vars(ppn).ASS((vars(ppn).SOUND2 == 1 & vars(ppn).TIMEPOINT == 50) | (vars(ppn).SOUND2 == 2 & vars(ppn).TIMEPOINT == 175)) = 1;
    
    ass{1} = [175 75; 75 175];
    ass{2} = [75 175; 175 75];
    va(1,1) = sum(vars(ppn).SOUND2 == 1 & vars(ppn).TIMEPOINT == 50);
    va(2,1) = sum(vars(ppn).SOUND2 == 2 & vars(ppn).TIMEPOINT == 50);
    va(1,2) = sum(vars(ppn).SOUND2 == 1 & vars(ppn).TIMEPOINT == 175);
    va(2,2) = sum(vars(ppn).SOUND2 == 2 & vars(ppn).TIMEPOINT == 175);
    if sum(ass{1}-va) == 0
        vars(ppn).ass = 1;
    elseif sum(ass{2}-va) == 0
        vars(ppn).ass = 2;
    else
        vars(ppn).ass = 0;
    end
end
save(['D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Behdata.mat'], 'ppnames', 'vars');
leg = {'NonassoTP', 'AssoTP'};

%% see over all blocks
clear proplp amtr meanv dif p
tp = [50 175];
blc = {[0 1]; [2 3 4]};

for pp = 1:length(vars)
    for bl = 1:size(blc,1)
        for fr = 1:2
            for t = 1:2
                inx = find(vars(pp).SOUND2 == fr & vars(pp).ASS == t-1 & ismember(vars(pp).BLOCK, blc{bl}));            
                amlp = sum(vars(pp).RESC(inx));
                allv.proplp(pp,bl,fr,t) = amlp./length(inx);
                meanv(pp,bl) = mean(vars(pp).currentvalue(inx));
                allv.amtr(pp,bl,fr,t) = length(inx);
                allv.pp(pp,bl,fr,t) = pp;
                allv.bl(pp,bl,fr,t) = bl;
                allv.t(pp,bl,fr,t) = t;
                allv.s(pp,bl,fr,t) = fr; 
            end
        end
    end
end

% check for outliers
% outliers overall
d = squeeze(mean(mean(mean(allv.proplp,4),3),2));
ol1 = [find(d > median(d) + 1.5*iqr(d)); find(d < median(d) - 1.5*iqr(d))];

% check for outliers for average of main effect (average over both sounds)
v = squeeze(allv.proplp(:,2,:,2)-allv.proplp(:,2,:,1));
d = mean(v,2);
ol = unique([ol1; find(d > median(d) + 1.5*iqr(d)); find(d < median(d) - 1.5*iqr(d))]);
ppto = setdiff(1:length(vars),ol);

figure
tl = {'First half';'Second half'};
for bl = 1:2
    subplot(1,2,bl)
    bar(squeeze(mean(allv.proplp(ppto,bl,:,:),1)));
    hold on
    errorbar([0.85 1.15; 1.85 2.15], squeeze(mean(allv.proplp(ppto,bl,:,:),1)), squeeze(std(allv.proplp(ppto,bl,:,:),[],1))./sqrt(length(ppto)), 'k','linestyle','none', 'linewidth',3)
    set(gca, 'ylim', [0.5 0.65], 'xticklabel', {'Sound A';'Sound B'})
    if bl == 1
        legend(leg, 'location', 'southeast');
    end
    ylabel('proportion correct');
    title(tl{bl})
    set(gcf, 'position', [500 200 600 250]);
end
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S01_02_BehResults'], 'Color', 'rgb');
print('D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S01_02_BehResults','-dpng')

%% do multilevel regression at the end
bigmat = [];
bigmat2 = [];
varnames = {'accuracy';'sound';'ass';'block';'pp';'curval'; 'dev'; 'Freq';'FreqDif';'ExcDis'};
saveloc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\';
load('D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Eff_SSVEP', 'freq','int');

for pp = 1:length(vars)
    vars(pp).DEV = linspace(0,1, length(vars(pp).SOUND));
    for bl = 1:size(blc,1)
        for fr = 1:2
            for t = 1:2
                inx = find(vars(pp).SOUND2 == fr & vars(pp).ASS == t-1 & ismember(vars(pp).BLOCK, blc{bl}));            
                inxT = size(bigmat,2) + 1:size(bigmat,2) + length(inx);  
                bigmat(1,inxT) = vars(pp).RESC(inx);                
                bigmat(2,inxT) = fr;
                bigmat(3,inxT) = t;
                bigmat(4,inxT) = bl;
                bigmat(5,inxT) = pp;
                bigmat(6,inxT) = vars(pp).currentvalue(inx); 
                bigmat(7,inxT) = vars(pp).DEV(inx); % 0 = beginning
                bigmat(8,inxT) = freq(pp);
                bigmat(9,inxT) = abs(4-freq(pp));
                bigmat(10,inxT) = 0.16./(1./freq(pp));
                if 0.16./(1./freq(pp)) > 0.5
                    bigmat(10,inxT) = 1-(0.16./(1./freq(pp)));
                end
                  
                bigmat2(1,end+1) = mean(vars(pp).RESC(inx)); 
                bigmat2(2,end) = fr;
                bigmat2(3,end) = t;
                bigmat2(4,end) = bl;
                bigmat2(5,end) = pp;
            end
        end
    end
end
bigmat(:,ismember(bigmat(5,:),ol)) = [];
bigmat(6,:) = zscore(bigmat(6,:)); % zscore difficulty
bigmat = bigmat';
bigmat2 = bigmat2';
save([saveloc 'S01_02_BIGMAT'], 'bigmat');
save([saveloc 'S01_02_BIGMAT_varnames'], 'varnames');
save([saveloc 'S01_02_outlierdef'], 'ol');

%% do the RM anova
restoredefaultpath();
addpath('D:\Other\Programs\Matlab\extrausefullScripts\');
ft.propcor = logit(allv.proplp(ppto,:,:,:));
ft.pp = allv.pp(ppto,:,:,:);
ft.block = allv.bl(ppto,:,:,:);
ft.associated = allv.t(ppto,:,:,:);
ft.sound = allv.s(ppto,:,:,:);

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

%% visualization of individual pp
cb(1,:) = [1 0 0];
cb(2,:) = [0 0 0];
dc = ol;
c = ones(1,size(v,1));
c(dc) = 0;
close all
scatter(v(:,1), v(:,2),200,c','.'); 
colormap(cb); 
set(gca, 'xlim', [-0.3 0.3], 'ylim', [-0.3 0.3]);
lx = get(gca, 'xlim');
ly = get(gca, 'ylim');
hold on
plot([0 0], ly, 'k:');
plot(lx, [0 0], 'k:');
plot(lx([2 1]), ly, 'k', 'linewidth',2);
xlabel('NonAsso better <    Sound A  > Asso better');ylabel('NonAsso better <   Sound B  > Asso better');
title(['Proportion difference NonAsso-Asso']);

exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S01_02_BehScatterResults'], 'Color', 'rgb');
print('D:\Experiments\StartTimeExp\AA_EEG\05_Figures\AVTIME\S01_02_BehScatterResults','-dpng')
