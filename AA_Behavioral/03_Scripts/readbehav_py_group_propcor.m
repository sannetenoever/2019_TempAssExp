%% read the files
clear 
addpath(genpath('D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231'));
addpath('D:\Other\Programs\Matlab\extrausefullScripts');
loclog = 'D:\Experiments\StartTimeExp\AA_Behavioral\02_Data\DataSC\';
cd(loclog)
vall = dir('*own.txt');
ppnames = {'pp4';'pp5';'pp6';'pp7';'pp8';'pp9';'pp10';'pp11';'pp12';'pp13'; 'pp14';'pp15';'pp16';'pp17';'pp18';'pp19';'20';'pp21';'pp22';'pp23';'ppR1';'ppR2';'ppR3';'ppR4';'ppR5';'ppR6';'ppR7';'ppR8';'ppR9';'ppR10'};
AS = [1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 2 1 1 1 2 1 2 2 1 1 2 2 2 1 1]; % association of the sounds
af = {[1:100], [102:201]};
AN = [1 2; 2 1];
for ppn = 1:length(ppnames)
    x = find(cellfun(@(x) ~isempty(findstr([ppnames{ppn} '_T'], x)), {vall.name}));
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
end
leg = {'NonassoTP', 'AssoTP'};

%% see over all blocks
clear proplp amtr meanv dif p
allfr = {[1:100], [102:201]};
tp = [50 175];
blc = {[0 1 2]; [3 4 5]};

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
v = squeeze(allv.proplp(:,2,:,2)-allv.proplp(:,2,:,1));
d = v(:,2);
ol = [find(d > mean(d) + 2.5*std(d)) find(d < mean(d) - 2.5*std(d))];
d = v(:,1);
ol = [ol find(d > mean(d) + 2.5*std(d)) find(d < mean(d) - 2.5*std(d))]

ppto = setdiff(1:length(vars),ol);

figure
tl = {'first part';'second part'};
for bl = 1:2
    subplot(1,2,bl)
    bar(squeeze(mean(allv.proplp(ppto,bl,:,:),1)));
    hold on
    errorbar([0.85 1.15; 1.85 2.15], squeeze(mean(allv.proplp(ppto,bl,:,:),1)), squeeze(std(allv.proplp(ppto,bl,:,:),[],1))./sqrt(length(ppto)), 'k','linestyle','none', 'linewidth',3)
    set(gca, 'ylim', [0.55 0.7], 'xticklabel', {'SOUND A';'SOUND B'})
    if bl == 1
        legend(leg, 'location', 'southeast');
    end
    ylabel('proportion correct');
    title(tl{bl})
    set(gcf, 'position', [500 200 600 250]);
end
exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_Behavioral\05_Figures\Results'], 'Color', 'rgb');
print('D:\Experiments\StartTimeExp\AA_Behavioral\05_Figures\Results','-dpng')

%% do the RM anova
restoredefaultpath();
addpath('D:\Other\Programs\Matlab\extrausefullScripts\');
ft.propcor = allv.proplp(ppto,:,:,:);
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

%% follow up on the block*ass interaction.
% same analysis but for the blocks seperately:
F2 = [1 1 2 2]';
F3 = [1 2 1 2]';
withdes = table(F2,F3, 'VariableNames', {'ASS';'SOUND'});
withdes.ASS = categorical(withdes.ASS);
withdes.SOUND = categorical(withdes.SOUND);

%first 3 blocks
rm = fitrm(tb,'BL1ASS1S1-BL1ASS2S2~1','WithinDesign',withdes);
ranovatbl = ranova(rm, 'WithinModel','ASS*SOUND');

%last 3 blocks
rm = fitrm(tb,'BL2ASS1S1-BL2ASS2S2~1','WithinDesign',withdes);
ranovatbl = ranova(rm, 'WithinModel','ASS*SOUND');

% just to ensure that both seperate are significant (although no
% interaction so strickly not necessary)
F2 = [1 2]';
withdes = table(F2, 'VariableNames', {'ASS'});
withdes.ASS = categorical(withdes.ASS);

% tp 1 associated:
rm = fitrm(tb,'BL2ASS1S1, BL2ASS2S1~1','WithinDesign',withdes);
ranovatbl = ranova(rm, 'WithinModel','ASS');
% tp 2 associated:
rm = fitrm(tb,'BL2ASS1S2, BL2ASS2S2~1','WithinDesign',withdes);
ranovatbl = ranova(rm, 'WithinModel','ASS');

%% visualization of individual pp
cb(1,:) = [1 0 0];
cb(2,:) = [0 0 0];
dc = ol;
c = ones(1,size(v,1));
c(dc) = 0;
close all
scatter(v(:,1), v(:,2),200,c','.'); 
colormap(cb); 
set(gca, 'xlim', [-0.25 0.25], 'ylim', [-0.25 0.25]);
lx = get(gca, 'xlim');
ly = get(gca, 'ylim');
hold on
plot([0 0], ly, 'k:');
plot(lx, [0 0], 'k:');
plot(lx([2 1]), ly, 'k', 'linewidth',2);
xlabel('NonAsso better <    Sound A  > Asso better');ylabel('NonAsso better <   Sound B  > Asso better');
title(['Proportion difference NonAsso-Asso']);

exportfig(gcf, ['D:\Experiments\StartTimeExp\AA_Behavioral\05_Figures\ScatterResults'], 'Color', 'rgb');
print('D:\Experiments\StartTimeExp\AA_Behavioral\05_Figures\ScatterResults','-dpng')
