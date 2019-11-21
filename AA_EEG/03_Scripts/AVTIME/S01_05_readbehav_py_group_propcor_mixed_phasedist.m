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

%% calculate PD
varnames = {'sound';'timepoint';'ass';'pp';'curval';'dev';'phasedist';'accuracy'};
saveloc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\';
load('D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Eff_SSVEP', 'freq','int');

for pp = 1:length(vars)
    % phase:
    c = 1./freq(pp);
    P(pp) = 0.16./c;
    if P(pp) > 1
        P(pp) = P(pp)-1;
    elseif P(pp) > 2
        P(pp) = P(pp)-2;
    end
    
    % phasedist
    c = 1./freq(pp);
    D = abs(0.16-c);
    %PD(pp) = D;
    PD(pp) = D./c;
    if PD(pp) > 0.5
        PD(pp) = 1-PD(pp);
    end
    OPD(pp) = PD(pp)+0.5;
    if OPD(pp) > 0.5
        OPD(pp) = 1-OPD(pp);
    end
end


%% do multilevel regression at the end
bigmat = [];
bigmat2 = [];
varnames = {'accuracy';'sound';'ass';'pp';'curval'; 'dev'; 'Freq';'FreqDif';'phasedist';'timepoint'};
saveloc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\SaveD\';
load('D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\Eff_SSVEP', 'freq','int');

for pp = 1:length(vars)
    vars(pp).DEV = linspace(0,1, length(vars(pp).SOUND));    
    for fr = 1:2
        for t = 1:2
            inx = find(vars(pp).SOUND2 == fr & vars(pp).ASS == t-1);
            inxT = size(bigmat,2) + 1:size(bigmat,2) + length(inx);
            bigmat(1,inxT) = vars(pp).RESC(inx);
            bigmat(2,inxT) = fr;
            bigmat(3,inxT) = t;
            bigmat(4,inxT) = pp;
            bigmat(5,inxT) = vars(pp).currentvalue(inx);
            bigmat(6,inxT) = vars(pp).DEV(inx); % 0 = beginning
            bigmat(7,inxT) = freq(pp);
            bigmat(8,inxT) = abs(4-freq(pp));
            bigmat(9,inxT) = PD(pp);
            bigmat(10,inxT) = vars(pp).TIMEPOINT(inx);
            
            bigmat2(1,end+1) = mean(vars(pp).RESC(inx));
            bigmat2(2,end) = fr;
            bigmat2(3,end) = t;
            bigmat2(5,end) = pp;
        end
    end    
end
%bigmat(:,ismember(bigmat(5,:),ol)) = [];
bigmat(6,:) = zscore(bigmat(6,:)); % zscore difficulty
bigmat = bigmat';
bigmat2 = bigmat2';
save([saveloc 'S01_05_BIGMAT'], 'bigmat');
save([saveloc 'S01_05_BIGMAT_varnames'], 'varnames');

%%
for pp = 1:length(ppnames)
    for st = 1:2
        x = find(bigmat(:,4) == pp & bigmat(:,2) == st);
        me(st,pp) = mean(bigmat(x,1));
    end
end
[x ia] = sort(PD);
plot(me(2,ia)-me(1,ia))

plot(me(:,ia)')



%% see over all blocks
clear proplp amtr meanv dif p allv
tp = [50 175];
blc = {[0 1]; [2 3 4]};

for pp = 1:length(vars)
    for bl = 1:size(blc,1)
        for fr = 1:2
            for t = 1:2
                %inx = find(vars(pp).SOUND2 == fr & vars(pp).ASS == t-1 & ismember(vars(pp).BLOCK, blc{bl}));
                inx = find(vars(pp).SOUND2 == fr & vars(pp).TIMEPOINT == tp(t) & ismember(vars(pp).BLOCK, blc{bl}));
                amlp = sum(vars(pp).RESC(inx));
                allv.proplpA(pp,bl,fr,t) = amlp./length(inx);
                allv.CVA(pp,bl,fr,t) = mean(vars(pp).currentvalue(inx));
            end
            allv.proplp(pp,bl,fr) = allv.proplpA(pp,bl,fr,2)-allv.proplpA(pp,bl,fr,1);
            allv.amtr(pp,bl,fr) = length(inx);
            allv.pp(pp,bl,fr) = pp;
            allv.bl(pp,bl,fr) = bl;
            allv.t(pp,bl,fr) = t;
            allv.s(pp,bl,fr) = fr;
            allv.pd(pp,bl,fr) = PD(pp);
            allv.CV(pp,bl,fr) = squeeze(mean(allv.CVA(pp,bl,fr,:),4));
        end
    end
end

ppto = 1:length(vars);

%% put sound back in

%% do the RM anova
restoredefaultpath();
addpath('D:\Other\Programs\Matlab\extrausefullScripts\');
ft.propcor = allv.proplp(ppto,:,:);
ft.propcor = allv.CV(ppto,:,:);
ft.pp = allv.pp(ppto,:,:);
ft.pd = allv.pd(ppto,:,:);
ft.block = allv.bl(ppto,:,:);
ft.fr = allv.s(ppto,:,:);

tb = table(ft.pp(:), ft.pd(:), ft.block(:),ft.fr(:), ft.propcor(:), 'VariableNames', {'pp', 'PD', 'BL','SND','pCOR'});
tb = uni2multiTB(tb, [3 4], 5, [1 2]);
F1 = [1 1 2 2]';
F2 = [1 2 1 2]';
withdes = table(F1,F2, 'VariableNames', {'BLOCK';'SND';});
withdes.BLOCK = categorical(withdes.BLOCK);
withdes.SND = categorical(withdes.SND);

rm = fitrm(tb,'BL1SND1-BL2SND2~1+PD','WithinDesign',withdes);
ranovatbl = ranova(rm, 'WithinModel','BLOCK*SND');

% independent of block
rm = fitrm(tb,'BL1SND1-BL2SND2~1+PD','WithinDesign',withdes);
ranovatbl = ranova(rm, 'WithinModel','SND');

%% only for first timepoint
F1 = [1 2]';
withdes = table(F1, 'VariableNames', {'BLOCK'});
withdes.BLOCK = categorical(withdes.BLOCK);

%tp 1
rm = fitrm(tb,'BL1TP1,BL2TP1~1+PD','WithinDesign',withdes);
ranovatbl = ranova(rm, 'WithinModel','BLOCK');

%tp 2
rm = fitrm(tb,'BL1TP2,BL2TP2~1+PD','WithinDesign',withdes);
ranovatbl = ranova(rm, 'WithinModel','BLOCK');

% just correlation
[CV p] = corr(PD', squeeze(mean(allv.proplp(ppto,1,:),3)))
scatter(PD', squeeze(mean(allv.proplp(ppto,1,:),3)))

[CV p] = corr(PD', squeeze(mean(allv.CV(ppto,1,:),3)))

