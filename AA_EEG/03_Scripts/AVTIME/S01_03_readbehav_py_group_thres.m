%% read the files
clear
addpath('D:\Other\Programs\Matlab\extrausefullScripts');
loclog = 'D:\Experiments\StartTimeExp\AA_EEG\02_Data\Logs\';
cd(loclog)

vall{1} = dir('*AVTIME*own.txt');
vall{2} = dir('*SSVEP*own.txt');
ppnames = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'41';'42';'43'};
AS = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1 2 1];

for ppn = 1:length(ppnames)
    for v = 1:2        
        x = find(cellfun(@(x) strcmp([ppnames{ppn} '_'], x(1:length(ppnames{ppn})+1)), {vall{v}.name}));
        %x = find(cellfun(@(x) ~isempty(findstr([ppnames{ppn} '_'], x(1:length(ppnames{ppn})+1))), {vall{v}.name}));
        fid = fopen(vall{v}(x).name);
        linesinf = fgetl(fid);
        s = textscan(linesinf, '%s');
        nums = fscanf(fid, '%f');
        nums = reshape(nums, [length(s{1}), length(nums)/length(s{1})]);
        for i = 1:length(s{1})
            vars(v,ppn).(s{1}{i}) = nums(i,:);
        end
        fclose(fid);
    end
end
leg = {'50ms, LP ass', '175ms, HP ass'};
AN = [1 2; 2 1];

%%
cntfig = 0;
figure
for ppn = 1:length(ppnames)
    cntfig = cntfig + 1;   
    subplot(4,4,cntfig)    
    x = [vars(2,ppn).currentvalue vars(1,ppn).currentvalue];
    plot(x);
    mv(ppn) = mean(vars(1,ppn).currentvalue(end-50:end)==1);
    mv(ppn) = var(vars(1,ppn).currentvalue(end-100:end));
    mv(ppn) = min(vars(1,ppn).currentvalue(end-100:end))-max(vars(1,ppn).currentvalue(end-100:end));
    title(ppnames{ppn});
    if cntfig == 16
       cntfig = 0; 
       figure
    end
end
[x y] = sort(mv)

%exportfig(gcf, ['F:\Experiments\StartTimeExp\Figure\Thres'], 'Color', 'rgb');

