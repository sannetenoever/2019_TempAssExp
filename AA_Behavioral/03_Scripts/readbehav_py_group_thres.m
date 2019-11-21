%% read the files
clear vars
loclog = 'F:\Experiments\StartTimeExp\DataSC\';
cd(loclog)
vall = dir('*own.txt');
ppnames = {'pp4';'pp5';'pp6';'pp7';'pp8';'pp9';'pp10';'pp11';'pp12';'pp13'; 'pp14';'pp15';'pp16_T';'pp17_T';'pp18_T';'pp19_T';'20_T';'pp21_T';'pp22_T';'pp23_T'};
AS = [1 1 1 1 1 1 2 2 2 2 2 2 2 2 1 2 1 1 1 2];
for ppn = 1:length(ppnames)
    x = find(cellfun(@(x) ~isempty(findstr(ppnames{ppn}, x)), {vall.name}));
    fid = fopen(vall(x).name);
    linesinf = fgetl(fid);
    s = textscan(linesinf, '%s');
    nums = fscanf(fid, '%f');
    nums = reshape(nums, [length(s{1}), length(nums)/length(s{1})]);
    for i = 1:length(s{1})
        vars(ppn).(s{1}{i}) = nums(i,:);
    end;
    fclose(fid);
end;
leg = {'50ms, LP ass', '175ms, HP ass'};
AN = [1 2; 2 1];

%%
figure
for ppn = 1:length(ppnames)
    subplot(4,5,ppn)
    plot(vars(ppn).currentvalue);
    mv(ppn) = mean(vars(ppn).currentvalue(end-50:end)==1);
end
exportfig(gcf, ['F:\Experiments\StartTimeExp\Figure\Thres'], 'Color', 'rgb');