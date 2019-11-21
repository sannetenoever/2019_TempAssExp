function [dat] = AddLFinfo(namefile, dat)

fid = fopen(namefile);
%linesinf = fgetl(fid);linesinf = fgetl(fid);
linesinf = fgetl(fid);
s = textscan(linesinf, '%s');
nums = fscanf(fid, '%f');
nums = reshape(nums, [length(s{1}), length(nums)/length(s{1})]);
for i = 1:length(s{1})
    vars.(s{1}{i}) = nums(i,:);    
end;
relnum = [1 3 4 5]';
for i = 1:length(relnum)
    dat.trialinfo(:,i+1) = nums(relnum(i),:);    
end;
dat.infoattrl = s{1}(relnum);
dat.allinfo = vars;
fclose(fid);


