function S01_01_EEGAV_preprocessing(pp, expname, first)

namefile = [pp '_' expname '.vhdr'];
loc = 'D:\Experiments\StartTimeExp\AA_EEG\04_MidData\AVTIME\PreProcData\';

% preprocessing for one file
cd D:\Experiments\StartTimeExp\AA_EEG\02_Data\EEG\

% first load the layout
load(['D:\Other\Programs\Matlab\fieldtrip-20120904\fieldtrip-20161231\template\layout\easycapM1']);

%% the first things
if first == 1
    if strcmp(pp, '25') & (strcmp(expname, 'AVTIME') | strcmp(expname, 'SSVEP'));
        namefile = [pp '_SSVEP.vhdr'];
        eventdata = ft_read_event(namefile);
        % find blocks:
        difer = diff([eventdata.sample]);
        [dif tp] = sort(difer);
        if strcmp(expname, 'AVTIME')
            eventdata = eventdata(tp(end-5):end);
        elseif strcmp(expname, 'SSVEP');
            eventdata = eventdata(1:tp(end-5));
        end
    else        
        % create the eventdata
        eventdata = ft_read_event(namefile);
    end
    % get the trialdefinition, one for 12 (and for 16 if it is the tACS
    % dataset)
    cfg = [];
    cfg.event = eventdata;
    cfg.dataset = namefile;
    cfg.trialdef.eventvalue = {'S 15'};
    cfg.trialdef.prestim = 12/3+2;
    cfg.trialdef.poststim = 3;
    cfg.trialdef.eventtype = 'Stimulus';
    cfg1 = ft_definetrial(cfg);
    trlVis = cfg1.trl; % trial definition for events
    
    if strcmp(pp, '38') || strcmp(pp, '39')
       for cnttr = 1:length(eventdata)
           if length(eventdata(cnttr).value) > 0
               if str2double(eventdata(cnttr).value(end)) < 5
                   eventdata(cnttr).value = ['S 1' eventdata(cnttr).value(end)];
               end
           end
       end
    end
    
    cfg = [];
    cfg.event = eventdata;
    cfg.dataset = namefile;
    cfg.trialdef.eventvalue = {'S 11';'S 12';'S 13';'S 14'};
    cfg.trialdef.prestim = 3;
    cfg.trialdef.poststim = 3;
    cfg.trialdef.eventtype = 'Stimulus';
    cfg1 = ft_definetrial(cfg);
    trlAud = cfg1.trl; % trial definition for events
    
    % now the preprocessing demean and notch filter over the whole data !!!!
    cfg = [];
    if pp == '37'
        trlAud = trlAud(1:end-1,:);    
    end
    cfg.trl = [trlAud(1,1)-3000 trlAud(end,2)+3000 0];
    cfg.dataset = namefile;    
    cfg.demean = 'yes';
    cfg.reref = 'yes';
    cfg.refchannel = 'all'; 
    datA = ft_preprocessing(cfg);

    [x indx] = intersect(lower(lay.label), lower(datA.label));
    laynew = struct;
    laynew.pos = lay.pos(indx,:);
    laynew.width = lay.width(indx,:);
    laynew.height = lay.height(indx);
    laynew.label = lay.label(indx);
    laynew.outline = lay.outline;
    laynew.mask = lay.mask;
   
    save([loc '\firstdat_' expname pp], 'datA', 'trlAud', 'trlVis');
    save([loc '\eventdata_' expname pp], 'eventdata'); 
end
