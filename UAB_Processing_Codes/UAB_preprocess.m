function [out] = UAB_preprocess(dirname,ccep_data_path,out_path)
%code to re-reference and epoch data; saved in an out file for each

% addpath '/data/project/NSPMlab/hbriny99/'
path = string(append(ccep_data_path,'/CCEPs/',dirname));
cd(path);
stimlist = dir(append('*','and','*'));
for n = 1:length(stimlist) %each stimulus
    stim = stimlist(n).name;
    newpath = string(append(path,'/',stim));
    cd(newpath);
    load('data.mat');
    dataTable = data;
    VariableNames = data.Properties.VariableNames;
    for k = 1:length(VariableNames)
        if contains(VariableNames(k), '0') 
            Eidx(k) = 1;
        elseif contains(VariableNames(k), '1')
            Eidx(k) = 1;
        elseif contains(VariableNames(k), '2')
            Eidx(k) = 1;
        elseif contains(VariableNames(k), '3')
            Eidx(k) = 1;
        elseif contains(VariableNames(k), '4')
            Eidx(k) = 1;
        elseif contains(VariableNames(k), '5')
            Eidx(k) = 1;
        elseif contains(VariableNames(k), '6')
            Eidx(k) = 1;
        elseif contains(VariableNames(k), '7')
            Eidx(k) = 1;
        elseif contains(VariableNames(k), '8')
            Eidx(k) = 1;
        elseif contains(VariableNames(k), '9')
            Eidx(k) = 1;
        else
            Eidx(k) = 0;
        end
    end
    pt_name = string(dirname);
    pt_name = strsplit(pt_name,' ');
    pt_name = pt_name(1);
    out.name = pt_name;
    electrodenames = VariableNames(Eidx==1);
    Data = data.Variables;
    sz1 = size(Data,2);
    s = length(Data{1,1});
    sz2 = s*size(Data,1);
    megaData = zeros(sz2,sz1);
    for j = 1:size(Data,1)
        jf = j*s;
        js = jf-(s-1);
        for m = 1:size(Data,2)
            megaData(js:jf,m) = Data{j,m};
        end
    end
    megaData = megaData(:,Eidx==1);
    Stim = strsplit(stim,' and ');
    Stim = string(Stim);
    %check for ' (primes)
    if contains(Stim(1),"'")
        Stim(1) = strrep((Stim(1)),"'","p");
        Stim(2) = strrep((Stim(2)),"'","p");
    end
    Stim = append(string(Stim(1)),'_',string(Stim(2)));
    OGdata.(Stim) = megaData; %Somewhat organized
    %% re-reference
    for o = 1:length(electrodenames)%for each stim channel
        object = electrodenames(o);
        object = regexprep(object,'\d','');
        object = erase(object,'_');
        prelimElectrodes(o) = object;
    end
    electrodes = unique(prelimElectrodes,'stable');
    stim_chans = fieldnames(OGdata);
    removeidx = zeros(1,length(prelimElectrodes));
    for remove = 1:length(prelimElectrodes)
        label = prelimElectrodes(remove);
        if ismember(label,'C')
            removeidx(remove) = 1;
        elseif ismember(label,'DC')
            removeidx(remove) = 1;
        end
    end
    removeidx = find(removeidx == 1);
    removeidx2c = find(ismember(electrodes, "C"));
    removeidx2dc = find(ismember(electrodes, "DC"));
    removeidx2 = [removeidx2c, removeidx2dc];
    prelimElectrodes(removeidx) = [];
    electrodes(removeidx2) = [];
    electrodenames(removeidx) = [];
    resp_e_names = electrodenames;
    for r = 1:length(resp_e_names)
        resp_e_names(r) = erase(resp_e_names(r),'_');
    end
    for l = 1:length(stim_chans)%each stim
        SC = string(stim_chans(l));
        Ordered_data = OGdata.(SC);
        Ordered_data(:,removeidx) = [];%removes non recording channels C and DC
        for p = 1:length(electrodes)%each response 
            data2grabidx.(char(electrodes(p))) = find(ismember(prelimElectrodes,electrodes(p)));
            Electrodes.(char(electrodes(p))) = Ordered_data(:,data2grabidx.(char(electrodes(p))));%Breaks up data by electrodes
        end
        z = 0;
        y = 0;
        for q = 1:length(electrodes)
            TempData = Electrodes.(char(electrodes(q)));
            for n = 1:(size(TempData,2)-1)
                RerefData(:,n) = TempData(:,n+1) - TempData(:,n); % "L2 - L1"
                str2 = resp_e_names(y+n+1);
                str1 = resp_e_names(y+n);
                BipolarRespLabels(z+n) = append(str2,'_',str1);
            end
            z = z+n;
            y = y+size(TempData,2);
            BipolarMontage.(char(electrodes(q))) = RerefData;
            clearvars TempData RerefData;
        end
    end
    BipolarMontageTable = struct2table(BipolarMontage);
    rereferenceddata.(Stim) = BipolarMontageTable.Variables;
    out.bipolarresplabels = BipolarRespLabels;
    clearvars BipolarMontageTable;
end
%%
% Extract Epochs
pt_name = erase(pt_name," CCEP");
path3 = string(append(out_path,'/',pt_name));
cd(path3)
load('info.mat')
dp = info.NumSamples(1);
duration = seconds(info.DataRecordDuration);
fs = dp/duration;
out.fs = fs;
stimbank = fieldnames(rereferenceddata);
out.stimlabels = stimbank;
for i = 1:length(stimbank) %for each stim
    stimatmoment = string(stimbank(i));
    datatoepoch = rereferenceddata.(stimatmoment);%matrix of all data recorded for this stim
    %Get EpochIdx
    [EpochIdx, allbad] = getepochidx(i, datatoepoch);
    for nn = 1:size(datatoepoch,2)%each response
        if allbad == 1
            continue; %skip this stim
        end
        for j = 1:length(EpochIdx)%each trial
            datato_epoch = rereferenceddata.(stimatmoment)(:,nn);%stim-response
            b4 = round(EpochIdx(j)-(.5*fs)); % 500 ms before stim
            after = round(EpochIdx(j)+(.9*fs)); %900 ms after stim
            if b4 < 0
                continue;
            end
            if after > length(datato_epoch)
                continue;
            end
            epoch = datato_epoch((b4:after));
            baseline1 = round(.450*fs);
            baseline2 = round(.05*fs);
            baseline = datato_epoch((EpochIdx(j)-baseline1):(EpochIdx(j)-baseline2));
            MeanBaseline = mean(baseline); 
            epoch = epoch- MeanBaseline;
            out.epoch_data.(stimatmoment).(char(BipolarRespLabels(nn)))(j,:) = epoch;
            clearvars epoch
        end
    end
end %for each stim 
save('out',"out",'-v7.3');
save('reref',"rereferenceddata",'-v7.3');