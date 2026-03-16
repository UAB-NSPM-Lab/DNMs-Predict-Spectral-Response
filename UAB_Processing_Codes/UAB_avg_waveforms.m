function [avg_waveforms] = UAB_avg_waveforms(dirname, filtered_data, out_path)

path = string(append(out_path,'/',dirname));
cd(path);
avg_waveforms.name = filtered_data.name;
Snames = fieldnames(filtered_data.filteredmeans);
Rnames = filtered_data.bipolarresplabels;
fs = filtered_data.fs;
avg_waveforms.resplabels = Rnames;
avg_waveforms.fs = filtered_data.fs;
f = ceil(1.4*fs);
for i = 1:length(Snames) %each stim
    S = string(Snames(i));
    Rlist = fieldnames(filtered_data.filteredmeans.(S));
    for j = 1:length(Rnames) %each resp
        R = string(Rnames(j));
        if contains(R,Rlist)
           avg_waveforms.average_waveforms.(S)(:,j) = filtered_data.filteredmeans.(S).(R)(1:f)';
        else
           avg_waveforms.average_waveforms.(S)(:,j) = nan((f),1);
        end
    end
end
savename = string(append(dirname,'_avg_waveforms'));
save(savename,"avg_waveforms",'-v7.3');