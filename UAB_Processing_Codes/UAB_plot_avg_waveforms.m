function NYU_plot_avg_waveforms(avg_waveforms)
%% This code is to plot the average waveforms 

r_length = length(fieldnames(avg_waveforms.average_waveforms));
plot_sz = ceil(sqrt(r_length));
stims = fieldnames(avg_waveforms.average_waveforms);
resplist = avg_waveforms.resplabels;
figure;
for i = 1:length(stims)%each stim
    S = string(stims(i));
    subplot(plot_sz,plot_sz,i);
    for ii = 1:length(resplist)
            plot(avg_waveforms.average_waveforms.(S)(:,ii));
            hold on;
    end
    T = string(S);
    title(T,Interpreter="none");
    T2 = string(avg_waveforms.name);
    sgtitle(T2);
end