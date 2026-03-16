function plotBode(freqs_log,singValues_log,colorLine)

switch nargin
    case 2
        figure;
        plot(freqs_log,singValues_log,'Color','b','LineWidth',3);
        xt = linspace(-3,7,11);
        xticks(xt);
        xticklabels(round(exp(xt)./(2*pi)));
        [maxSingValue,indMaxSingValue] = max(singValues_log);
%         text(freqs_log(indMaxSingValue),maxSingValue,chanLabels(stimChan1));
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');
    case 3
        figure;
        plot(freqs_log,singValues_log,'Color',colorLine,'LineWidth',3);
        xt = linspace(-3,7,11);
        xticks(xt);
        xticklabels(round(exp(xt)./(2*pi)));
        [maxSingValue,indMaxSingValue] = max(singValues_log);
%         text(freqs_log(indMaxSingValue),maxSingValue,chanLabels(stimChan1));
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');
end

end