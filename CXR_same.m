%Calculates CXRsame
clearvars
cd('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/CCSR_Resonance_Comparison')
load("SplitCorr.mat")
pt_list = string(fieldnames(SplitCorr));
pt_list_ignore = 2;
pt_list(pt_list_ignore) = [];
for i = 1:length(pt_list)
    a = 0;
    pt = pt_list(i);
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt));
    load("base.mat");
    load("ArtifactualChannels.mat")
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt,'/TFMsBode_figures2/Bode-CCSR Freq Comparison'))
    load("ActualDomFreqs.mat");
    load("AllPeaks_Bode.mat");
    s_list = string(fieldnames(Arts));
    r_list = base.bipolarresplabels;
    S = fieldnames(Arts);
    M_avg = zeros(length(S),length(Arts.(string(S(1)))));
    for m = 1:length(S)
        M_avg(m,:) = Arts.(string(S(m)));
    end
    F = find(mean(M_avg,1) == 1);
    r_list(F) = [];
    for j = 1:length(s_list)
        S = s_list(j);
        r_cand_list = fieldnames(AllPeaks_Bode.(S));
        for k = 1:length(r_list)
            K = k*3-2;
            R = string(r_list(k));
            if ~ismember(R,r_cand_list)
                FreqMatches.(pt).Avg{j,k} = "NaN";
                avglayer(j,K) = nan;
                avglayer(j,K+1) = nan;
                avglayer(j,K+2) = nan;
                FreqMatches.(pt).N1{j,k} = "NaN";
                N1layer(j,K) = nan;
                N1layer(j,K+1) = nan;
                N1layer(j,K+2) = nan;
                FreqMatches.(pt).N2{j,k} = "NaN";
                N2layer(j,K) = nan;
                N2layer(j,K+1) = nan;
                N2layer(j,K+2) = nan;
                FreqMatches.(pt).Rand{j,k} = "NaN";
                Randlayer(j,K) = nan;
                Randlayer(j,K+1) = nan;
                Randlayer(j,K+2) = nan;
                FreqMatches.(pt).Null{j,k} = "NaN";
                Nulllayer(j,K) = nan;
                Nulllayer(j,K+1) = nan;
                Nulllayer(j,K+2) = nan;
                continue;
            end
            bode_freqs = round(AllPeaks_Bode.(S).(R).freq_Hz);%in Hz
            if isnan(ActualDomFreqs.Avg{j,k})
                FreqMatches.(pt).Avg{j,k} = "NaN";
                avglayer(j,K) = nan;
                avglayer(j,K+1) = nan;
                avglayer(j,K+2) = nan;
                FreqMatches.(pt).N1{j,k} = "NaN";
                N1layer(j,K) = nan;
                N1layer(j,K+1) = nan;
                N1layer(j,K+2) = nan;
                FreqMatches.(pt).N2{j,k} = "NaN";
                N2layer(j,K) = nan;
                N2layer(j,K+1) = nan;
                N2layer(j,K+2) = nan;
                FreqMatches.(pt).Rand{j,k} = "NaN";
                Randlayer(j,K) = nan;
                Randlayer(j,K+1) = nan;
                Randlayer(j,K+2) = nan;
                FreqMatches.(pt).Null{j,k} = "NaN";
                Nulllayer(j,K) = nan;
                Nulllayer(j,K+1) = nan;
                Nulllayer(j,K+2) = nan;
                continue;
            end
            a = a+1;
            avg_freqs = round(ActualDomFreqs.Avg{j,k});%in Hz
            n1_freqs = round(ActualDomFreqs.N1{j,k});
            n2_freqs = round(ActualDomFreqs.N2{j,k});
            rand_freqs = round(ActualDomFreqs.Rand{j,k});
            null_freqs = round(ActualDomFreqs.Null{j,k});

            a_length.(pt)(a) = length(avg_freqs);
            n1_length.(pt)(a) = length(n1_freqs);
            n2_length.(pt)(a) = length(n2_freqs);
            r_length.(pt)(a) = length(rand_freqs);
            n_length.(pt)(a) = length(null_freqs);
            
            I_avg = intersect(avg_freqs,bode_freqs);
            if ~isempty(I_avg)
                FreqMatches.(pt).Avg{j,k} = I_avg;
                avglayer(j,K) = I_avg(1);
                if length(I_avg)>1
                    avglayer(j,K+1) = I_avg(2);
                else
                    avglayer(j,K+1) = nan;
                end
                if length(I_avg)>2
                    avglayer(j,K+2)=I_avg(3);
                else
                    avglayer(j,K+2) = nan;
                end
            else
                FreqMatches.(pt).Avg{j,k} = "NaN";
                avglayer(j,K) = nan;
                avglayer(j,K+1) = nan;
                avglayer(j,K+2) = nan;
                M_avg(j,k) = 0.5;
            end
            clearvars I_avg avg_freqs
            
            I_n1 = intersect(n1_freqs,bode_freqs);
            if ~isempty(I_n1)
                FreqMatches.(pt).N1{j,k} = I_n1;
                N1layer(j,K) = I_n1(1);
                if length(I_n1)>1
                    N1layer(j,K+1) = I_n1(2);
                else
                    N1layer(j,K+1) = nan;
                end
                if length(I_n1)>2
                    N1layer(j,K+2)=I_n1(3);
                else
                    N1layer(j,K+2) = nan;
                end
            else
                FreqMatches.(pt).N1{j,k} = "NaN";
                N1layer(j,K) = nan;
                N1layer(j,K+1) = nan;
                N1layer(j,K+2) = nan;
                M_n1(j,k) = 0.5;
            end
            clearvars I_n1 n1_freqs
        
            I_n2 = intersect(n2_freqs,bode_freqs);
            if ~isempty(I_n2)
                FreqMatches.(pt).N2{j,k} = I_n2;
                N2layer(j,K) = I_n2(1);
                if length(I_n2)>1
                    N2layer(j,K+1) = I_n2(2);
                else
                    N2layer(j,K+1) = nan;
                end
                if length(I_n2)>2
                    N2layer(j,K+2)=I_n2(3);
                else
                    N2layer(j,K+2) = nan;
                end
            else
                FreqMatches.(pt).N2{j,k} = "NaN";
                N2layer(j,K) = nan;
                N2layer(j,K+1) = nan;
                N2layer(j,K+2) = nan;
                M_n2(j,k) = 0.5;
            end
            clearvars I_n2 n2_freqs

            I_rand = intersect(rand_freqs,bode_freqs);
            if ~isempty(I_rand)
                FreqMatches.(pt).Rand{j,k} = I_rand;
                Randlayer(j,K) = I_rand(1);
                if length(I_rand)>1
                    Randlayer(j,K+1) = I_rand(2);
                else
                    Randlayer(j,K+1) = nan;
                end
                if length(I_rand)>2
                    Randlayer(j,K+2)=I_rand(3);
                else
                    Randlayer(j,K+2) = nan;
                end
            else
                FreqMatches.(pt).Rand{j,k} = "NaN";
                Randlayer(j,K) = nan;
                Randlayer(j,K+1) = nan;
                Randlayer(j,K+2) = nan;
                M_rand(j,k) = 0.5;
            end
            clearvars I_rand rand_freqs

            I_null = intersect(null_freqs,bode_freqs);
            if ~isempty(I_null)
                FreqMatches.(pt).Null{j,k} = I_null;
                Nulllayer(j,K) = I_null(1);
                if length(I_null)>1
                    Nulllayer(j,K+1) = I_null(2);
                else
                    Nulllayer(j,K+1) = nan;
                end
                if length(I_null)>2
                    Nulllayer(j,K+2)=I_null(3);
                else
                    Nulllayer(j,K+2) = nan;
                end
            else
                FreqMatches.(pt).Null{j,k} = "NaN";
                Nulllayer(j,K) = nan;
                Nulllayer(j,K+1) = nan;
                Nulllayer(j,K+2) = nan;
                M_null(j,k) = 0.5;
            end
            clearvars I_null null_freqs bode_freqs
        end
    end
    Match_avg = find(M_avg==0);
    NonMatch_avg = find(M_avg==0.5);
    CXR.(pt).avg = length(Match_avg)/(length(Match_avg)+length(NonMatch_avg));
    CXR.All.avg(i) = length(Match_avg)/(length(Match_avg)+length(NonMatch_avg));
    Match_n1 = find(M_n1==0);
    NonMatch_n1 = find(M_n1==0.5);
    CXR.(pt).n1 = length(Match_n1)/(length(Match_n1)+length(NonMatch_n1));
    CXR.All.n1(i) = length(Match_n1)/(length(Match_n1)+length(NonMatch_n1));
    Match_n2 = find(M_n2==0);
    NonMatch_n2 = find(M_n2==0.5);
    CXR.(pt).n2 = length(Match_n2)/(length(Match_n2)+length(NonMatch_n2));
    CXR.All.n2(i) = length(Match_n2)/(length(Match_n2)+length(NonMatch_n2));
    Match_rand = find(M_rand==0);
    NonMatch_rand = find(M_rand==0.5);
    CXR.(pt).rand = length(Match_rand)/(length(Match_rand)+length(NonMatch_rand));
    CXR.All.rand(i) = length(Match_rand)/(length(Match_rand)+length(NonMatch_rand));
    Match_null = find(M_null==0);
    NonMatch_null = find(M_null==0.5);
    CXR.(pt).null = length(Match_null)/(length(Match_null)+length(NonMatch_null));
    CXR.All.null(i) = length(Match_null)/(length(Match_null)+length(NonMatch_null));
    clearvars -except FreqMatches i pt_list CustomColormap CXR a_length n1_length n2_length r_length n_length
end
% cd /data/project/NSPMlab/hbriny99/UAB_CCEPdata/CCSR_Resonance_Comparison
% save("CXR","CXR",'-v7.3');
% save("FreqMatches","FreqMatches",'-v7.3');

%% CXRsame figure 
cd('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/CCSR_Resonance_Comparison')
load("SplitCorr.mat")
pt_list = string(fieldnames(SplitCorr));
pt_list_ignore = 2;
pt_list(pt_list_ignore) = [];
figure; 
hold on; 
x1 = ones(1,length(pt_list));
x2 = x1*2;
x3 = x1*3;
x4 = x1*4;
y1 = CXR.All.null*100;
y2 = CXR.All.n1*100; 
y3 = CXR.All.n2*100; 
y4 = CXR.All.avg*100;
boxplot([y1 y2 y3 y4],[x1 x2 x3 x4],"Notch","on","Colors",[0 0 0],"PlotStyle","traditional");
scatter(1,y1,"MarkerEdgeColor",'k','MarkerFaceColor','k'); 
scatter(2,y2,"MarkerEdgeColor",'k','MarkerFaceColor','b'); 
scatter(3,y3,"MarkerEdgeColor",'k','MarkerFaceColor','r'); 
scatter(4,y4,"MarkerEdgeColor",'k','MarkerFaceColor','g');
xlim([0.5 4.5]);
ylim([0 60]);
xticklabels(["Baseline" "N1" "N2" "Average"]);
xlabel("CCSR Curves")
ylabel("Frequency Co-expression (%)")
title("Percentage of Frequency Co-expression")

%% CXRsame vs CXRrandom Figure
cd('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/CCSR_Resonance_Comparison')
load("MatchFraction.mat")
TrueMF = CXR; %CXRsame
clearvars MatchFraction
load("SplitCorr.mat")
pt_list = string(fieldnames(SplitCorr));
pt_list_ignore = [2];
pt_list(pt_list_ignore) = [];

figure; 
hold on; 
%N1
x1 = ones(1,length(pt_list));
x2 = x1*2;
for i = 1:length(pt_list)
    pt = string(pt_list(i));
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt,'/TFMsBode_figures2/Bode-CCSR Freq Comparison'));
    load("RandFreqMatches.mat");
    RandLocMatchFraction = RandFreqMatches;
    y1(i) = TrueMF.(pt).n1;%same
    y2(i) = median(RandLocMatchFraction.(pt).n1);%random
end
xdata = [x1' x2' nan(numel(x1),1)];
ydata = [y1' y2' nan(numel(y1),1)];
subplot(1,3,1);
hold on
for i = 1:length(x1)
    plot([xdata(i,1) xdata(i,2)],[ydata(i,1) ydata(i,2)],'k');
end
scatter(x1,y1,"MarkerEdgeColor",'k','MarkerFaceColor','b');
hold on;
scatter(x2,y2,"MarkerEdgeColor",'k','MarkerFaceColor',[.7 .7 1]);
xlim([0 3]);
ylim([0 1])
w = ranksum(y2,y1);
disp(w);
%N2
x1 = ones(1,length(pt_list));
x2 = x1*2;
for i = 1:length(pt_list)
    pt = string(pt_list(i));
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt,'/TFMsBode_figures2/Bode-CCSR Freq Comparison'));
    load("RandFreqMatches.mat");
    RandLocMatchFraction = RandFreqMatches;
    y1(i) = TrueMF.(pt).n2;%same
    y2(i) = median(RandLocMatchFraction.(pt).n2);%random
end
xdata = [x1' x2' nan(numel(x1),1)];
ydata = [y1' y2' nan(numel(y1),1)];
subplot(1,3,2);
hold on
for i = 1:length(x1)
    plot([xdata(i,1) xdata(i,2)],[ydata(i,1) ydata(i,2)],'k');
end
scatter(x1,y1,"MarkerEdgeColor",'k','MarkerFaceColor','r');
hold on;
scatter(x2,y2,"MarkerEdgeColor",'k','MarkerFaceColor',[1 .7 .7]);
xlim([0 3]);
ylim([0 1])
w = ranksum(y2,y1);
disp(w);
%Average
x1 = ones(1,length(pt_list));
x2 = x1*2;
for i = 1:length(pt_list)
    pt = string(pt_list(i));
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt,'/TFMsBode_figures2/Bode-CCSR Freq Comparison'));
    load("RandFreqMatches.mat");
    RandLocMatchFraction = RandFreqMatches;
    y1(i) = TrueMF.(pt).avg;%same
    y2(i) = median(RandLocMatchFraction.(pt).avg);%random
end
xdata = [x1' x2' nan(numel(x1),1)];
ydata = [y1' y2' nan(numel(y1),1)];
subplot(1,3,3);
hold on
for i = 1:length(x1)
    plot([xdata(i,1) xdata(i,2)],[ydata(i,1) ydata(i,2)],'k');
end
scatter(x1,y1,"MarkerEdgeColor",'k','MarkerFaceColor','g');
hold on;
scatter(x2,y2,"MarkerEdgeColor",'k','MarkerFaceColor',[.7 1 .7]);
xlim([0 3]);
ylim([0 1])
w = ranksum(y2,y1);
disp(w);
