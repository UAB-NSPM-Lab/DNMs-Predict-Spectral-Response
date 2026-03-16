clearvars -except CCEPs
cd /data/project/NSPMlab/hbriny99/UAB_CCEPdata/CCSR_Resonance_Comparison
load("CCEPs.mat");
cd('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output');
pt_list = string(fieldnames(CCEPs));
for i = 1:length(pt_list)%each pt
    pt = pt_list(i);
    cd(append('/data/project/NSPMlab/hbriny99/UAB_CCEPdata/Output/',pt,'/TFMsBode_figures2/Bode-CCSR Freq Comparison'))
    load("CorrAvg.mat")
    load("CorrN1.mat")
    load("CorrN2.mat")
    load("CorrNrand.mat")
    load("CorrNull.mat")
    load("meanRandCorrAvg.mat")
    load("meanRandCorrN1.mat")
    load("meanRandCorrN2.mat")
    load("meanRandCorrNrand.mat")
    index = CCEPs.(pt).index;
    m = nanmean(index);
    f = find(isnan(m));
    if ~isempty(f)
        index(:,f) = [];
    end
    RccepCorrAvg = CorrAvg(index==1);
    NRccepCorrAvg = CorrAvg(index==0);
    SplitCorr.(pt).Avg.R = RccepCorrAvg;
    SplitCorr.(pt).Avg.NR = NRccepCorrAvg;

    RccepCorrN1 = CorrN1(index==1);
    NRccepCorrN1 = CorrN1(index==0);
    SplitCorr.(pt).N1.R = RccepCorrN1;
    SplitCorr.(pt).N1.NR = NRccepCorrN1;

    RccepCorrN2 = CorrN2(index==1);
    NRccepCorrN2 = CorrN2(index==0);
    SplitCorr.(pt).N2.R = RccepCorrN2;
    SplitCorr.(pt).N2.NR = NRccepCorrN2;

    RccepCorrNrand = CorrNrand(index==1);
    NRccepCorrNrand = CorrNrand(index==0);
    SplitCorr.(pt).Rand.R = RccepCorrNrand;
    SplitCorr.(pt).Rand.NR = NRccepCorrNrand;

    RccepCorrNull = CorrNull(index==1);
    NRccepCorrNull = CorrNull(index==0);
    SplitCorr.(pt).Null.R = RccepCorrNull;
    SplitCorr.(pt).Null.NR = NRccepCorrNull;

    RccepRandCorrAvg = meanRandCorrAvg(index==1);
    NRccepRandCorrAvg = meanRandCorrAvg(index==0);
    SplitCorr.(pt).RandAvg.R = RccepRandCorrAvg;
    SplitCorr.(pt).RandAvg.NR = NRccepRandCorrAvg;

    RccepRandCorrN1 = meanRandCorrN1(index==1);
    NRccepRandCorrN1 = meanRandCorrN1(index==0);
    SplitCorr.(pt).RandN1.R = RccepRandCorrN1;
    SplitCorr.(pt).RandN1.NR = NRccepRandCorrN1;

    RccepRandCorrN2 = meanRandCorrN2(index==1);
    NRccepRandCorrN2 = meanRandCorrN2(index==0);
    SplitCorr.(pt).RandN2.R = RccepRandCorrN2;
    SplitCorr.(pt).RandN2.NR = NRccepRandCorrN2;

    RccepRandCorrNrand = meanRandCorrNrand(index==1);
    NRccepRandCorrNrand = meanRandCorrNrand(index==0);
    SplitCorr.(pt).RandRand.R = RccepRandCorrNrand;
    SplitCorr.(pt).RandRand.NR = NRccepRandCorrNrand;

    SplitCorr.Medians.Avg.R(i) = nanmedian(RccepCorrAvg);
    SplitCorr.Medians.Avg.NR(i) = nanmedian(NRccepCorrAvg);

    SplitCorr.Medians.N1.R(i) = nanmedian(RccepCorrN1);
    SplitCorr.Medians.N1.NR(i) = nanmedian(NRccepCorrN1);

    SplitCorr.Medians.N2.R(i) = nanmedian(RccepCorrN2);
    SplitCorr.Medians.N2.NR(i) = nanmedian(NRccepCorrN2);

    SplitCorr.Medians.Rand.R(i) = nanmedian(RccepCorrNrand);
    SplitCorr.Medians.Rand.NR(i) = nanmedian(NRccepCorrNrand);

    SplitCorr.Medians.Null.R(i) = nanmedian(RccepCorrNull);
    SplitCorr.Medians.Null.NR(i) = nanmedian(NRccepCorrNull);

    SplitCorr.Medians.RandAvg.R(i) = nanmedian(RccepRandCorrAvg);
    SplitCorr.Medians.RandAvg.NR(i) = nanmedian(NRccepRandCorrAvg);

    SplitCorr.Medians.RandN1.R(i) = nanmedian(RccepRandCorrN1);
    SplitCorr.Medians.RandN1.NR(i) = nanmedian(NRccepRandCorrN1);

    SplitCorr.Medians.RandN2.R(i) = nanmedian(RccepRandCorrN2);
    SplitCorr.Medians.RandN2.NR(i) = nanmedian(NRccepRandCorrN2);

    SplitCorr.Medians.RandRand.R(i) = nanmedian(RccepRandCorrNrand);
    SplitCorr.Medians.RandRand.NR(i) = nanmedian(NRccepRandCorrNrand);
end
z = zeros(i,1);
o = ones(i,1);
for j = 1:i
    t(j) = 2; 
end
t = t';

%% RCCEP
pt_list = string(fieldnames(SplitCorr));
pt_list(2) = [];
figure; hold on;
x1 = ones(1,length(pt_list));
x2 = x1*2;
x3 = x1*3;
x4 = x1*4; 
y1 = SplitCorr.Medians.Null.R;
y2 = SplitCorr.Medians.N1.R; 
y3 = SplitCorr.Medians.N2.R; 
y4 = SplitCorr.Medians.Avg.R;
boxplot([y1 y2 y3 y4],[x1 x2 x3 x4],"Notch","on","Colors",[0 0 0],"PlotStyle","traditional");
scatter(x1,y1,"MarkerEdgeColor",'k','MarkerFaceColor','k')
hold on
scatter(x2,y2,"MarkerEdgeColor",'k','MarkerFaceColor','b')
scatter(x3,y3,"MarkerEdgeColor",'k','MarkerFaceColor','r')
scatter(x4,y4,"MarkerEdgeColor",'k','MarkerFaceColor','g')
xlim([0.5 4.5]);
ylim([-1 1])
xticklabels(["Baseline","N1","N2","Average"])    
title("Delayed Responsivity Better Captured by Models")
ylabel("Median Correlation")
xlabel("Responsive CCSR Curves")

%% NRCCEPs/CCSRs for Supplemental Figure
pt_list = string(fieldnames(SplitCorr));
pt_list(2) = [];
figure; hold on;
x1 = ones(1,length(pt_list));
x2 = x1*2;
x3 = x1*3;
x4 = x1*4; 
y1 = SplitCorr.Medians.Null.NR; 
y2 = SplitCorr.Medians.N1.NR;
y3 = SplitCorr.Medians.N2.NR;
y4 = SplitCorr.Medians.Avg.NR;
boxplot([y1 y2 y3 y4],[x1 x2 x3 x4],"Notch","on","Colors",[0 0 0],"PlotStyle","traditional");
scatter(x1,y1,"MarkerEdgeColor",'k','MarkerFaceColor','k')
hold on
scatter(x2,y2,"MarkerEdgeColor",'k','MarkerFaceColor','b')
scatter(x3,y3,"MarkerEdgeColor",'k','MarkerFaceColor','r')
scatter(x4,y4,"MarkerEdgeColor",'k','MarkerFaceColor','g')
xlim([0.5 4.5]);
ylim([-0.2 1])
xticklabels(["Baseline","N1","N2","Average"])    
ylabel("Median Correlation")
xlabel("Non-responsive CCSR Curves")
