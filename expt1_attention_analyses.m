% expt1_attention_analyses.m
%
% This script completes analyses on the data from Experiment 1 of 
% deBettencourt et al PB&R using the sustained attention phase
% 
% Instructions to execute script:
% - place above the folder containing the data ("data")
% - it also requires the matlab exchange function shaded error bar for some plots  
%
% Author: Megan deBettencourt, debetten@uchicago.edu

%all of the subjects included in the experiment 
%all_subj = [51 52 53 54 55 56 57 58 59 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83];
%Note: I started counting at 51... because after piloting things like experiment length I rounded to 50
%      I know it looks weird like I picked the most random number! It's not I swear
%excludeID_list = [4,10,12,21,25]; %@ these are the IDs to be excluded
all_subj = [1,2,3,5,6,7,8,9,11,13,14,15,16,17,18,19,20,22,23,24,26,27,28,29,30,31,32,33,34,35,36,37];
nsubj = numel(all_subj);


%%

for isubj = 1:nsubj
    
    subjectNum = all_subj(isubj);
    
    %load data
    fn = deblank(ls(['data/' num2str(subjectNum) '/attndata_*']));
    %fn = ['data/', num2str(subjectNum), '/', fn];
    load(fn);
    
    %number of trials in the attention phase
    ntrials_attn(isubj) = numel(attnData.trial);
    
    % find outdoor and indoor trials
    ind_outdoor = find(attnData.categs==1);
    ind_indoor = find(attnData.categs==2);
    
    % indices of frequent and rare images 
    if numel(ind_outdoor)>numel(ind_indoor)
        ind_freq{isubj} = ind_outdoor;
        ind_rare{isubj} = ind_indoor;
    else
        ind_freq{isubj} = ind_indoor;
        ind_rare{isubj} = ind_outdoor;
    end
    
    % which category (outdoor or indoor) was frequent and rare 
    if median(attnData.categs)==1
        freq_categ(isubj) = 1;
        rare_categ(isubj) = 2;
    else
        freq_categ(isubj) = 2;
        rare_categ(isubj) = 1;
    end
    
    % calculate hit rates/false alarm rates
    hit_rate(isubj) = sum(attnData.accs(ind_freq{isubj}))/numel(ind_freq{isubj});
    miss_rate(isubj) = sum(attnData.accs(ind_freq{isubj})==0)/numel(ind_freq{isubj});
    corrrej_rate(isubj) = sum(attnData.accs(ind_rare{isubj})==1)/numel(ind_rare{isubj});
    falsealarm_rate(isubj) = sum(attnData.accs(ind_rare{isubj})==0)/numel(ind_rare{isubj});
    
    % calculate a prime
    aprime(isubj) = .5 + ((hit_rate(isubj)-falsealarm_rate(isubj))*(1+hit_rate(isubj)-falsealarm_rate(isubj)))/(4*hit_rate(isubj)*(1-falsealarm_rate(isubj)));
    
    % overall accuracy to rare images
    rare_acc(isubj,:) = attnData.accs(ind_rare{isubj});
    
    %calculate a linear fit to the RTs
    attnData.linfit = polyfit(attnData.trial(~isnan(attnData.rts)),attnData.rts(~isnan(attnData.rts)),1);
    attnData.rts_est = attnData.linfit(1).*attnData.trial+attnData.linfit(2);
    attnData.rts_resid = attnData.rts-attnData.rts_est; %regress out linear RT trend across trials
    
    % calculate the RTs surrounding the rare category
    for iLure =1:(numel(ind_rare{isubj}))
        tpts = ind_rare{isubj}(iLure)+(-6:6);
        %exclude the rare trials that were in the first 6 or last 6 images
        if (ind_rare{isubj}(iLure)>6) && (ind_rare{isubj}(iLure)<(attnData.trialsPerRun-6))
            %find the surrounding timepoints
            
            %if you want to use RTs - commented out 
            %rts_aroundraretrials(iLure,:) = blockData.rts(tpts);
            
            %if you want to use RT residuals, after regressing linear trend
            rts_aroundraretrials(iLure,:) = attnData.rts_resid(tpts);
        else
            rts_aroundraretrials(iLure,:) = NaN(1,numel(tpts));
        end
    end
    
    subj_rare_acc{subjectNum} = rare_acc(isubj,:);
    subj_rts_aroundraretrials{subjectNum} = rts_aroundraretrials;
    subj_rts_aroundraretrials_correct(isubj,:) = nanmean(rts_aroundraretrials(rare_acc(isubj,:)==1,:));
    subj_rts_aroundraretrials_error(isubj,:) = nanmean(rts_aroundraretrials(rare_acc(isubj,:)==0,:));
    subj_aprime(isubj) = aprime(isubj);
    subj_rt_mean(isubj) = nanmean(attnData.rts);
    subj_rt_std(isubj) = nanstd(attnData.rts);
    
    % quartile based analyses
    for iquartile = 1:4
        t = 125*(iquartile-1)+(1:125); % probably to-do: replace magic numbers with more generalisable code that queries trials per run from the attention data
        temp_accs = attnData.accs(t);
        temp_rts = attnData.rts(t);
        temp_rare = find(attnData.categs(t)==rare_categ(isubj));
        temp_freq = find(attnData.categs(t)==freq_categ(isubj));
        
        quartile_commission_errors(isubj,iquartile) = mean(temp_accs(temp_rare)==0)*100;
        quartile_ommission_errors(isubj,iquartile) = mean(temp_accs(temp_freq)==0)*100;
        quartile_rt(isubj,iquartile) = nanmean(temp_rts)*1000;
        quartile_rtvar(isubj,iquartile) = nanstd(temp_rts)*1000;
    end
    
    % sliding window analyses
    for islidingwindow = 1:(attnData.trial(end)-119)
        t = islidingwindow+(0:119);
        ifreq{isubj}{islidingwindow} = find(attnData.categs(t)==freq_categ(isubj));
        irare{isubj}{islidingwindow} = find(attnData.categs(t)==rare_categ(isubj));
        temp_accs = attnData.accs(t);
        temp_rts = attnData.rts(t);
        lapse_rate(isubj,islidingwindow) = 1-mean(temp_accs(irare{isubj}{islidingwindow}));
        coeff_var(isubj,islidingwindow) = nanstd(temp_rts(temp_accs==1))/nanmean(temp_rts(temp_accs==1));
        omissionerrs_overtime(isubj,islidingwindow) = 1-mean(temp_accs(ifreq{isubj}{islidingwindow}));
        rts_overtime_correct(isubj,islidingwindow) = nanmean(temp_rts(temp_accs==1));
        rts_overtime(isubj,islidingwindow) = nanmean(temp_rts);
    end
end

%% Figure 1; Miss rate vs false alarm rate

figure;
hold on;
h = bar(1,mean(miss_rate),'g');
h = bar(2,mean(falsealarm_rate),'r');
errorbar([1 2],[mean(miss_rate) mean(falsealarm_rate)],[std(miss_rate) std(falsealarm_rate)]/sqrt(nsubj),'k+');
set(gca,'xlim',[.5 2.5],'xtick',[1,2],'xticklabel',{'miss rate','false alarm rate'},'ylim',[0 .4],'ytick',[0 .1 .2 .3 .4]);
title('Frequent trials')

%% Figure 2

figure;
subplot(1,3,1);
hold on;
h = bar(1,mean(hit_rate),'g');
h = bar(2,mean(miss_rate),'r');
errorbar([1 2],[mean(hit_rate) mean(miss_rate)],[std(hit_rate) std(miss_rate)],'k+');
set(gca,'xlim',[.5 2.5],'xtick',[],'ylim',[0 1],'ytick',[0 .25 .5 .75 1]);
title('Frequent trials')

subplot(1,3,2);
hold on;
bar(1,mean(corrrej_rate),'g');
bar(2,mean(falsealarm_rate),'r');
errorbar([1 2],[mean(corrrej_rate) mean(falsealarm_rate)],[std(corrrej_rate) std(falsealarm_rate)],'k+');
set(gca,'xlim',[.5 2.5],'xtick',[],'ylim',[0 1],'ytick',[0 .25 .5 .75 1]);
title('Rare trials')

subplot(1,3,3);
hold on;
h = bar(mean(aprime),'w');
errorbar(1,mean(aprime),std(aprime),'k+');
plot([0 3],[.5 .5],'--','color',[.5 .5 .5]);
set(gca,'xlim',[.5 1.5],'xtick',[],'ylim',[0 1],'ytick',[0 .25 .5 .75 1]);
title('A prime')

%% Figure 3; all subjects RTs around lure

figure;
hold on;
errorbar((1:13)-7,mean(subj_rts_aroundraretrials_correct),std(subj_rts_aroundraretrials_correct)/sqrt(nsubj),'g');
errorbar((1:13)-7,mean(subj_rts_aroundraretrials_error),std(subj_rts_aroundraretrials_error)/sqrt(nsubj),'r');
set(gca,'xlim',[-7,7],'xtick',(1:13)-7);%,'xticklabel',{'-3','-2','-1','0','1','2','3'})
%set(gca,'ytick',[.3 .4 .5 .6 .7])
ylabel('RT resid [s]');
xlabel('Trial from rare subcategory');
title(['Average RT around rare trial']);

%% Figure 4; individual subjects RTs around lure

figure;
for isubj = 1:nsubj
    subjectNum = all_subj(isubj);
    subplot(4,8,isubj);
    hold on;
    errorbar(-6:6,nanmean(subj_rts_aroundraretrials{subjectNum}(subj_rare_acc{subjectNum}==1,:)),nanstd(subj_rts_aroundraretrials{subjectNum}(subj_rare_acc{subjectNum}==1,:))/sqrt(sum(subj_rare_acc{subjectNum}==1)),'g')
    errorbar(-6:6,nanmean(subj_rts_aroundraretrials{subjectNum}(subj_rare_acc{subjectNum}==0,:)),nanstd(subj_rts_aroundraretrials{subjectNum}(subj_rare_acc{subjectNum}==0,:))/sqrt(sum(subj_rare_acc{subjectNum}==0)),'r')
    set(gca,'xlim',[-7,7],'ylim',[-.2 .2],'ytick',[],'xtick',[]);
    title([num2str(isubj)],'fontsize',12);
end


%% Figure 5 (figure 3 from rosenberg ap&p)

figure;

%figure 3a from rosenberg ap&p
subplot(2,2,1);
hold on;
shadedErrorBar(1:4,mean(quartile_commission_errors),std(quartile_commission_errors)/sqrt(nsubj),'lineprops','k');
plot(1:4,mean(quartile_commission_errors),'k');
set(gca,'ylim',[25,40],'ytick',[25 30 35 40]);
ylabel('Commission errors')
xlabel('Quartile')

%figure 3b from rosenberg ap&p
subplot(2,2,2);
hold on;
shadedErrorBar(1:4,mean(quartile_rtvar),std(quartile_rtvar)/sqrt(nsubj),'lineprops','k');
plot(1:4,mean(quartile_rtvar),'k');
ylabel('RT variability')
xlabel('Quartile')

%figure 3c from rosenberg ap&p
subplot(2,2,3);
hold on;
shadedErrorBar(1:4,mean(quartile_ommission_errors),std(quartile_ommission_errors)/sqrt(nsubj),'lineprops','k');
plot(1:4,mean(quartile_ommission_errors),'k');
set(gca,'ylim',[0,5],'ytick',[1 2 3 4 5]);
ylabel('Omission errors')
xlabel('Quartile')

%figure 3d from rosenberg ap&p
subplot(2,2,4);
hold on;
shadedErrorBar(1:4,mean(quartile_rt),std(quartile_rt)/sqrt(nsubj),'lineprops','k');
plot(1:4,mean(quartile_rt),'k');
set(gca,'xlim',[.5,4.5],'xtick',[1 2 3 4]);
ylabel('RT')
xlabel('Quartile')

%% Figure 6 (figure 1b from esterman cer cor)

%lapse rate vs coefficient of variation (RT SD/mean)

figure;
hold on;
scatter(falsealarm_rate,subj_rt_std./subj_rt_mean,200,[.5 .5 .5],'filled');
xlabel('Lapse Rate');
ylabel('Coefficient of Variation (RT SD/mean)');
[r,p] = corr(falsealarm_rate',transpose(subj_rt_std./subj_rt_mean))
P = polyfit(falsealarm_rate,subj_rt_std./subj_rt_mean,1);
plot([min(falsealarm_rate) max(falsealarm_rate)],P(1)*[min(falsealarm_rate) max(falsealarm_rate)]+P(2),'k--');
set(gca,'xlim',[0 .75],'xtick',[0 .25 .5 .75],'ylim',[0 .35],'ytick',[0 .1 .2 .3])

%% Figure 7 (figure 1b but with means and not SD/mean)

figure;
hold on;
scatter(falsealarm_rate,subj_rt_mean,200,[.5 .5 .5],'filled');
xlabel('Lapse Rate');
ylabel('Mean RT (s)');
[r,p] = corr(falsealarm_rate',subj_rt_mean')
P = polyfit(falsealarm_rate,subj_rt_mean,1);
plot([min(falsealarm_rate) max(falsealarm_rate)],P(1)*[min(falsealarm_rate) max(falsealarm_rate)]+P(2),'k--');
set(gca,'xlim',[0 .75],'xtick',[0 .25 .5 .75],'ylim',[.3 .6],'ytick',[.3 .4 .5 .6])


%% Figure 8 (figures 1c-f from esterman cer cor)

%figure 1c from esterman cer cor: sliding time window vs lapse rate
figure;
set(gcf,'Position',[0 0 1600 600]);
subplot(1,4,1);
hold on;
shadedErrorBar(1:381,nanmean(lapse_rate(:,1:381)),nanstd(lapse_rate(:,1:381))./sqrt(sum(~isnan(lapse_rate(:,1:381)))));
xlabel('Sliding time window (2min)');
ylabel('Lapse Rate');
set(gca,'ylim',[.2 .45],'ytick',[.2 .25 .3 .35 .4 .45])
set(gca,'xlim',[1 381],'xtick',[1:60:381],'xticklabel',{'0-2','1-3','2-4','3-5','4-6','5-7','6-8'})

%figure 1d from esterman cer cor
subplot(1,4,2);
hold on;
shadedErrorBar(1:381,nanmean(coeff_var(:,1:381)),nanstd(coeff_var(:,1:381))./sqrt(sum(~isnan(coeff_var(:,1:381)))));
xlabel('Sliding time window');
ylabel('Coefficient of variation');
set(gca,'ylim',[.2 .3],'ytick',[.2 .225 .25 .275 .3])
set(gca,'xlim',[1 381],'xtick',[1:60:381],'xticklabel',{'0-2','1-3','2-4','3-5','4-6','5-7','6-8'})

%figure 1e from esterman cer cor
subplot(1,4,3);
hold on;
shadedErrorBar(1:381,nanmean(omissionerrs_overtime(:,1:381)),nanstd(omissionerrs_overtime(:,1:381))./sqrt(sum(~isnan(omissionerrs_overtime(:,1:381)))));
xlabel('Sliding time window');
ylabel('Omission error rate');
set(gca,'xlim',[1 381],'xtick',[1:60:381],'xticklabel',{'0-2','1-3','2-4','3-5','4-6','5-7','6-8'})

%figure 1f from esterman cer cor
subplot(1,4,4);
hold on;

shadedErrorBar(1:381,nanmean(rts_overtime(:,1:381)),nanstd(rts_overtime(:,1:381))./sqrt(nsubj));
xlabel('Sliding time window');
ylabel('RTs');
set(gca,'xlim',[1 381],'xtick',[1:60:381],'xticklabel',{'0-2','1-3','2-4','3-5','4-6','5-7','6-8'})
set(gca,'ylim',[.45 .55]);




