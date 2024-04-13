% just a first check that sampling from an F distribution produces
% something that resembles typical alpha power values,
% non-essential for the code below
rng(123)
num_samples = 10000; % Number of random samples
random_values1 = random('F', 10, 20, num_samples, 1);
random_values2 = random('F', 10, 20, num_samples, 1)*1.2;

% % Plot a histogram of the generated random values
% histogram(random_values1,100, 'Normalization', 'probability');
% hold on
% histogram(random_values2,100, 'Normalization', 'probability');
% xlabel('Random Values');
% ylabel('Probability');
% title('Histogram of Random Values from F Distribution');

%% use this for alpha power
clearvars

permNo = 100; % number of "studies"

rng(123) % to produce replicable results
num_samples = 20; % Number of participants per study

% all the effect sizes with different measures - pre-allocate
t_sim_ES = zeros(permNo,1);
t_AMI_ES = zeros(permNo,1);
F_ES = zeros(permNo,1);
t_sim_log_ES = zeros(permNo,1);
t_AMI_log_ES = zeros(permNo,1);
t_sim_bl_ES = zeros(permNo,1);
t_AMI_bl_ES = zeros(permNo,1);
t_raw_ALI_ES = zeros(permNo,1);
t_norm_ALI_ES = zeros(permNo,1);
t_log_ALI_ES = zeros(permNo,1);
t_bl_ALI_ES  = zeros(permNo,1);

for i_p = 1:permNo

% generate power distrib
pow_F=random('F', 10, 20, num_samples, 1)+.5; % power distribution
% assume same effect L&R
mean_eff=0.1;
sd_eff=2*mean_eff;
att_eff=random("normal",mean_eff,sd_eff,num_samples,1);
% in this case, the expected effect (Cohen's d) is given by mean_eff/sd_eff

HemiL     = pow_F+(rand(num_samples,1)-.5)*.2; % power distribution + noise
HemiLCueR = HemiL-att_eff/2; % att vis field = lower power
HemiLCueL = HemiL+att_eff/2; % ign vis field = higher power
HemiR     = pow_F+(rand(num_samples,1)-.5)*.2;
HemiRCueL = HemiR-att_eff/2;
HemiRCueR = HemiR+att_eff/2;

% HemiL & HemiR can be used as "baseline" or "neutral cue condition"

alphaPOW=[HemiLCueR,HemiLCueL,HemiRCueL,HemiRCueR];
% figure,boxplot(alphaPOW)

%% calculate hemi-specific effects
hemiSimpleEff=[HemiLCueL-HemiLCueR,HemiRCueR-HemiRCueL];
% normalised
hemiAMI      =[(HemiLCueL-HemiLCueR)./(HemiLCueL+HemiLCueR),...
               (HemiRCueR-HemiRCueL)./(HemiRCueR+HemiRCueL)];

% meanEffectSize(hemiSimpleEff(:,1),Mean=0,Effect="cohen")
% meanEffectSize(HemiLCueL,HemiLCueR,Paired=true,Effect="cohen")
% %meanEffectSize(hemiSimpleEff(:,2),0,Effect="cohen")
% meanEffectSize(hemiAMI(:,1),0,Effect="cohen")
% %meanEffectSize(hemiAMI(:,2),0,Effect="cohen")

% for before R2021b because I do not have meanEffectSize
t_sim_ES(i_p) = mean(hemiSimpleEff(:,1))/std(hemiSimpleEff(:,1));
t_AMI_ES(i_p) = mean(hemiAMI(:,1))/std(hemiAMI(:,1));

% observations
% 1) differs between measures - normalised smaller(?)
% 2) noise inflation by normalisation?
% (could run simulation to check/quantify)

%% eff size for cue x hemi interaction

% by ANOVA
depvar=     [HemiRCueL;HemiRCueR;HemiLCueL;HemiLCueR];
% design
subvar=repmat((1:num_samples).',4,1);
hemi  =[repmat(ones(num_samples,1),2,1);2*repmat(ones(num_samples,1),2,1)];
cue   =[ones(num_samples,1);2*ones(num_samples,1);ones(num_samples,1);2*ones(num_samples,1)];

[p,tbl,stats]=anovan(depvar,{subvar,hemi,cue},...
    "Model","interaction",...
    "varnames",["sub","hemi","cue"],...
    "random",1,...
    "display",'off');

Fval = tbl{7,6};
dof1 = tbl{7,3};
dof2 = tbl{8,3};
EST = Fval.*dof1./(Fval.*dof1+dof2);
F_ES(i_p) = sqrt((num_samples-1)/num_samples * EST/(1-EST)); % convert to Cohen's d
clear Fval dof1 dof2

% % extract eff size by ssq ratio from table ...
% % for an ANOVA:
% EScorr=Fval.*dof1./(Fval.*dof1+dof2);
% % for a t-test (dof1 = 1):
% EScorr=tval^2./(tval^2+dof2);

% % by fitting a repeated measures model
% datatab=table(subvar,hemi,cue,depvar,'VariableNames',{'subj','hemi','cue','pow'});
% 
% lme = fitlme(datatab,'pow~hemi+cue+hemi*cue+(1|subj)');

%% same but for log vals
% compute all kinds of derivatives
logPOW=10*log10(alphaPOW);
% figure, boxplot(logPOW)

hemiLogSimpleEff=[logPOW(:,2)-logPOW(:,1),logPOW(:,4)-logPOW(:,3)];

% normalised - not sth that should be done, but maybe shows up in the lit
hemiLogAMI   =[(logPOW(:,2)-logPOW(:,1))./(logPOW(:,2)+logPOW(:,1)),...
               (logPOW(:,4)-logPOW(:,3))./(logPOW(:,4)+logPOW(:,3))];
 
% meanEffectSize(hemiLogSimpleEff(:,1),0,Effect="cohen")
% meanEffectSize(hemiLogSimpleEff(:,2),0,Effect="cohen")
% meanEffectSize(hemiLogAMI(:,1),0,Effect="cohen")
% meanEffectSize(hemiLogAMI(:,2),0,Effect="cohen")

% for before R2021b because I do not have meanEffectSize
t_sim_log_ES(i_p) = mean(hemiLogSimpleEff(:,1))/std(hemiLogSimpleEff(:,1));
t_AMI_log_ES(i_p) = mean(hemiLogAMI(:,1))/std(hemiLogAMI(:,1));

%% same but for baseline-normalised values

blPOW=[(HemiLCueR-HemiL)./HemiL,...
       (HemiLCueL-HemiL)./HemiL,...
       (HemiRCueL-HemiR)./HemiR,...
       (HemiRCueR-HemiR)./HemiR];

hemiBlSimpleEff=[blPOW(:,2)-blPOW(:,1),blPOW(:,4)-blPOW(:,3)];

% the below doesn't make sense because it's a "double normalisation"
% producing inf values, only kept for symmetry
hemiBlAMI   =[(blPOW(:,2)-blPOW(:,1))./(blPOW(:,2)+blPOW(:,1)),...
               (blPOW(:,4)-blPOW(:,3))./(blPOW(:,4)+blPOW(:,3))];

% for before R2021b because I do not have meanEffectSize
t_sim_bl_ES(i_p) = mean(hemiBlSimpleEff(:,1))/std(hemiBlSimpleEff(:,1));
t_AMI_bl_ES(i_p) = mean(hemiBlAMI(:,1))/std(hemiBlAMI(:,1));

%% compute ALI

% should average L & R hemi data here, but will only look at:
% L=contra & R=ipsi for now to mirror the selection of L hemi data above
% for effect size calculation
% (alphaPOW=[HemiLCueR,HemiLCueL,HemiRCueL,HemiRCueR];)
% (rawALI=(HemiLCueR-HemiRCueR)./(HemiLCueR+HemiRCueR);)

rawALI=(alphaPOW(:,1)-alphaPOW(:,4));
normALI=(alphaPOW(:,1)-alphaPOW(:,4))./(alphaPOW(:,1)+alphaPOW(:,4));
logALI=(logPOW(:,1)-logPOW(:,4));
blALI =(blPOW(:,1)-blPOW(:,4));

% make effect sizes positive
t_raw_ALI_ES(i_p) = abs(mean(rawALI(:,1))/std(rawALI(:,1)));
t_norm_ALI_ES(i_p) = abs(mean(normALI(:,1))/std(normALI(:,1)));
t_log_ALI_ES(i_p) = abs(mean(logALI(:,1))/std(logALI(:,1)));
t_bl_ALI_ES (i_p) = abs(mean(blALI(:,1))/std(blALI(:,1)));

end

all_ES = [t_sim_ES,t_AMI_ES,F_ES,...
          t_sim_log_ES,t_AMI_log_ES,...
          t_sim_bl_ES,t_AMI_bl_ES,...
          t_raw_ALI_ES,t_norm_ALI_ES,t_log_ALI_ES,t_bl_ALI_ES];
figure, boxplot(all_ES)
estimLabel={'RawSimple';'RawAMI';'RawCueHemiInter';...
            'LogSimple';'LogAMI';...
            'BaselineSimple';'BaselineAMI';...
            'ALI-Raw';'ALI-norm';'ALI-Log';'ALI-Baseline'};
set(gca,'xticklabel',estimLabel);

%% alpha lateralisation

% use raw values, (contra-ipsi/contra+ipsi)
