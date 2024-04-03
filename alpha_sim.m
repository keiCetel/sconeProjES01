% just a first check that sampling from an F distribution produces
% something that resembles typical alpha power values
% non-essential for the code below
rng(123)
num_samples = 10000; % Number of random samples
random_values1 = random('F', 10, 20, num_samples, 1);
random_values2 = random('F', 10, 20, num_samples, 1)*1.2;

% Plot a histogram of the generated random values
histogram(random_values1,100, 'Normalization', 'probability');
hold on
histogram(random_values2,100, 'Normalization', 'probability');
xlabel('Random Values');
ylabel('Probability');
title('Histogram of Random Values from F Distribution');

%% use this for alpha power
clearvars
rng(123)
num_samples = 100; % Number of random samples
% generate power distrib
pow_F=random('F', 10, 20, num_samples, 1)+.5; % power distribution
% assume same effect L&R
mean_eff=0.1;
sd_eff=2*mean_eff;
att_eff=random("normal",mean_eff,sd_eff,num_samples,1);
% the expected effect (Cohen's d) is then mean_eff/sd_eff

HemiLCueR = pow_F+(rand(num_samples,1)-.5)*.2; % power distribution + noise
HemiLCueL = HemiLCueR+att_eff; % ign vis field = higher power
HemiRCueL = pow_F+(rand(num_samples,1)-.5)*.2;
HemiRCueR = HemiRCueL+att_eff; % ign vis field = higher power

alphaPOW=[HemiLCueR,HemiLCueL,HemiRCueL,HemiRCueR];
figure,boxplot(alphaPOW)

%% calculate hemi-specific effects
hemiSimpleEff=[HemiLCueL-HemiLCueR,HemiRCueR-HemiRCueL];
% normalised
hemiAMI      =[(HemiLCueL-HemiLCueR)./(HemiLCueL+HemiLCueR),...
               (HemiRCueR-HemiRCueL)./(HemiRCueR+HemiRCueL)];

meanEffectSize(hemiSimpleEff(:,1),0,Effect="cohen")
meanEffectSize(hemiSimpleEff(:,2),0,Effect="cohen")
meanEffectSize(hemiAMI(:,1),0,Effect="cohen")
meanEffectSize(hemiAMI(:,2),0,Effect="cohen")
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
    "random",1);

% extract eff size by ssq ratio from table ...
% for an ANOVA:
EScorr=Fval.*dof1./(Fval.*dof1+dof2);
% for a t-test (dof1 = 1):
EScorr=tval^2./(tval^2+dof2);

% by fitting a repeated measures model
datatab=table(subvar,hemi,cue,depvar,'VariableNames',{'subj','hemi','cue','pow'});

lme = fitlme(datatab,'pow~hemi+cue+hemi*cue+(1|subj)');


%% same but for log vals
% compute all kinds of derivatives
logPOW=10*log10(alphaPOW);
figure, boxplot(logPOW)

hemiLogSimpleEff=[logPOW(:,2)-logPOW(:,1),logPOW(:,4)-logPOW(:,3)];
% normalised:
hemiLogAMI   =[(logPOW(:,2)-logPOW(:,1))./(logPOW(:,2)+logPOW(:,1)),...
                  (logPOW(:,4)-logPOW(:,3))./(logPOW(:,4)+logPOW(:,3))];

meanEffectSize(hemiLogSimpleEff(:,1),0,Effect="cohen")
meanEffectSize(hemiLogSimpleEff(:,2),0,Effect="cohen")
meanEffectSize(hemiLogAMI(:,1),0,Effect="cohen")
meanEffectSize(hemiLogAMI(:,2),0,Effect="cohen")
