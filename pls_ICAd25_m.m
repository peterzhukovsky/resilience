clear
load('D:\Canada_2020\OASIS\data\fMRI\clinical_and_demographic\matched_vars.mat');
subjects=fmri_finalA.Subs;

%amyloid_final=readtable('amyA.csv'); 
clinical_final.Age=clinical_final.ageAtEntry+amyloid_final.time_since_start/365;
clinical_final.delta_fmri=fmri_finalA.Ses-amyloid_final.time_since_start;

%clinical_final.HYPERTEN(clinical_final.test==0)=NaN;clinical_final.HYPERCHO(clinical_final.test==0)=NaN;clinical_final.bmi=(clinical_final.weight*0.453592)./((clinical_final.height*2.54/100).^2);
cognitive_final.MEMUNITS(cognitive_final.Animal==0)=NaN;
                                                full_icad25=full_icad25A;
                                                part_icad25=part_icad25A;
                                                prop_signal=fmri_finalA.prop_signal ; mean_fd=fmri_finalA.mean_fd;
cd D:\Canada_2020\OASIS\reports\fmri

ICs={'IC1';'IC2';'IC3';'IC4';'IC5';'IC6';'IC7';'IC8';'IC9';'IC10';'IC11';'IC12';'IC13';'IC14';'IC15';'IC16';'IC17';'IC18';'IC19';'IC20';'IC21'};
combs=allcomb(ICs, ICs);D=21;
for i=1:441; comb(i)=cellstr(strcat(combs{i,2}, combs{i,1})); end
ICs=reshape(comb, [21,21]); clear comb combs i
ICs_vector=ICs(triu(ones(D),1)==1)';

ix=abs(clinical_final.delta_fmri)<365 & fmri_finalA.Motion_QC<=1 & clinical_final.dx2num_1neuro_2PD_3thyroid_4CVD_5medicineCogDys_6B12def_7Seizur==0;%ix=logical(ones(675,1)); 

%% setting the missing values to NaN instead of 0
cognitive_final.TRAILB(cognitive_final.Animal==0)=NaN;
cognitive_final.TRAILBRR(cognitive_final.Animal==0)=NaN;
cognitive_final.tmt_cor=cognitive_final.TRAILB+5*cognitive_final.TRAILBRR;
cognitive_final.BOSTON(cognitive_final.Animal==0)=NaN;
cognitive_final.MEMUNITS(cognitive_final.Animal==0)=NaN;
cognitive_final.aLOGIMEM(cognitive_final.Animal==0)=NaN;
%% pls setup
ids=cellfun(@isempty,  strfind(ICs_vector, 'IC4')) & cellfun(@isempty,  strfind(ICs_vector, 'IC8')) & cellfun(@isempty,  strfind(ICs_vector, 'IC19')) & cellfun(@isempty,  strfind(ICs_vector, 'IC15'));% & cellfun(@isempty,  strfind(ICs_vector, 'IC2IC'));
%ispresent = cellfun(@(s) ~isempty(strfind(string, s)), elements)

Y=[cognitive_final.tmt_cor,cognitive_final.WAIS ];%   clinical_final.Age
Y=[cognitive_final.MEMUNITS,clinical_final.Age];%  %(isnan(Y)')' ,cognitive_final.aLOGIMEM,, 
X=part_icad25A;%(:,ids==1); %X(96,:)=NaN;

%%% remove nans - 
naninx=sum(isnan(X)')'>0| sum(isnan(Y)')' >0  |clinical_final.cdr>1.5 |~ix; Y=Y(naninx==0,:);  X=X(naninx==0,:); clin=clinical_final(naninx==0,:); mean_fd=mean_fd(naninx==0,:); prop_signal=prop_signal(naninx==0,:); %|clinical_final.Age<50 |mean_fd>0.1
clinical_pls=clinical_final(naninx==0,:); 
Y=zscore(Y); X=zscore(X);x=X;
ncomp=4
%max(partialcorr(X,Y, [clin.Age,prop_signal,mean_fd, clin.delta_fmri]))

%%% regress out sex and age and potentially other covariates
%mdl = fitlm(table(clin.Age,  clin.M_F, clin.delta_fmri, Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
%mdl = fitlm(table(clin.Age,  clin.M_F, clin.delta_fmri, Y(:,2)) ); Y(:,2)=mdl.Residuals.Raw; %clin.Age, 
%clear x; for i=1:length(X(1,:)); mdl = fitlm(table(clin.M_F,mean_fd, X(:,i)) ); x(:,i)=mdl.Residuals.Raw;end

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR
%%% explore the PLS components and their correlation with Y
%figure;imagesc(corr(XS,Y)); 
 XS=XS*-1; c=corr(XS, Y) %1.TMT  2.Tower  3.PAL  4.DSST
combs=allcomb([1:3], [1:2]);figure(13);
for i=1:length(combs)
hold off; ix1=combs(i,1);ix2=combs(i,2); subplot(3,2,i); 
scatter(XS(clin.cdr==0,ix1), Y(clin.cdr==0,ix2), 5,'b', 'filled'); 
hold on;scatter(XS(clin.cdr>0,ix1), Y(clin.cdr>0,ix2), 5, 'r', 'filled')
mdl= fitlm(XS(:, ix1), Y(:,ix2)); [ypred,yci] = predict(mdl,XS(:, ix1), 'Alpha',0.001); hold on
yci=[yci,XS(:, ix1)]; yci=sortrows(yci, 3);
plot(XS(:, ix1), ypred, 'k', 'LineWidth', 2);
plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);
plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5); xlim([-.22 .22]); 
text(-0.2, -2, strcat('R=',num2str(round(c(ix1, ix2),3))) ); %annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'string', num2str(c(ix1, ix2)))
end; clear mdl ypred yci mdl ix1 ix2 combs c
anova1(XS(:,1),clin.cdr)
%% permutation testing
permutations=5000;   
allobservations=Y; 
for ncomp=1:4
    
    for n = 1:permutations
    % selecting either next combination, or random permutation
    permutation_index = randperm(length(allobservations));
    % creating random sample based on permutation index
    randomSample = allobservations(permutation_index,:);
    % running the PLS for this permutation
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,randomSample,ncomp);
    Rsq(n) = sum(PCTVAR(2,:));
    Rsq1(n) = sum(PCTVAR(1,:));
    end
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);
    p(ncomp)=sum(sum(PCTVAR(2,:))<Rsq')/permutations
    p_1(ncomp)=sum(sum(PCTVAR(1,:))<Rsq1')/permutations
end
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR

%% bootstrapping to get the func connectivity weights for PLS1, 2 and 3
dim=4
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,dim);PCTVAR
PLS3w=stats.W(:,3);

bootnum=5000;
PLS3weights=[];

parfor i=1:bootnum
    i;
    myresample = randsample(size(x,1),size(x,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=x(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    newW=stats.W(:,3);%extract PLS2 weights
    if corr(PLS3w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS3weights=[PLS3weights,newW]; %store (ordered) weights from this bootstrap run    
end

PLS3sw=std(PLS3weights');
plsweights3=PLS3w./PLS3sw';

%plsweights=zeros(1,210); plsweights(ids)=plsweights3;
%plsweights3=plsweights;
ICs_vector(plsweights3>3)
ICs_vector(plsweights3<-3)
c
%doublecheck:%ICs_vector_cut=ICs_vector(ids);ICs_vector_cut(plsweights3>3)ICs_vector_cut(plsweights3<-3)


%% visuals for individual brain connectivities
rvalues=corr(X, Y(:,1), 'rows', 'pairwise')
figure; scatter(rvalues, plsweights3)

ids=1:210;figure(2);hold off
for i=ids(plsweights3>3)
    j=find(i==ids(plsweights3>3))
    allobservations=full_icad25(:,i);
    subplot(3,3,j); scatter(X(clin.cdr==0,i), Y(clin.cdr==0,1), 3,'filled', 'k'); hold on; scatter(X(clin.cdr==0.5,i), Y(clin.cdr==0.5,1), 3,'filled', 'k');
    m=fitlm(X(:,i),Y(:,1));
    [ypred,yci] = predict(m,X(:,i), 'Alpha',0.001); hold on
    yci=[yci,X(:,i)]; yci=sortrows(yci, 3);
    plot(X(:,i), ypred, 'k', 'LineWidth', 2);
    plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);
    plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5); 
    r=round(corr(X(:,i), Y(:,1), 'rows', 'pairwise'),2)
    %xlim([7 25]); 
    %text(8, 5, strcat('R=',num2str(r)) ); %ylim([4 14]);
    %text(-4, min(Y(:,1)), strcat('R=',num2str(r)) ); %ylim([-1 1]);
    xlabel(ICs_vector{i}); clear r; 
end
amy=[1.09400000000000;0.899000000000000;0.970000000000000;0.853000000000000;1.08100000000000;1.00400000000000;1.10500000000000;0.943000000000000;1.34000000000000;3.47800000000000;1.32500000000000;1.01200000000000;3.15200000000000;1.06000000000000;1.02100000000000;1.08500000000000;2.87800000000000;1.14800000000000;1.10300000000000;1.22400000000000;1.07300000000000;1.22000000000000;0.750000000000000;0.998000000000000;0.946000000000000;1.04000000000000;1.33100000000000;3.98400000000000;0.808000000000000;2.03000000000000;3.24000000000000;1.16100000000000;0.977000000000000;0.918000000000000;1.40900000000000;1.01000000000000;2.56000000000000;1.00500000000000;0.846000000000000;1.06000000000000;1.02100000000000;1.18600000000000;2.46900000000000;0.881000000000000;0.942000000000000;2.50900000000000;0.986000000000000;4.18300000000000;1.08100000000000;1.07500000000000;0.989000000000000;1.58100000000000;2.03300000000000;1.06400000000000;1.76500000000000;0.882000000000000;1.27900000000000;1.01400000000000;1.86800000000000;1.01300000000000;2.23500000000000;1.21800000000000;0.545000000000000;0.749000000000000;1.01300000000000;2.41600000000000;0.805000000000000;1.14500000000000;1.85900000000000;2.14200000000000;3.24500000000000;3.54200000000000;2.86800000000000;2.24800000000000;3.39200000000000;0.987000000000000;1.00700000000000;0.963000000000000;1.06600000000000;1.10500000000000;0.948000000000000;1.11100000000000;1.78100000000000;1.07500000000000;1.43700000000000;1.02600000000000;2.47600000000000;1.07100000000000;1.99100000000000;0.827000000000000;3.71200000000000;1.16900000000000;1.02000000000000;1.14200000000000;3.33600000000000;1.09700000000000;1.65800000000000;1.29700000000000;0.936000000000000;0.897000000000000;0.958000000000000;1.06000000000000;1.03700000000000;3.72300000000000;2.19800000000000;1.24800000000000;2.96200000000000;1.19400000000000;2.03200000000000;1.07400000000000;2.95200000000000;1.24200000000000;1.80100000000000;0.744000000000000;2.15000000000000;0.960000000000000;0.594000000000000;0.737000000000000;1.05500000000000;2.29700000000000;0.987000000000000;2.23300000000000;1.91600000000000;1.84500000000000;1.49400000000000;2.58000000000000;1.03700000000000;1.56200000000000;0.975000000000000;1.06000000000000;1.12500000000000;3.81900000000000;0.824000000000000;1.34900000000000;2.87700000000000;0.760000000000000;1.21600000000000;2.36600000000000;0.985000000000000;0.738000000000000;1.67900000000000;0.987000000000000;2.34300000000000;1.57300000000000;2.24400000000000;2.53100000000000;2.04200000000000;1.06000000000000;2.19100000000000;1.39000000000000;1.04900000000000;2.51900000000000;4.21400000000000;2.38500000000000;2.54500000000000;1.06700000000000;0.977000000000000;1.02700000000000;0.998000000000000;0.889000000000000;2.59900000000000;5.41800000000000;1.06300000000000;2.25800000000000;1.22100000000000;0.966000000000000;2.40000000000000;2.04200000000000;1.02100000000000;1.81900000000000;2.27000000000000;0.827000000000000;0.848000000000000;1.84700000000000;1.70300000000000;1.02400000000000;0.985000000000000;0.659000000000000;1.72300000000000;2.77800000000000;1.01800000000000;3.27700000000000;1.20700000000000;2.19600000000000;1.96100000000000;2.89700000000000;2.98100000000000;1.01100000000000;2.61600000000000;1.03600000000000;3.25100000000000;1.22300000000000;1.03400000000000;1.05500000000000;1.03600000000000;2.45800000000000;2.49200000000000;1.62600000000000;2.73000000000000;0.999000000000000;1.04600000000000;1.31600000000000;1.11000000000000;0.961000000000000;1.39500000000000;3.80400000000000;1.10600000000000;1.03400000000000;3.61500000000000;3.47900000000000;1.65500000000000;0.989000000000000;0.991000000000000;1.09100000000000;0.994000000000000;1.05900000000000;0.939000000000000;1.90900000000000;1.99000000000000;1.20800000000000;2.81300000000000;0.999000000000000;1.38300000000000;3.10300000000000;1.02700000000000;0.966000000000000;2.68800000000000;1.00900000000000;0.980000000000000;2.48200000000000;0.855000000000000;1.91400000000000;0.995000000000000;1.06300000000000;0.722000000000000;3.04400000000000;1.96000000000000;1.70600000000000;1.04000000000000;1.03100000000000;0.893000000000000;1.01100000000000;2.39700000000000;0.856000000000000;2.32200000000000;0.923000000000000;0.994000000000000;2.34800000000000;1.12100000000000;0.949000000000000;0.768000000000000;0.947000000000000;1.13600000000000;0.857000000000000;2.48300000000000;0.882000000000000;1.17000000000000;2.87600000000000;3.96700000000000;1.08900000000000;2.04100000000000;1.28400000000000;0.924000000000000;1.10600000000000;0.755000000000000;2.32700000000000;0.963000000000000;2.01500000000000;2.38900000000000;1.08100000000000;0.920000000000000;0.984000000000000;2.10800000000000;0.949000000000000;0.996000000000000;0.795000000000000;0.891000000000000;3.41200000000000;2.02100000000000;1.26800000000000;1.27100000000000;1.09800000000000;1.01000000000000;0.983000000000000;1.95000000000000;1.07500000000000;3.07600000000000;0.736000000000000;0.933000000000000;0.946000000000000;1.00700000000000;1.43500000000000;1.12500000000000;2.27700000000000;1.01400000000000;1.55700000000000;1.05000000000000;1.15300000000000;2.14000000000000;2.30600000000000;3.05200000000000;1.06200000000000;2.55600000000000;0.899000000000000;1.34200000000000;0.960000000000000;1.01200000000000;2.37000000000000;1.55200000000000;1.63700000000000;2.73500000000000;1.28900000000000;1.55500000000000;2.48800000000000;2.16000000000000;1.11800000000000;0.973000000000000;0.987000000000000;0.901000000000000;1.53000000000000;1.93500000000000;0.886000000000000;0.805000000000000;0.818000000000000;2.13600000000000;1.07600000000000;2.00700000000000;0.905000000000000;2.75900000000000;1.24200000000000;0.949000000000000;0.994000000000000;1.06000000000000;1.15300000000000;1.38300000000000;1.03300000000000;2.29700000000000;1.00900000000000;2.92500000000000;1.58500000000000;2.63700000000000;3.00800000000000;1.89300000000000;1.68500000000000;1.09200000000000;3.41800000000000;1.03200000000000;3.41700000000000;2.09200000000000;1.04600000000000;1.03800000000000;0.952000000000000;0.862000000000000;2.41200000000000;2.00300000000000;0.921000000000000;0.801000000000000;1.10700000000000;2.53200000000000;0.992000000000000;0.949000000000000;1.08400000000000;3.75900000000000;0.926000000000000;1.02800000000000;1.44000000000000;1.04400000000000;0.933000000000000;0.870000000000000;1.21400000000000;1.40900000000000;0.853000000000000;1.22100000000000;4.01500000000000;1.05400000000000;1.39600000000000;1.14900000000000;1.02700000000000;0.900000000000000;0.808000000000000;2.33600000000000;3.49800000000000;1.92700000000000;0.651000000000000;0.927000000000000;1.41700000000000;1.12800000000000;1.04300000000000;1.00800000000000;1.03200000000000;2.74900000000000;0.923000000000000];
figure(3);hold off
for i=ids(plsweights3<-4)
    j=find(i==ids(plsweights3<-4))
    allobservations=part_icad25(:,i);
    subplot(3,3,j); scatter(amy, X(:,i), 3,'filled', 'k'); hold on; 
    m=fitlm(amy, X(:,i));
    [ypred,yci] = predict(m,amy, 'Alpha',0.001); hold on
    yci=[yci,amy]; yci=sortrows(yci, 3);
    plot(amy, ypred, 'k', 'LineWidth', 2);
    plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);
    plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5); 
    r=round(corr(X(:,i), amy, 'rows', 'pairwise'),2)%xlim([7 25]); 
    %text(0, min(Y(:,1)), strcat('R=',num2str(r)) ); %ylim([-1 1]);
end

%% visualize the weights using circularGraph(x);
component=1; 
ica2yeo7=readtable('D:\Canada_2020\UK_biobank\reports\ica2yeo7.csv');
nodes=dlmread('D:\Canada_2020\UK_biobank\reports\icad25_nodesordered.csv'); tmp=nanmean(X(clin.cdr==0,:))';
%nodes=[plsweights3,(1:210)',tmp, nodes]; nodes=sortrows(nodes, component);

myLabel=ica2yeo7.Yeo7N;figure
upperlim=3; lowerlim=-3; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=plsweights3; Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=plsweights3; Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

%%% average connectivity visual
figure; 
Weights=tmp;
Weights(Weights>-0.04)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

Weights=tmp;
Weights(Weights<0.03)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

