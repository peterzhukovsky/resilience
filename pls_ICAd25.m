clear
fmri_metrics=readtable('D:\Canada_2020\OASIS\data\fMRI\fmri_preproc_data.xlsx');
load('D:\Canada_2020\OASIS\reports\fmri\ICAd25_rmats_July2021.mat');
load('D:\Canada_2020\OASIS\reports\fmri\OASIS_rsfc_evaluation_July2021.mat');
subjects=subjects{1};

cd D:\Canada_2020\OASIS\data\fMRI\clinical_and_demographic
amyloidA=readtable('amyA.csv'); 
clinicalA=readtable('clinA.csv'); clinicalA.Age=clinicalA.ageAtEntry+fmri_metrics.RunA_days/365;
cognitiveA=readtable('cogA.csv'); clinicalA.delta_AmyloidA=fmri_metrics.RunA_days-amyloidA.time_since_start;

cd D:\Canada_2020\OASIS\reports\fmri

for sub=1:length(subjects)
rmat=reshape(Fnetmats_all(sub, :,:), [21 21]);D=21;
full_icad25(sub,:)=rmat(triu(ones(D),1)==1)';
rmat=reshape(Pnetmats_all(sub, :,:), [21 21]);D=21;
part_icad25(sub,:)=rmat(triu(ones(D),1)==1)';
end; clear rmat %% reshape the netmats to one line from the square form

clinicalA(isnan(Fnetmats_all(:,1,1)),:)=[];cognitiveA(isnan(Fnetmats_all(:,1,1)),:)=[];amyloidA(isnan(Fnetmats_all(:,1,1)),:)=[];
subjects(isnan(Fnetmats_all(:,1,1)),:)=[];prop_signal(isnan(Fnetmats_all(:,1,1)),:)=[]; mean_fd(isnan(Fnetmats_all(:,1,1)),:)=[]; 
part_icad25(isnan(Fnetmats_all(:,1,1)),:)=[]; full_icad25(isnan(Fnetmats_all(:,1,1)),:)=[];
fmri_metrics(isnan(Fnetmats_all(:,1,1)),:)=[]; %% cut down the sheets to only include the non-NaN values, ie those with amyloid data
%%%%% didnt pass registration QC
full_icad25(78,:)=NaN;full_icad25(93,:)=NaN;full_icad25(99,:)=NaN;full_icad25(253,:)=NaN;full_icad25(262,:)=NaN; full_icad25(429,:)=NaN; full_icad25(534,:)=NaN; full_icad25(569,:)=NaN; 
part_icad25(78,:)=NaN;part_icad25(93,:)=NaN;part_icad25(99,:)=NaN;part_icad25(253,:)=NaN;part_icad25(262,:)=NaN; part_icad25(429,:)=NaN; part_icad25(534,:)=NaN; part_icad25(569,:)=NaN; 

clinicalA.delta_cog_t=abs(cognitiveA.Days_since_start-fmri_metrics.RunA_days);
ix=amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN~=0 & abs(clinicalA.delta_cog_t)<365 & clinicalA.dx2num==0;% & (clinicalA.apoe==34|clinicalA.apoe==44|clinicalA.apoe==33);%ix=logical(ones(675,1)); 
%for i=1:210;    full_icad25(:,i)=clean(full_icad25(:,i),3); end & mean_fd<0.1 

ICs={'IC1';'IC2';'IC3';'IC4';'IC5';'IC6';'IC7';'IC8';'IC9';'IC10';'IC11';'IC12';'IC13';'IC14';'IC15';'IC16';'IC17';'IC18';'IC19';'IC20';'IC21'};
combs=allcomb(ICs, ICs);
for i=1:441; comb(i)=cellstr(strcat(combs{i,2}, combs{i,1})); end
ICs=reshape(comb, [21,21]); clear comb combs i
ICs_vector=ICs(triu(ones(D),1)==1)';

mean_rms=dlmread('mean_rms_boldfiles.txt');
subs_rms=readtable('mean_rms_boldfiles_ids.txt', 'ReadVariableNames', 0);subs_rms.Properties.VariableNames={'Subject'};
subs_rms.mean_rms=mean_rms;
clinicalA=innerjoin(subs_rms, clinicalA);

%% setting the missing values to NaN instead of 0
cognitiveA.TRAILB(cognitiveA.Animal==0)=NaN;
cognitiveA.TRAILBRR(cognitiveA.Animal==0)=NaN;
cognitiveA.tmt_cor=cognitiveA.TRAILB+5*cognitiveA.TRAILBRR;
cognitiveA.BOSTON(cognitiveA.Animal==0)=NaN;
cognitiveA.MEMUNITS(cognitiveA.Animal==0)=NaN;
cognitiveA.aLOGIMEM(cognitiveA.Animal==0)=NaN;
%% pls setup
Y=[cognitiveA.tmt_cor,cognitiveA.WAIS ];%   clinicalA.Age
Y=[cognitiveA.MEMUNITS ];%  %(isnan(Y)')' , clinicalA.Age
X=part_icad25;

%%% remove nans
%naninx=sum(isnan(X)')'>0| sum(isnan(Y)')'>0 | clinicalA.mean_rms>0.4  |clinicalA.cdr>0.5 |~ix; Y=Y(naninx==0,:);  X=X(naninx==0,:); clin=clinicalA(naninx==0,:); mean_fd=mean_fd(naninx==0,:); prop_signal=prop_signal(naninx==0,:); 
naninx=sum(isnan(X)')'>0| isnan(Y) | clinicalA.mean_rms>0.3  |clinicalA.cdr>1 |~ix; Y=Y(naninx==0,:);  X=X(naninx==0,:); clin=clinicalA(naninx==0,:);amy=amyloidA(naninx==0,:); mean_fd=mean_fd(naninx==0,:); prop_signal=prop_signal(naninx==0,:); 

Y=zscore(Y); X=zscore(X);x=X;
clinical_pls_bl=clinicalA(naninx==0,:); 
ncomp=3
%max(partialcorr(X,Y, [clin.Age,prop_signal,mean_fd, clin.delta_cog_t]))

%%% regress out sex and age and potentially other covariates
%mdl = fitlm(table( clin.M_F, clin.delta_cog_t, Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
mdl = fitlm(table(clin.Age, clin.M_F, clin.delta_cog_t, Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
clear x; for i=1:length(X(1,:)); mdl = fitlm(table(clin.M_F,clin.Age, clin.mean_rms, X(:,i)) ); x(:,i)=mdl.Residuals.Raw;end

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR
%%% explore the PLS components and their correlation with Y
figure;imagesc(corr(XS,Y)); c=corr(XS, Y) %1.TMT  2.Tower  3.PAL  4.DSST
combs=allcomb([1:3], [1]);figure(13); 
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
for ncomp=1:3
    
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
figure; histogram(Rsq'); sum(PCTVAR')
%% bootstrapping to get the func connectivity weights for PLS1, 2 and 3
dim=3
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,dim);PCTVAR
PLS1w=stats.W(:,1);PLS2w=stats.W(:,2);PLS3w=stats.W(:,3);

bootnum=5000;
PLS1weights=[];PLS2weights=[];PLS3weights=[];

parfor i=1:bootnum
    i;
    myresample = randsample(size(x,1),size(x,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=x(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
    
    newW=stats.W(:,1);%extract PLS1 weights
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW]; %store (ordered) weights from this bootstrap run    
    
    newW=stats.W(:,2);%extract PLS2 weights
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
  
    newW=stats.W(:,3);%extract PLS3 weights
    if corr(PLS3w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS3weights=[PLS3weights,newW]; %store (ordered) weights from this bootstrap run    
end

PLS1sw=std(PLS1weights');
plsweights1=PLS1w./PLS1sw';
PLS2sw=std(PLS2weights');
plsweights2=PLS2w./PLS2sw';
PLS3sw=std(PLS3weights');
plsweights3=PLS3w./PLS3sw';

ICs_vector(plsweights1>3)
ICs_vector(plsweights1<-3)

ICs_vector(plsweights2>3)
ICs_vector(plsweights2<-3)

ICs_vector(plsweights3>3)
ICs_vector(plsweights3<-3)


%% visuals for individual brain connectivities
rvalues=corr(X, Y(:,1), 'rows', 'pairwise')
figure; scatter(rvalues, plsweights3)

ids=1:210;figure;hold off
for i=ids(plsweights1>3)
    j=find(i==ids(plsweights1>3))
    allobservations=part_icad25(:,i);
    subplot(3,3,j); scatter(x(clin.cdr==0,i), Y(clin.cdr==0,1), 3,'filled', 'k'); hold on; scatter(x(clin.cdr==0.5,i), Y(clin.cdr==0.5,1), 3,'filled', 'k');
    m=fitlm(x(:,i),Y(:,1));
    [ypred,yci] = predict(m,x(:,i), 'Alpha',0.001); hold on
    yci=[yci,X(:,i)]; yci=sortrows(yci, 3);
    plot(x(:,i), ypred, 'k', 'LineWidth', 2);
    plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);
    plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5); 
    r=round(corr(x(:,i), Y(:,1), 'rows', 'pairwise'),2);%xlim([7 25]); 
    %text(8, 5, strcat('R=',num2str(r)) ); %ylim([4 14]);
    text(-4, min(Y(:,1)), strcat('R=',num2str(r)) ); %ylim([-1 1]);
    xlabel(ICs_vector{i}); clear r; 
end

ids=1:210;figure;hold off
for i=ids(plsweights1<-3)
    j=find(i==ids(plsweights1<-3))
    allobservations=part_icad25(:,i);
    subplot(3,4,j); scatter(amy.PET_fSUVR_rsf_TOT_CORTMEAN, x(:,i), 3,'filled', 'k'); hold on; 
    m=fitlm(amy.PET_fSUVR_rsf_TOT_CORTMEAN, x(:,i));
    [ypred,yci] = predict(m,amy.PET_fSUVR_rsf_TOT_CORTMEAN, 'Alpha',0.001); hold on
    yci=[yci,amy.PET_fSUVR_rsf_TOT_CORTMEAN]; yci=sortrows(yci, 3);
    plot(amy.PET_fSUVR_rsf_TOT_CORTMEAN, ypred, 'k', 'LineWidth', 2);
    plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);
    plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5); 
    r=round(corr(x(:,i), amy.PET_fSUVR_rsf_TOT_CORTMEAN, 'rows', 'pairwise'),2)%xlim([7 25]); 
    text(0, min(Y(:,1)), strcat('R=',num2str(r)) ); %ylim([-1 1]);
    xlabel(ICs_vector{i}); clear r; 
end

%% visualize the weights using circularGraph(x);
component=1; 
ica2yeo7=readtable('D:\Canada_2020\UK_biobank\reports\ica2yeo7.csv');
nodes=dlmread('D:\Canada_2020\UK_biobank\reports\icad25_nodesordered.csv'); tmp=nanmean(X(clin.cdr==0,:))';
nodes=[plsweights3,(1:210)',tmp, nodes]; nodes=sortrows(nodes, component);

ic_mapping=readtable('D:\Canada_2020\OASIS\reports\fmri\ic_order_mapping.csv');
tmp=ic_mapping.original_order_num;

myLabel=ica2yeo7.Yeo7N;myLabel=myLabel([1 7 9 13 14 20 21 5 6 16 10 11 12 17 3 2 4 8 19 15 18]);
myLabel=repmat({''}, [21,1]); figure
upperlim=3; lowerlim=-3; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=plsweights1(tmp); Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=plsweights1(tmp); Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

myLabel=ica2yeo7.Yeo7N;figure
upperlim=3; lowerlim=-3; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=plsweights1; Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=plsweights1; Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

%% compare BL vs matched
load('D:\Canada_2020\OASIS\reports\fmri\HBM_revised\baseline_pls.mat','clinical_pls_bl');
load('D:\Canada_2020\OASIS\reports\fmri\HBM_revised\FU_pls.mat','clinical_pls');
for i=1:length(clinical_pls_bl.Subject)
idx=strcmp(clinical_pls_bl.Subject{i}, clinical_pls.Subject);
if sum(idx)>0
tmp(i)=sum(clinical_pls_bl.time_d(i)==clinical_pls.time_d(idx));
else
    tmp(i)=0;
end
end
sum(~tmp)/length(tmp)


for i=1:length(clinical_pls.Subject)
idx=strcmp(clinical_pls.Subject{i}, clinical_pls_bl.Subject);
tmp(i)=sum(idx);
end
sum(tmp)

clear tmp
for i=1:length(clinical_pls.Subject)
idx=strcmp(clinical_pls.Subject{i}, clinical_pls_bl.Subject);
if sum(idx)>0
tmp(i)=sum(clinical_pls.time_d(i)==clinical_pls_bl.time_d(idx));
else
    tmp(i)=0;
end
end
sum(tmp)/length(tmp)

sum(~tmp)


mean(clinical_pls_bl.Age)
mean(clinical_pls.Age)


%% extract data for mediation bootstrap 
load('D:\Canada_2020\OASIS\reports\fmri\HBM_revised\baseline_pls.mat');
load('D:\Canada_2020\OASIS\reports\fmri\HBM_revised\FU_pls.mat');
ICs_vector(plsweights1<-3)
Y
X(:,plsweights1<-3)
x(:,plsweights1<-3)
XS(:,1)
amy.PET_fSUVR_rsf_TOT_CORTMEAN

%% predict FU data based on baseline
clear
load('D:\Canada_2020\OASIS\reports\fmri\HBM_revised\baseline_pls.mat');
BETA_train=BETA;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,dim);PCTVAR;
r=corr(XS, Y)
load('D:\Canada_2020\OASIS\reports\fmri\HBM_revised\FU_pls.mat');
x=x(tmp==0,:);Y=Y(tmp==0);
y_pred = [ones(size(x ,1),1) x]*BETA_train;
[r, p]=corr(Y, y_pred)
figure;scatter(Y, y_pred, 'filled', 'k'); lsline;

figure(11); hold off; for i=1:4; subplot(4,1,i); scatter(y_pred(:,i), Y_test(:,i), 'x', 'k'); 
    p = polyfit(y_pred(:,i),Y_test(:,i),1); pred = polyval(p,y_pred(:,i)); hold on; plot(y_pred(:,i),pred,'r','LineWidth',3); set(gca,'xtick',[]); set(gca,'ytick',[]); end
