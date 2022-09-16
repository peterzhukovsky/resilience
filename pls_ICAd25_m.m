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

mean_rms=dlmread('mean_rms_boldfiles_rep.txt');
subs_rms=readtable('mean_rms_boldfiles_ids_rep.txt', 'ReadVariableNames', 0);subs_rms.Properties.VariableNames={'Subject'};
subs_rms.mean_rms=mean_rms;
clinical_final=innerjoin(subs_rms, clinical_final);

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
Y=[cognitive_final.MEMUNITS];%  %(isnan(Y)')' ,cognitive_final.aLOGIMEM,, ,clinical_final.Age
X=part_icad25A;%(:,ids==1); %X(96,:)=NaN;

%%% remove nans - 
%naninx=sum(isnan(X)')'>0| sum(isnan(Y)')' >0  |clinical_final.cdr>1 |~ix; Y=Y(naninx==0,:);  X=X(naninx==0,:); clin=clinical_final(naninx==0,:); mean_fd=mean_fd(naninx==0,:); prop_signal=prop_signal(naninx==0,:); %|clinical_final.Age<50 |mean_fd>0.1
naninx=sum(isnan(X)')'>0| isnan(Y) >0  |clinical_final.cdr>1 |~ix | clinical_final.mean_rms>0.3; Y=Y(naninx==0,:);  X=X(naninx==0,:); clin=clinical_final(naninx==0,:); amy=amyloid_final(naninx==0,:);mean_fd=mean_fd(naninx==0,:); prop_signal=prop_signal(naninx==0,:); %|clinical_final.Age<50 |mean_fd>0.1
clinical_pls=clinical_final(naninx==0,:); 
Y=zscore(Y); X=zscore(X);x=X;
ncomp=3
%max(partialcorr(X,Y, [clin.Age,prop_signal,mean_fd, clin.delta_fmri]))

%%% regress out sex and age and potentially other covariates
%mdl = fitlm(table(clin.Age,  clin.M_F, clin.delta_fmri, Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
mdl = fitlm(table(clin.Age, clin.M_F, clin.delta_fmri, Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
clear x; for i=1:length(X(1,:)); mdl = fitlm(table(clin.M_F,clin.Age, clin.mean_rms, X(:,i)) ); x(:,i)=mdl.Residuals.Raw;end

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR
%%% explore the PLS components and their correlation with Y
%figure;imagesc(corr(XS,Y)); 
%XS=XS*-1; 
c=corr(XS, Y) %1.TMT  2.Tower  3.PAL  4.DSST
combs=allcomb([1:3], [1:1]);figure(13);
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
figure; histogram(Rsq'); sum(PCTVAR')

%% bootstrapping to get the func connectivity weights for PLS1, 2 and 3
dim=4
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

%plsweights=zeros(1,210); plsweights(ids)=plsweights3;
%plsweights3=plsweights;
ICs_vector(plsweights1>3)
ICs_vector(plsweights1<-3)

ICs_vector(plsweights2>3)
ICs_vector(plsweights2<-3)

%doublecheck:%ICs_vector_cut=ICs_vector(ids);ICs_vector_cut(plsweights3>3)ICs_vector_cut(plsweights3<-3)

[coef,scores,~,~,explained] = pca(x);
corr(coef(:,1), plsweights1)
%% visuals for individual brain connectivities
rvalues=corr(X, Y(:,1), 'rows', 'pairwise')
figure; scatter(rvalues, plsweights3)

ids=1:210;figure(2);hold off
for i=ids(plsweights1>3)
    j=find(i==ids(plsweights1>3))
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

ids=1:210;figure;hold off
for i=ids(plsweights1>3)
    j=find(i==ids(plsweights1>3))
    allobservations=part_icad25(:,i);
    subplot(3,3,j); scatter(amy.PET_fSUVR_rsf_TOT_CORTMEAN, X(:,i), 3,'filled', 'k'); hold on; 
    m=fitlm(amy.PET_fSUVR_rsf_TOT_CORTMEAN, X(:,i));
    [ypred,yci] = predict(m,amy.PET_fSUVR_rsf_TOT_CORTMEAN, 'Alpha',0.001); hold on
    yci=[yci,amy.PET_fSUVR_rsf_TOT_CORTMEAN]; yci=sortrows(yci, 3);
    plot(amy.PET_fSUVR_rsf_TOT_CORTMEAN, ypred, 'k', 'LineWidth', 2);
    plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);
    plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5); 
    r=round(corr(X(:,i), amy.PET_fSUVR_rsf_TOT_CORTMEAN, 'rows', 'pairwise'),2)%xlim([7 25]); 
    text(0, min(Y(:,1)), strcat('R=',num2str(r)) ); %ylim([-1 1]);
     xlabel(ICs_vector{i}); clear r; 
end

%% visualize the weights using circularGraph(x);
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
