clear
fmri_metrics=readtable('D:\Canada_2020\OASIS\data\fMRI\fmri_preproc_data.xlsx');
load('D:\Canada_2020\OASIS\reports\fmri\ICAd25_rmats_July2021.mat');
load('D:\Canada_2020\OASIS\reports\fmri\OASIS_rsfc_evaluation_July2021.mat');
subjects=subjects{1};

cd D:\Canada_2020\OASIS\data\fMRI\clinical_and_demographic
amyloidA=readtable('amyA.csv'); 
clinicalA=readtable('clinA.csv'); clinicalA.Age=clinicalA.ageAtEntry+fmri_metrics.RunA_days/365;
cognitiveA=readtable('cogA.csv'); clinicalA.delta_AmyloidA=fmri_metrics.RunA_days-amyloidA.time_since_start;
clinicalA.HYPERTEN(clinicalA.test==0)=NaN;clinicalA.HYPERCHO(clinicalA.test==0)=NaN;clinicalA.bmi=(clinicalA.weight*0.453592)./((clinicalA.height*2.54/100).^2);
cognitiveA.MEMUNITS(cognitiveA.Animal==0)=NaN;

cd D:\Canada_2020\OASIS\reports\fmri

for sub=1:length(subjects)
rmat=reshape(Fnetmats_all(sub, :,:), [21 21]);D=21;
full_icad25(sub,:)=rmat(triu(ones(D),1)==1)';
rmat=reshape(Pnetmats_all(sub, :,:), [21 21]);D=21;
part_icad25(sub,:)=rmat(triu(ones(D),1)==1)';
end; clear rmat %% reshape the netmats to one line from the square form
%%%%% didnt pass registration QC
full_icad25(78,:)=NaN;full_icad25(93,:)=NaN;full_icad25(99,:)=NaN;full_icad25(253,:)=NaN;full_icad25(262,:)=NaN; full_icad25(429,:)=NaN; full_icad25(534,:)=NaN; full_icad25(569,:)=NaN; 


clinicalA(isnan(Fnetmats_all(:,1,1)),:)=[];cognitiveA(isnan(Fnetmats_all(:,1,1)),:)=[];amyloidA(isnan(Fnetmats_all(:,1,1)),:)=[];
subjects(isnan(Fnetmats_all(:,1,1)),:)=[];prop_signal(isnan(Fnetmats_all(:,1,1)),:)=[]; mean_fd(isnan(Fnetmats_all(:,1,1)),:)=[]; 
part_icad25(isnan(Fnetmats_all(:,1,1)),:)=[]; full_icad25(isnan(Fnetmats_all(:,1,1)),:)=[];
fmri_metrics(isnan(Fnetmats_all(:,1,1)),:)=[]; %% cut down the sheets to only include the non-NaN values, ie those with amyloid data

%sum(~isnan(Fnetmats_all(:,1,1)))
clinicalA.cdrcat=ones(length(subjects),1);clinicalA.cdrcat(clinicalA.cdrcat==1)=NaN; clinicalA.cdrcat(clinicalA.cdr==0)=0;clinicalA.cdrcat(clinicalA.cdr==0.5)=0.5;clinicalA.cdrcat(clinicalA.cdr==1)=1;
clinicalA.risk=zeros(length(subjects),1); 

tmp=table(amyloidA.tracer, amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN);mdl=fitlm(tmp); 
%figure; scatter(mdl.Residuals.Raw,amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN)
clinicalA.risk=mdl.Residuals.Raw;

ICs={'IC1';'IC2';'IC3';'IC4';'IC5';'IC6';'IC7';'IC8';'IC9';'IC10';'IC11';'IC12';'IC13';'IC14';'IC15';'IC16';'IC17';'IC18';'IC19';'IC20';'IC21'};
combs=allcomb(ICs, ICs);
for i=1:441; comb(i)=cellstr(strcat(combs{i,2}, combs{i,1})); end
ICs=reshape(comb, [21,21]); clear comb combs i
ICs_vector=ICs(triu(ones(D),1)==1)';

mean_rms=dlmread('mean_rms_boldfiles.txt');
subs_rms=readtable('mean_rms_boldfiles_ids.txt', 'ReadVariableNames', 0);subs_rms.Properties.VariableNames={'Subject'};
subs_rms.mean_rms=mean_rms;
clinicalA=innerjoin(subs_rms, clinicalA);


%% Education analysis - step 2
ix=amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN~=0 & abs(clinicalA.delta_AmyloidA)<365 & clinicalA.dx2num==0 & mean_fd<0.1 & (clinicalA.apoe==34|clinicalA.apoe==44|clinicalA.apoe==33);%ix=logical(ones(675,1)); 
ix=amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN~=0 & abs(clinicalA.delta_AmyloidA)<365 & clinicalA.dx2num==0 & clinicalA.mean_rms <0.3;% & (clinicalA.apoe==34|clinicalA.apoe==44|clinicalA.apoe==33);%ix=logical(ones(675,1)); 

permutations=1000;

for i=1:length(part_icad25(1,:))
        allobservations=clean(part_icad25(:,i),3); clin=clinicalA(~isnan(allobservations),:);p=prop_signal(~isnan(allobservations));m=clinicalA.mean_rms(~isnan(allobservations)); a=allobservations(~isnan(allobservations)); 
        i
        parfor n = 1:permutations; 
        permutation_index = randperm(length(a));
        randomSample = a(permutation_index,:);
        T=table(clin.Education, clin.risk, clin.cdrcat, clin.Age, clin.M_F,p,m, abs(clin.delta_AmyloidA), randomSample ); %clear clin p m a
        T.Properties.VariableNames={'edu', 'risk', 'cdrcat', 'Age', 'M_F','prop_signal', 'mean_fd', 'abs_delta_AmyloidA', 'FC'};T=T(ix(~isnan(allobservations)),:);
        mdl = fitlm(T,'FC~edu+risk+cdrcat+Age+M_F+prop_signal+mean_fd+abs_delta_AmyloidA+Age*M_F'); 
        anovamdl=anova(mdl); F(1,n)=anovamdl.F(1);
        end; 
        F_permedu(i,:)=F;
  clear T;T=table(clinicalA.Education, clinicalA.risk, clinicalA.cdrcat, clinicalA.Age, clinicalA.M_F,prop_signal,clinicalA.mean_rms, abs(clinicalA.delta_AmyloidA), clean(allobservations,3) ); 
  T.Properties.VariableNames={'edu', 'risk', 'cdrcat', 'Age', 'M_F','prop_signal', 'mean_fd', 'abs_delta_AmyloidA', 'FC'};T=T(ix,:);
  mdl = fitlm(T,'FC~edu+risk+cdrcat+Age+M_F+prop_signal+mean_fd+abs_delta_AmyloidA+Age*M_F'); anovamdl=anova(mdl);
  p_perm_edu(i)=1-sum(anovamdl.F(1)>F_permedu(i,:))/permutations; 
  Tedu(1,i)=mdl.Coefficients.tStat(2); Pedu(1,i)=mdl.Coefficients.pValue(2);
end
ICs_vector(p_perm_edu<0.05)
square_P_part = zeros(D);
square_P_part(triu(ones(D),1)>0) = p_perm_edu;
square_T_part = zeros(D); square_T_part(triu(ones(D),1)>0) = Tedu; square_T_part(square_P_part>0.05)=0;
square_T = square_T_part'; 
figure;imagesc(square_T); set(gca,'XTick',[1:21], 'YTick',[1:21]); %clear square_* 

ids=1:210;figure(131);
for i=ids(p_perm_edu<0.05)
    j=find(i==ids(p_perm_edu<0.05))
    allobservations=part_icad25(:,i);
    clear T;T=table(clinicalA.Education, clinicalA.risk, clinicalA.cdrcat, clinicalA.Age, clinicalA.M_F,prop_signal,clinicalA.mean_rms, abs(clinicalA.delta_AmyloidA), clean(allobservations,3) ); 
    T.Properties.VariableNames={'edu', 'risk', 'cdrcat', 'Age', 'M_F','prop_signal', 'mean_rms', 'abs_delta_AmyloidA', 'FC'}; T=T(ix,:);
    mdl = fitlm(T,'FC~edu+risk+cdrcat+Age+M_F+prop_signal+mean_rms+abs_delta_AmyloidA+Age*M_F'); anovamdl=anova(mdl);
    mdl_p(j)=mdl.coefTest;
    subplot(3,5,j); scatter(T.edu(T.cdrcat==0), mdl.Fitted(T.cdrcat==0), 5,'filled', 'b'); hold on; scatter(T.edu(T.cdrcat>=0.5), mdl.Fitted(T.cdrcat>=0.5), 5,'filled', 'r');
    m= fitlm(T.edu,mdl.Fitted);
    [ypred,yci] = predict(m,T.edu, 'Alpha',0.001); hold on
    yci=[yci,T.edu]; yci=sortrows(yci, 3);
    plot(T.edu, ypred, 'k', 'LineWidth', 2);
    plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);
    plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5); 
    r=round(corr(T.edu, mdl.Fitted, 'rows', 'pairwise'),2);xlim([7 25]); 
    %text(8, 5, strcat('R=',num2str(r)) ); ylim([4 14]); 
    text(8, min(ypred), strcat('R=',num2str(r)) ); %ylim([-1 1]);
    xlabel(ICs_vector{i}); clear r; 
end


figure;plotEffects(mdl)

  %% Memory
  
cognitiveA.TRAILB(cognitiveA.Animal==0)=NaN;
cognitiveA.TRAILBRR(cognitiveA.Animal==0)=NaN;
cognitiveA.tmt_cor=cognitiveA.TRAILB+5*cognitiveA.TRAILBRR;
cognitiveA.BOSTON(cognitiveA.Animal==0)=NaN;
cognitiveA.MEMUNITS(cognitiveA.Animal==0)=NaN;
cognitiveA.aLOGIMEM(cognitiveA.Animal==0)=NaN;
clinicalA.delta_cog_t=abs(cognitiveA.Days_since_start-fmri_metrics.RunA_days);

%ix=abs(clinicalA.delta_cog_t)<365;%ix=logical(ones(969,1)); 
ix=amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN~=0 & abs(clinicalA.delta_cog_t)<365 & clinicalA.dx2num==0 %& (clinicalA.apoe==34|clinicalA.apoe==44|clinicalA.apoe==33);%ix=logical(ones(675,1)); 

  clear T;T=table(clinicalA.apoe, cognitiveA.tmt_cor, cognitiveA.BOSTON,cognitiveA.MEMUNITS, cognitiveA.aLOGIMEM, clinicalA.Education, clinicalA.risk, clinicalA.cdrcat, clinicalA.Age, clinicalA.M_F, amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN, abs(clinicalA.delta_cog_t),  abs(clinicalA.delta_AmyloidA) ); 
  T.Properties.VariableNames={'apoe','TMT', 'BOSTON', 'MEMUNITS', 'LOGIMEM','education', 'risk', 'cdrcat', 'Age', 'M_F','amyloid', 'abs_delta_cogA', 'abs_delta_AmyloidA'};T=T(ix,:);

 mdl = fitlm(T,'LOGIMEM~education+amyloid+cdrcat+Age+M_F+abs_delta_AmyloidA+abs_delta_cogA+abs_delta_cogA+Age*M_F')
  subplot(2,2,3);
  scatter(T.education(T.cdrcat==0), mdl.Fitted(T.cdrcat==0), 5,'filled', 'b'); hold on; scatter(T.education(T.cdrcat==1), mdl.Fitted(T.cdrcat==1), 5,'filled', 'r');
    m= fitlm(T.education,mdl.Fitted);
    [ypred,yci] = predict(m,T.education, 'Alpha',0.001); hold on
    yci=[yci,T.education]; yci=sortrows(yci, 3);    plot(T.education, ypred, 'k', 'LineWidth', 2);    plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);    plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5);     r=round(corr(T.education, mdl.Fitted, 'rows', 'pairwise'),2);xlim([7 25]);     text(8, min(mdl.Fitted), strcat('R=',num2str(r)) ); %ylim([4 14]);
      
  mdl = fitlm(T,'MEMUNITS~education+amyloid+apoe+cdrcat+Age+M_F+abs_delta_AmyloidA+abs_delta_cogA+Age*M_F')
  subplot(2,2,4);
  scatter(T.education(T.cdrcat==0), mdl.Fitted(T.cdrcat==0), 5,'filled', 'b'); hold on; scatter(T.education(T.cdrcat==1), mdl.Fitted(T.cdrcat==1), 5,'filled', 'r');
    m= fitlm(T.education,mdl.Fitted);
    [ypred,yci] = predict(m,T.education, 'Alpha',0.001); hold on
    yci=[yci,T.education]; yci=sortrows(yci, 3);    plot(T.education, ypred, 'k', 'LineWidth', 2);    plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);    plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5);     r=round(corr(T.education, mdl.Fitted, 'rows', 'pairwise'),2);xlim([7 25]);     text(8, min(mdl.Fitted), strcat('R=',num2str(r)) ); %ylim([4 14]);
  

  %% demographics
  

%clinicalA.risk( ~(clinicalA.apoe==34 | clinicalA.apoe==24 | clinicalA.apoe==44) & ~amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN<1.19)=0;
ix=amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN~=0 & clinicalA.cdr<1 & abs(clinicalA.delta_AmyloidA)<365 & clinicalA.dx2num==0 & mean_fd<0.1 & (clinicalA.apoe==34|clinicalA.apoe==44|clinicalA.apoe==33);%ix=logical(ones(675,1)); 

clinicalA.sex10=zeros(969,1); clinicalA.sex10(ismember(clinicalA.M_F, 'M'))=1;
clin=clinicalA(ix~=0,:);
clin.m=mean_fd(ix~=0,:);
clin.p=prop_signal(ix~=0,:);
amyA=amyloidA(ix~=0,:);

tmp=table(clin.sex10,clin.Education, clin.Age,clin.mmse,amyA.PET_fSUVR_rsf_TOT_CORTMEAN, clin.m, clin.cdrcat, clin.risk, 'VariableNames', {'sex','Education','Age','mmse', 'amyloid','meanfd', 'cdr','risk'})

statsm=grpstats(tmp, {'cdr'}, {'mean'})         %statsm=statsm([1 2 5 3 4], :)

sum(strcmp(amyA.tracer, 'PIB'))
[tbl,chi2,p]=crosstab(clin.sex10, clin.risk, clin.cdrcat)
[tbl,chi2,p]=crosstab(clin.sex10, clin.risk)
[tbl,chi2,p]=crosstab(clin.sex10, clin.cdrcat)

[tbl,chi2,p]=crosstab(clin.Race, clin.risk, clin.cdrcat)

anova1(tmp.Age, tmp.risk)
anova1(tmp.Age, tmp.cdr)

anova1(tmp.Education, tmp.risk)
anova1(tmp.Education, tmp.cdr)

anova1(tmp.mmse, tmp.risk)
anova1(tmp.mmse, tmp.cdr)

anova1(tmp.amyloid, tmp.risk)
anova1(tmp.amyloid, tmp.cdr)

anova1(tmp.meanfd, tmp.risk)
anova1(tmp.meanfd, tmp.cdr)


%% 
ICs={'IC1';'IC7';'IC9';'IC13';'IC14';'IC20';'IC21';'IC5';'IC6';'IC16';'IC10';'IC11';'IC12';'IC17';'IC3';'IC2';'IC4';'IC8';'IC19';'IC15';'IC18'};
combs=allcomb(ICs, ICs);
for i=1:441; comb(i)=cellstr(strcat(combs{i,2}, combs{i,1})); end
ICs=reshape(comb, [21,21]); clear comb combs i
ICs_vector=ICs(triu(ones(D),1)==1)';

ic_mapping=readtable('ic_order_mapping.csv');
tmp=ic_mapping.original_order_num;

square_P_full = zeros(D); square_P_full(triu(ones(D),1)>0) = p_perm_edu(tmp); 
square_T_full = zeros(D); square_T_full(triu(ones(D),1)>0) = Tedu(tmp); square_T_full(square_P_full>0.05)=0;
figure(1);imagesc(square_T_full); set(gca,'XTick',[1:21], 'YTick',[1:21]); %clear square_* 

square_P_full = zeros(D); square_P_full(triu(ones(D),1)>0) = p_perm_edu(tmp); 
square_T_full = zeros(D); square_T_full(triu(ones(D),1)>0) = Tedu(tmp); square_T_full(square_P_full>0.05)=0;
figure(2);imagesc(square_T_full); set(gca,'XTick',[1:21], 'YTick',[1:21]); %clear square_* 

%% alternative
ica2yeo7=readtable('D:\Canada_2020\UK_biobank\reports\ica2yeo7.csv');
nodes=dlmread('D:\Canada_2020\UK_biobank\reports\icad25_nodesordered.csv'); 

plsweights3=Tedu; plsweights3(p_perm_edu>0.05)=0;
myLabel=ica2yeo7.Yeo7N;figure
upperlim=1; lowerlim=-1; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=plsweights3; Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=plsweights3; Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);



%% %
for i=1:length(part_icad25(1,:))
        allobservations=clean(full_icad25(:,i),3); clin=clinicalA(~isnan(allobservations),:);p=prop_signal(~isnan(allobservations));m=mean_fd(~isnan(allobservations)); a=allobservations(~isnan(allobservations)); 
        i
        parfor n = 1:permutations; 
        permutation_index = randperm(length(a));
        randomSample = a(permutation_index,:);
        T=table(clin.Education, clin.risk, clin.cdrcat, clin.Age, clin.M_F,p,m, abs(clin.delta_AmyloidA), randomSample ); %clear clin p m a
        T.Properties.VariableNames={'edu', 'risk', 'cdrcat', 'Age', 'M_F','prop_signal', 'mean_fd', 'abs_delta_AmyloidA', 'FC'};T=T(ix(~isnan(allobservations)),:);
        mdl = fitlm(T,'FC~edu+risk+cdrcat+Age+M_F+prop_signal+mean_fd+abs_delta_AmyloidA+Age*M_F'); 
        anovamdl=anova(mdl); F(1,n)=anovamdl.F(1);
        end; 
        F_permedu(i,:)=F;
  clear T;T=table(clinicalA.Education, clinicalA.risk, clinicalA.cdrcat, clinicalA.Age, clinicalA.M_F,prop_signal,mean_fd, abs(clinicalA.delta_AmyloidA), clean(allobservations,3) ); 
  T.Properties.VariableNames={'edu', 'risk', 'cdrcat', 'Age', 'M_F','prop_signal', 'mean_fd', 'abs_delta_AmyloidA', 'FC'};T=T(ix,:);
  mdl = fitlm(T,'FC~edu+risk+cdrcat+Age+M_F+prop_signal+mean_fd+abs_delta_AmyloidA+Age*M_F'); anovamdl=anova(mdl);
  p_perm_edu(i)=1-sum(anovamdl.F(1)>F_permedu(i,:))/permutations; 
  Tedu(1,i)=mdl.Coefficients.tStat(2); Pedu(1,i)=mdl.Coefficients.pValue(2);
end
square_P_full = zeros(D); square_P_full(triu(ones(D),1)>0) = p_perm_edu; 
square_T_full = zeros(D); square_T_full(triu(ones(D),1)>0) = Tedu; square_T_full(square_P_full>0.05)=0;
ICs_vector(p_perm_edu<0.05)
