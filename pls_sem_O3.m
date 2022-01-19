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
ix=amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN~=0 & abs(clinicalA.delta_cog_t)<365 & clinicalA.dx2num==0 & mean_fd<0.1 & (clinicalA.apoe==34|clinicalA.apoe==44|clinicalA.apoe==33);%ix=logical(ones(675,1)); 
%for i=1:210;    full_icad25(:,i)=clean(full_icad25(:,i),3); end

ICs={'IC1';'IC2';'IC3';'IC4';'IC5';'IC6';'IC7';'IC8';'IC9';'IC10';'IC11';'IC12';'IC13';'IC14';'IC15';'IC16';'IC17';'IC18';'IC19';'IC20';'IC21'};
combs=allcomb(ICs, ICs);
for i=1:441; comb(i)=cellstr(strcat(combs{i,2}, combs{i,1})); end
ICs=reshape(comb, [21,21]); clear comb combs i
ICs_vector=ICs(triu(ones(D),1)==1)';

%% setting the missing values to NaN instead of 0
cognitiveA.TRAILB(cognitiveA.Animal==0)=NaN;
cognitiveA.TRAILBRR(cognitiveA.Animal==0)=NaN;
cognitiveA.tmt_cor=cognitiveA.TRAILB+5*cognitiveA.TRAILBRR;
cognitiveA.BOSTON(cognitiveA.Animal==0)=NaN;
cognitiveA.MEMUNITS(cognitiveA.Animal==0)=NaN;

%% pls setup
%Y=[cognitiveA.tmt_cor,cognitiveA.WAIS ];%   clinicalA.Age
%ids=[0	1	0	0	0	0	1	0	1	0	1	0	1	0	1	1	0	1	0	1	1	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	1	0	1	0	1	1	1	0	1	1	0	1	0	1	1	1	0	1	1	1	0	1	0	1	1	1	0	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	1	1	1	1	1	1	0	1	0	1	0	1	1	1	0	1	1	1	1	1	1	0	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
%ids=[0	1	0	0	0	0	1	0	1	0	1	0	1	0	1	1	0	1	0	1	1	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	1	0	0	1	0	1	0	1	1	1	0	1	0	0	1	1	0	1	0	1	1	1	0	1	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	1	0	0	1	1	1	0	1	0	1	0	1	1	1	0	1	0	0	1	1	1	0	1	1	0	1	0	1	1	1	0	1	0	0	1	1	1	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0]
%ids=[0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
ids=[0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
load('part_pls.mat', 'plsweights3')%load('permedu_part.mat', 'p_perm_edu'); %load('perm_part.mat', 'p_perm_inter')
ids=abs(plsweights3)>3;%ids=p_perm_edu<0.05
Y=[cognitiveA.MEMUNITS];%  %(isnan(Y)')' 
%tmp=table(amyloidA.tracer, amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN);mdl=fitlm(tmp); mdl.Residuals.Raw
%amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN
X=[amyloidA.PET_fSUVR_rsf_TOT_CORTMEAN, part_icad25(:,ids==1)];clear mdl

%%% remove nans
naninx=sum(isnan(X)')'>0| isnan(Y)>0  |clinicalA.cdr>0.5 |~ix; Y=Y(naninx==0,:);  X=X(naninx==0,:); clin=clinicalA(naninx==0,:);amy=amyloidA(naninx==0,:); mean_fd=mean_fd(naninx==0,:); prop_signal=prop_signal(naninx==0,:); 
%Y=zscore(Y); X=zscore(X);x=X;
%ncomp=3
% mediation analysis 
%Y=[cognitiveA.MEMUNITS , clinicalA.Age];%  %(isnan(Y)')' 
%X=part_icad25;
%mdl=fitlm(table(amy.tracer, X(:,1)) );X(:,1)=mdl.Residuals.Raw; % doesnt make a difference

%% MV is the dataset 
MV=[X,Y]; % this is the matrix with all the data
LV_labels=[{'Amyloid'},{'FC'}, {'MEMUNITS'}]';
MV_labels=[{'Amy'}, ICs_vector(ids==1), {'Memory'}]';

%% DB - this is the "inner" design matrix which specifies how the latent variables are linked to each other
% The DB matrix can be obtained drawing the structural model through the
% function SMDesign_gui
load('StructuralModelDesign.mat')
DB=[0,1,1;0,0,1;0,0,0];
%% DP - this is the "outer" design matrix which specifies what variables each of the latent variables is summarizing
% DP is a p (number of MVs) x j (number of LVs) matrix which associates each
% MV to one or more LVs putting value 1 in the appropriate cells.
DP=zeros(length(MV(1,:)),3); % 212 variables and 3 LVs
DP(1,1)=1; DP(2:length(MV(1,:))-1,2)=1; DP(length(MV(1,:)),3)=1;
%DP(1,1)=1;DP(2,2)=1; DP(3:length(MV(1,:))-1,3)=1; %DP(length(MV(1,:)),4)=1;% for a 4 variable model

color='blue';
path_graph(DB,xy,2,color,LV_labels)


% for each of 3 blocks, we have to define the reflective (A) or the
% formative (B) measurement model.
Mmode=['A'; 'A'; 'A']; %Mmode=['A'; 'A'; 'A'; 'A'];
method='path';
max_iter=300;

[results]=plssem(MV,DP,DB,Mmode,method,max_iter);

dec=3; % decimals of numbers in tables
%
[ OM_Tables, IM_Tables, OVERALL_Tables ] = summary_plssem(results,MV_labels,LV_labels,dec );




nboot=5000;
ci_type='studentized';
conf_level=0.95;
bias=1;
[ results ] = boot_confint( results,nboot,ci_type,conf_level,bias)
[ CI_Tables ] = summary_bootstrap(results,MV_labels,LV_labels,dec )

plssem_publish('html')
%https://www.cscu.cornell.edu/news/Handouts/SEM_fit.pdf
% for help with interpreting the GOF measures see above

%% visualize the weights
%https://stats.stackexchange.com/questions/235809/pls-partial-least-squares-weights-loadings-and-scores-interpretations
%need to look at weights

CI_Tables.W.t_stat(CI_Tables.W.Pvalue< .00135)
Tstats_plssem=zeros(1,210);
Tstats_plssem(ids)=CI_Tables.W.t_stat(2:length(CI_Tables.W.t_stat)-1)


component=1; 
ica2yeo7=readtable('D:\Canada_2020\UK_biobank\reports\ica2yeo7.csv');
nodes=dlmread('D:\Canada_2020\UK_biobank\reports\icad25_nodesordered.csv'); tmp=nanmean(X(clin.cdr==0,:))';
nodes=[Tstats_plssem',(1:210)', nodes]; nodes=sortrows(nodes, component);

myLabel=ica2yeo7.Yeo7N;figure
upperlim=1.8; lowerlim=-1.8; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=Tstats_plssem'; Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=Tstats_plssem'; Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

