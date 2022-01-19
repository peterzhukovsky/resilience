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

%% pls setup
%Y=[cognitive_final.tmt_cor,cognitive_final.WAIS ];%   clinical_final.Age
%ids=[0	1	0	0	0	0	1	0	1	0	1	0	1	0	1	1	0	1	0	1	1	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	1	0	1	0	1	1	1	0	1	1	0	1	0	1	1	1	0	1	1	1	0	1	0	1	1	1	0	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	0	1	0	1	1	1	0	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	1	1	1	1	1	1	0	1	0	1	0	1	1	1	0	1	1	1	1	1	1	0	1	1	0	1	0	1	1	1	0	1	1	1	1	1	1	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
%ids=[0	1	0	0	0	0	1	0	1	0	1	0	1	0	1	1	0	1	0	1	1	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	1	0	0	1	0	1	0	1	1	1	0	1	0	0	1	1	0	1	0	1	1	1	0	1	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	1	1	1	0	1	0	0	1	1	1	0	1	0	1	0	1	1	1	0	1	0	0	1	1	1	0	1	1	0	1	0	1	1	1	0	1	0	0	1	1	1	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0]
%ids=[0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];
%ids=[0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];

load('part_pls_m.mat', 'plsweights3')%load('permedu_part.mat', 'p_perm_edu'); %load('perm_part.mat', 'p_perm_inter')
ids=abs(plsweights3)>4 ;%ids=p_perm_edu<0.05
%ids=zeros(1,210); ids(21)=1; ids(150)=1; ids(177)=1; 
%ids=cellfun(@isempty,  strfind(ICs_vector, 'IC4')) & cellfun(@isempty,  strfind(ICs_vector, 'IC8')) & cellfun(@isempty,  strfind(ICs_vector, 'IC19')) & cellfun(@isempty,  strfind(ICs_vector, 'IC15'));

Y=[cognitive_final.MEMUNITS];%  %(isnan(Y)')' 
X=[amyloid_final.PET_fSUVR_rsf_TOT_CORTMEAN, part_icad25A(:,ids==1)];

%%% remove nans
naninx=sum(isnan(X)')'>0| isnan(Y)>0 |clinical_final.cdr>1.5   |~ix; Y=Y(naninx==0,:);  X=X(naninx==0,:); clin=clinical_final(naninx==0,:); mean_fd=mean_fd(naninx==0,:); prop_signal=prop_signal(naninx==0,:); %|clinical_final.Age<50 |mean_fd>0.1
%Y=zscore(Y); X=zscore(X);x=X;
%ncomp=3
% mediation analysis 
%Y=[cognitive_final.MEMUNITS , clinical_final.Age];%  %(isnan(Y)')' 
%X=part_icad25;

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
cd replication
plssem_publish('html')
cd ..
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
upperlim=1.9; lowerlim=-1.9; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=Tstats_plssem'; Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=Tstats_plssem'; Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

