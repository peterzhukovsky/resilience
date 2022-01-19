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
cognitiveA.aLOGIMEM(cognitiveA.Animal==0)=NaN;
%% pls setup
Y=[cognitiveA.tmt_cor,cognitiveA.WAIS ];%   clinicalA.Age
Y=[cognitiveA.MEMUNITS , clinicalA.Age];%  %(isnan(Y)')' 
X=part_icad25;

%%% remove nans
naninx=sum(isnan(X)')'>0| sum(isnan(Y)')'>0  |clinicalA.cdr>0.5 |~ix; Y=Y(naninx==0,:);  X=X(naninx==0,:); clin=clinicalA(naninx==0,:); mean_fd=mean_fd(naninx==0,:); prop_signal=prop_signal(naninx==0,:); 
Y=zscore(Y); X=zscore(X);x=X;
clinical_pls_bl=clinicalA(naninx==0,:); 
ncomp=3
%max(partialcorr(X,Y, [clin.Age,prop_signal,mean_fd, clin.delta_cog_t]))

%%% regress out sex and age and potentially other covariates
%mdl = fitlm(table( clin.M_F, clin.delta_cog_t, Y(:,1)) ); Y(:,1)=mdl.Residuals.Raw;
mdl = fitlm(table(clin.Age, clin.M_F, clin.delta_cog_t, Y(:,2)) ); Y(:,2)=mdl.Residuals.Raw;
clear x; for i=1:length(X(1,:)); mdl = fitlm(table(clin.M_F,clin.Age, mean_fd, X(:,i)) ); x(:,i)=mdl.Residuals.Raw;end

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,Y,ncomp);PCTVAR
%%% explore the PLS components and their correlation with Y
figure;imagesc(corr(XS,Y)); c=corr(XS, Y) %1.TMT  2.Tower  3.PAL  4.DSST
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
dim=3
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

ICs_vector(plsweights3>3)
ICs_vector(plsweights3<-3)


%% visuals for individual brain connectivities
rvalues=corr(X, Y(:,1), 'rows', 'pairwise')
figure; scatter(rvalues, plsweights3)

ids=1:210;figure;hold off
for i=ids(plsweights3<-3)
    j=find(i==ids(plsweights3<-3))
    allobservations=part_icad25(:,i);
    subplot(3,3,j); scatter(X(clin.cdr==0,i), Y(clin.cdr==0,1), 3,'filled', 'k'); hold on; scatter(X(clin.cdr==0.5,i), Y(clin.cdr==0.5,1), 3,'filled', 'k');
    m=fitlm(X(:,i),Y(:,1));
    [ypred,yci] = predict(m,X(:,i), 'Alpha',0.001); hold on
    yci=[yci,X(:,i)]; yci=sortrows(yci, 3);
    plot(X(:,i), ypred, 'k', 'LineWidth', 2);
    plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);
    plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5); 
    r=round(corr(X(:,i), Y(:,1), 'rows', 'pairwise'),2);%xlim([7 25]); 
    %text(8, 5, strcat('R=',num2str(r)) ); %ylim([4 14]);
    text(-4, min(Y(:,1)), strcat('R=',num2str(r)) ); %ylim([-1 1]);
    xlabel(ICs_vector{i}); clear r; 
end
amy=[-0.873425504603953;-0.930408169219439;-0.647972353299203;-0.743356378851212;-0.618242267412862;-0.701238757178896;0.385148131251134;2.32131997459907;-0.418802941258660;-0.733446350222432;-0.673986178449751;-0.643017338984813;-0.726013828750847;-0.564975863533169;-0.769370204001760;0.527604792789849;-0.755743914637188;-0.615764760255667;-0.849889186610600;-0.584795920790729;1.18414418944654;-0.939079444269621;-0.804055304202491;-0.722297568015055;-0.517903227546463;2.00296030489950;-0.900678083333098;-0.895723069018708;1.12096775693806;-0.765653943265968;3.19464124751032;-0.655404874770788;-0.761937682530175;-0.0285955640004398;0.568483660883568;-0.669031164135361;-0.894484315440110;-0.402699144736893;1.76635837138738;-0.730968843065237;0.326926713057050;-0.732207596643835;0.781549276402342;-0.862276722396575;-0.610809745941277;-0.989868340992120;-0.116547068080864;0.315777930849672;0.264989034127174;1.40588108001550;1.03673251359343;0.797653072924110;-0.739640118115420;-0.726013828750847;-0.618242267412862;1.06398509232258;-0.610809745941277;0.219155151719065;0.304629148642295;-0.716103800122067;1.08008888884435;-0.660359889085178;2.68551352670674;0.479293403224546;-0.962615762262974;2.61118831199089;-0.391550362529515;-0.697522496443104;2.14541696643822;-0.628152296041643;-0.272630018984153;-0.380401580322137;-0.827591622195844;-0.875903011761148;-0.800339043466699;-0.673986178449751;-0.537723284804023;2.62481460135546;1.46781875894537;-0.0570868963081828;1.19405421807532;-0.507993198917683;1.66973559225677;-0.448533027145001;0.243930223291016;0.676255222221552;0.250123991184004;-0.797861536309504;-0.915543126276268;-0.724775075172250;-0.680179946342739;-0.764415189687370;1.48640006262433;1.23617183974763;-1.12241497390206;-0.525335749018048;-0.192111036375313;0.386386884829731;-1.16453259557437;-0.815204086409869;0.123771126167056;-0.702477510757494;-0.0521318819937927;-0.525335749018048;-0.673986178449751;-0.593467195840912;2.74373494490082;-0.483218127345732;1.57682907386195;-0.629391049620240;0.943825995198618;-0.766892696844565;-0.764415189687370;0.915334662890875;-0.421280448415856;0.792698058609720;0.735715393994234;0.229065180347846;-0.101682025137694;0.727044118944051;-0.687612467814324;1.13335529272404;0.725805365365454;1.16556288576758;-0.690089974971519;1.49011632336013;-0.680179946342739;-0.776802725473346;-0.714865046543469;1.23245557901184;-0.760698928951578;0.810040608710085;-0.474546852295549;-0.779280232630541;0.985943616870934;0.372760595465159;0.744386669044417;-0.722297568015055;0.824905651653256;-0.647972353299203;-0.411370419787075;-0.132650864602632;-0.766892696844565;0.147307444160409;1.45419246958080;2.07233050530096;-0.491889402395914;0.733237886837039;1.22873931827605;1.70565944603610;-0.771847711158956;-0.703716264336092;-0.671508671292556;2.04012291225743;-0.706193771493286;0.0271483470364487;1.55453150944720;-0.749550146744200;-0.577363399319144;-0.356865262328784;-0.612048499519874;-0.796622782730906;-0.259003729619580;-0.879619272496940;2.49102921486693;0.0630722008157772;-0.761937682530175;-0.202021065004093;-0.755743914637188;-0.753266407479993;-0.763176436108773;-0.0992045179804987;-0.490650648817317;1.49754884483171;-0.646733599720605;0.681210236535943;1.85678738262500;-0.790429014837918;1.34270464750702;1.35509218329300;-0.737162610958225;-0.753266407479993;1.48640006262433;-0.716103800122067;-0.0831007214587310;1.78370092148774;0.440892042288023;0.126248633324251;-0.698761250021701;-0.903155590490293;-0.862276722396575;0.225348919612053;-0.926691908483646;0.0692659687087647;-0.843695418717612;-0.698761250021701;0.842248201753621;0.921528430783863;-1.04932851276480;-0.598422210155302;-0.380401580322137;-0.405176651894088;1.57559032028336;2.92707047453326;-0.595944702998107;0.871978287639961;-0.396505376843905;-0.617003513834265;-0.722297568015055;0.895514605633315;-0.794145275573711;0.509023489110887;-0.686373714235726;-0.768131450423163;0.624227571920457;-1.00968839824968;-0.753266407479993;-0.884574286811330;-0.883335533232733;1.92615758302646;-0.470830591559757;-0.626913542463045;-0.712387539386274;-0.418802941258660;-0.666553656978166;2.05746546235779;-0.920498140590659;-0.603377224469692;-0.831307882931637;-0.815204086409869;-0.739640118115420;-0.209453586475679;-0.657882381927983;0.833576926703439;-0.0583256498867805;-0.686373714235726;1.00824118128569;1.79361095011652;-0.671508671292556;1.17918917513215;-0.873425504603953;2.03516789794304;-0.324657669285248;0.948781009513008;-0.135128371759827;0.0407746364010217;1.62018544911287;-0.698761250021701;-0.0608031570439755;1.36871847265757;0.153501212053397;-0.745833886008407;-0.768131450423163;-0.764415189687370;-0.870947997446757;0.409923202823085;-0.634346063934630;-0.732207596643835;-0.681418699921336;1.06522384590118;0.658912672121187;-0.517903227546463;-0.654166121192190;0.499113460482106;-0.865992983132367;-0.448533027145001;-0.811487825674076;-0.755743914637188;-0.558782095640181;-0.273868772562750;-0.650449860456398;0.858351998275389;-0.0236405496860497;1.27952821499855;1.73910579265823;0.357895552521988;-0.297405090556103;-0.634346063934630;3.15376237941660;-0.708671278650482;2.24575600630462;0.388864391986926;-0.691328728550116;-0.701238757178896;-0.807771564938284;-0.574885892161949;-0.718581307279262;-0.696283742864506;1.14945908924581;-0.747072639587005;-0.811487825674076;2.66940973018497;-0.839979157981819;0.172082515732360;-0.693806235707311;-0.719820060857859;0.0135220576718761;-0.851127940189197;2.98653064630594;-0.730968843065237;-0.257764976040982;-0.714865046543469;-0.872186751025355;-0.986152080256327;-0.838740404403222;-0.589750935105119;-0.738401364536822;-0.708671278650482;1.86174239693939];
ids=1:210;figure;hold off
for i=ids(plsweights3<-3)
    j=find(i==ids(plsweights3<-3))
    allobservations=part_icad25(:,i);
    subplot(3,3,j); scatter(amy, X(:,i), 3,'filled', 'k'); hold on; 
    m=fitlm(amy, X(:,i));
    [ypred,yci] = predict(m,amy, 'Alpha',0.001); hold on
    yci=[yci,amy]; yci=sortrows(yci, 3);
    plot(amy, ypred, 'k', 'LineWidth', 2);
    plot(yci(:,3), yci(:,1), 'k', 'LineWidth', 0.5);
    plot(yci(:,3), yci(:,2), 'k', 'LineWidth', 0.5); 
    r=round(corr(X(:,i), amy, 'rows', 'pairwise'),2)%xlim([7 25]); 
    text(0, min(Y(:,1)), strcat('R=',num2str(r)) ); %ylim([-1 1]);
end

%% visualize the weights using circularGraph(x);
component=1; 
ica2yeo7=readtable('D:\Canada_2020\UK_biobank\reports\ica2yeo7.csv');
nodes=dlmread('D:\Canada_2020\UK_biobank\reports\icad25_nodesordered.csv'); tmp=nanmean(X(clin.cdr==0,:))';
nodes=[plsweights3,(1:210)',tmp, nodes]; nodes=sortrows(nodes, component);

myLabel=ica2yeo7.Yeo7N;figure
upperlim=3; lowerlim=-3; %upperlim=nodes(205,component); lowerlim=nodes(5,component);
Weights=plsweights3; Weights(Weights<upperlim & Weights>lowerlim)=0;%Weights=stats.W(:,component); 
Weights(Weights<0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);
Weights=plsweights3; Weights(Weights<upperlim & Weights>lowerlim)=0;
Weights(Weights>0)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

%%% average connectivity visual
figure; 
Weights=tmp;
Weights(Weights>-0.04)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,3)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

Weights=tmp;
Weights(Weights<0.03)=0; Weights_square= zeros(21); Weights_square(triu(ones(21),1)>0) = abs(Weights);
hold on; myColorMap=zeros(21,3);myColorMap(:,1)=1; circularGraph(Weights_square, 'Colormap',myColorMap, 'Label',myLabel);

%% compare BL vs matched
load('clinical_pls.mat');
for i=1:length(clinical_pls_bl.Subject)
idx=strcmp(clinical_pls_bl.Subject{i}, clinical_pls.Subject);
if sum(idx)>0
tmp(i)=sum(clinical_pls_bl.time_d(i)==clinical_pls.time_d(idx));
else
    tmp(i)=0;
end
end
sum(tmp)/length(tmp)


for i=1:length(clinical_pls.Subject)
idx=strcmp(clinical_pls.Subject{i}, clinical_pls_bl.Subject);
tmp(i)=sum(idx);
end



for i=1:length(clinical_pls.Subject)
idx=strcmp(clinical_pls.Subject{i}, clinical_pls_bl.Subject);
if sum(idx)>0
tmp(i)=sum(clinical_pls.time_d(i)==clinical_pls_bl.time_d(idx));
else
    tmp(i)=0;
end
end
sum(tmp)/length(tmp)

mean(clinical_pls_bl.Age)
mean(clinical_pls.Age)
