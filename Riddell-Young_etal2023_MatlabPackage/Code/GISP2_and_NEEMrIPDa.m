%% Script for calculating the methane interpolar difference using the WAIS Divide, NEEM and GISP2 ice cores

% This script was originally written by Julia Rosen and was modified by Ben Riddell-Young
% The output of this code is published in Riddell-Young et al., 2023.
% "Tropical Sources mainly controlled last glacial maximum and deglacial
% methane variability." Nature Geoscience.

% This script uses discrete CH4 records from each core and published and
% unpublished tie-points to calculate the methane interpolar difference. It
% uses a monte-carlo analysis to propagate analytical and data processing
% uncertainties. 


%% Start with corrected NEEM record

tic
clear all
%clf
format long

%Number of Monte Carlo iterations:
nt=1000;

%Create cell empty array to save all data for figure plotting
gradient = cell(8,5);

for un = [2 8]
%Load WAIS chronology. 
WAIS_Chronology = csvread('WAIS_Chronology_Buizert2015.csv',1,0);
WAISdepth_Chronology = WAIS_Chronology(:,1);
WAISiceage_Chronology = WAIS_Chronology(:,2);
WAISgasage_Chronology = WAIS_Chronology(:,3);

%Load data:
NEEMWAIS_Data=csvread('WAISfinal_NEEM_may22.csv',1,0);
%Define variables:
NEEMWAIS_Depth=NEEMWAIS_Data(1:(end-1),1);
NEEMWAIS_CH4=NEEMWAIS_Data(1:(end-1),2);
NEEMWAIS_Age = interp1(WAISdepth_Chronology,WAISgasage_Chronology,NEEMWAIS_Depth);
%Average duplicates 
     [C2,ia2,idx2] = unique(NEEMWAIS_Depth,'stable');
     NEEMWAIS_Depth = accumarray(idx2,NEEMWAIS_Depth,[],@mean); 
     NEEMWAIS_CH4 = accumarray(idx2,NEEMWAIS_CH4,[],@mean); 
     NEEMWAIS_Age = accumarray(idx2,NEEMWAIS_Age,[],@mean); 

%Create evenly-spaced age scale for interpolation:
dt=50;
NEEMInterp_Age=[ceil(min(NEEMWAIS_Age)):dt:floor(max(NEEMWAIS_Age))]';

%Define initial tie points:
            %NEEM depth %WAIS age
NEEMStartTies= [1558.00 2480.375;          
            1540.88 2441.831;          
            1536.77	2432.261851;
            1533.21 2420.578218;
            1525.37 2382.548527;
            1523.49 2370.955072;
            1515.04 2323.285244;
            1513.09	2310.066215;
            1511.56 2300.893482;
            1509.29 2286.051225;
            1504.82 2259.764413;
            1496.74 2234.089041;
            1493.83 2227.404703;
            1489.49 2214.565606;
            1487.98 2207.309825;
            1483.96 2191.487283;
            1482.86 2185.604679;
            1478.41 2169.272781;
            1476.57 2162.709564;
            1471.76 2147.22058;
            1469.15 2141.862982;
            1466.13 2128.187257;
            1461.82 2111.609573;
            1460.63 2104.354008;
            1455.77 2083.761553;
            1450.17 2057.608585
            1444.49 2028.86036;
            1441.79 2017.700511;
            1435.18 1983.600412; 
            1432.53 1973.051413;
            1428.89 1960.239512;
            1425.54 1951.335983];

%Define variables:
[n,m]=size(NEEMStartTies);
TieDepthNEEM=NEEMStartTies(:,1);
NEEMTieDepthWAIS=NEEMStartTies(:,2);
NEEMTieAge = interp1 (WAISdepth_Chronology,WAISgasage_Chronology,NEEMTieDepthWAIS);

%Load data:  
NEEM_Data=csvread('NEEMfinal_nov22.csv',1,0);
%Average duplicates
     [Depth,ia2,idx2] = unique(NEEM_Data(:,1),'stable');
     val2 = accumarray(idx2,NEEM_Data(:,un),[],@mean);           %un changes between 2 and 8 for uncorrected and corrected GISP2 data, respectively.
     val3 = accumarray(idx2,NEEM_Data(:,7),[],@mean);
     NEEM_Data2 = [Depth val2 val3];   

%Define variables:
NEEM_Depth = NEEM_Data2(:,1);
NEEM_CH4 = NEEM_Data2(:,2); 
NEEM_Ca2 = NEEM_Data2(:,3);

%Interpolate NEEM data onto WAIS timescale for initial estimate:
NEEM_Age_Start=interp1(TieDepthNEEM,NEEMTieAge,NEEM_Depth);
NEEM_CH4_Start=interp1(NEEM_Age_Start,NEEM_CH4,NEEMWAIS_Age);

%Calculate starting R-value:
NEEMRstart=corrcoef(NEEMWAIS_CH4,NEEM_CH4_Start);

%Initialize solution arrays:
BestCorr=zeros(length(NEEM_CH4),1);
BestSmooth=zeros(length(NEEM_CH4),1);
NEEMR=zeros(nt,1);
Smoothness=zeros(nt,1);
BestIndex=zeros(nt,1);
%%
%Monte Carlo synchronization algorithm:
%Uncertainty in CH4xs correction (applied to each Monte Carlo Simulation,
%rather than each sample).
CH4xs_u = normrnd(0,0.0012,[1000,1]);
for i=1:nt   
    
    %Add analytical error to data (including error from CH4xs correction). 

    %Define constants for CH4xs uncertainty
    Pstp = 101.325;                                                             %Kpa
    Tstp = 273;                                                                 %Degrees K
    mCalcium = 40.078;                                                          %Molecular Weight of elemental calcium
    R = 8.314463;                                                               %Ideal Gas Constant (J/(mol*K) & L*kPa/(K*mol)))
    Cxs_unc = .0012;                                                            %standard deviation of regressions from binning experiments.
    mTACn = .092;                                                               %Rough mean of TAC. Does not impact uncertainty to calculate rough mean.
    %Calculatation of NEEM stdev 
       NEEM_stdev = 4;
    NEEM_CH4_Error = NEEM_CH4+normrnd(0,NEEM_stdev,[length(NEEM_CH4),1])+(NEEM_Ca2*R*Tstp*CH4xs_u(i)/(mCalcium*mTACn*Pstp));
    WAIS_stdev = 2;
    NEEMWAIS_CH4_Error = NEEMWAIS_CH4+normrnd(0,WAIS_stdev,[length(NEEMWAIS_CH4),1]);
    
    %Create evenly-spaced chronological anchors:   
    NEEMnties=50;
    NEEMTieAgeTemp=linspace(NEEMTieAge(1),NEEMTieAge(end),NEEMnties)';
    NEEMTieDepthTemp=interp1(NEEMTieAge,TieDepthNEEM,NEEMTieAgeTemp);
    
    %Define tie point age uncertainty:
    NEEMTieFix=25;
    NEEMTieError=[NEEMTieFix*randn(NEEMnties-2,1)];
    
    %Add uncertainty to all points except first and last:
    NEEMTieAgeError=[NEEMTieAgeTemp(1); NEEMTieAgeTemp(2:end-1)+NEEMTieError;...
        NEEMTieAgeTemp(end)];

    %Make sure tie points are still in decreasing order after adding 
    %age uncertainty:
    NEEMTest=NEEMTieAgeError(2:end)-NEEMTieAgeError(1:end-1);
    
    %If not, redo tie points until they all decrease:
    while max(NEEMTest)>0   
        NEEMTieError=[NEEMTieFix*randn(NEEMnties-2,1)];
        NEEMTieAgeError=[NEEMTieAgeTemp(1); NEEMTieAgeTemp(2:end-1)+NEEMTieError;...
        NEEMTieAgeTemp(end)];

        %Check again:
        NEEMTest=NEEMTieAgeError(2:end)-NEEMTieAgeError(1:end-1);
    end
    
    %Define tie points for this iteration:
    NEEMTieDepthTotal=NEEMTieDepthTemp;
    NEEMTieAgeTotal=NEEMTieAgeError;

    %Interpolate NEEM data with error onto new chronology:
    NEEM_Age=interp1(NEEMTieDepthTotal,NEEMTieAgeTotal,NEEM_Depth);
    NEEM_Interp_CH4=interp1(NEEM_Age,NEEM_CH4_Error,NEEMInterp_Age);
    
    %Interpolate WAIS data onto evenly-spaced chronology:
    NEEMWAIS_Interp_CH4=interp1(NEEMWAIS_Age,NEEMWAIS_CH4_Error,NEEMInterp_Age);

    %Calculate correlation between records:
    NEEMRs=corrcoef(NEEM_Interp_CH4,NEEMWAIS_Interp_CH4);
    NEEMR(i)=NEEMRs(2,1);

    %Calculate gas age annual layer thickness as a check:
    ALT=(NEEM_Depth(2:end)-NEEM_Depth(1:end-1))./...
        (NEEM_Age(2:end)-NEEM_Age(1:end-1));

    %Save NEEM age scale that results in best fit:
    if NEEMR(i)>=max(NEEMR(1:i-1))
        NEEMBestAll=NEEM_Age;
        NEEMRBestAll=NEEMR(i);
        NEEMTieBestAll=[NEEMTieDepthTotal NEEMTieAgeTotal];
        NEEMALTAll=ALT;
    end
    %%
    
    %Loop through different amounts of smoothing to test sensitivity:
    for j=5:2:15
    
        %Generate gaussian smoothing filter with width = w*dt (dt=50):
        w=j;
        window=gausswin(w);
        window=window/sum(window);
        
        %Smooth both sets of data with filter:
        NEEM_Smooth_CH4=conv(NEEM_Interp_CH4,window,'same');
        NEEMWAIS_Smooth_CH4=conv(NEEMWAIS_Interp_CH4,window,'same');
    
        %Interpolate both records onto WAIS age scale:
        NEEM_Smooth_onWAIS=interp1(NEEMInterp_Age,NEEM_Smooth_CH4,...
            NEEMWAIS_Age);
        NEEMWAIS_Smooth_onWAIS=interp1(NEEMInterp_Age,NEEMWAIS_Smooth_CH4,...
            NEEMWAIS_Age);
        
        %Transfer variable names for next calculations:
        NEEM_Calc=NEEM_Interp_CH4;
        NEEMWAIS_Calc=NEEMWAIS_Interp_CH4;
        NEEMCalc_Age=NEEMInterp_Age;
    
        %Calculate estimates of the IPD for this iteration and degree of smoothing:
        NEEMIPGAbs(:,i,j)=NEEM_Smooth_CH4-NEEMWAIS_Smooth_CH4;
        NEEMIPG(:,i,j)=((NEEM_Smooth_CH4./NEEMWAIS_Smooth_CH4)-1)*100;
        
        %Calculate interval-average variables:
        NEEMOD2_NH=mean(NEEM_Calc(find((NEEMCalc_Age<20415)&...
            (NEEMCalc_Age>17803))));
        NEEMOD2_NH_Std=std(NEEM_Calc(find((NEEMCalc_Age<20415)&...
            (NEEMCalc_Age>17803))));
        NEEMOD2_SH=mean(NEEMWAIS_Calc(find((NEEMCalc_Age<20415)&...
            (NEEMCalc_Age>17803))));
        NEEMOD2_SH_Std=std(NEEMWAIS_Calc(find((NEEMCalc_Age<20415)&...
            (NEEMCalc_Age>17803))));
             
        NEEMOD2_IPG_Abs(i,j)=NEEMOD2_NH-NEEMOD2_SH;
        NEEMOD2_IPG_Abs_Std(i,j)=sqrt(NEEMOD2_NH_Std^2+NEEMOD2_SH_Std^2);
        NEEMOD2_rIPG(i,j)=((NEEMOD2_NH/NEEMOD2_SH)-1)*100;
        
        NEEMHS1_NH=mean(NEEM_Calc(find((NEEMCalc_Age<17803)&...
            (NEEMCalc_Age>14705))));
        NEEMHS1_NH_Std=std(NEEM_Calc(find((NEEMCalc_Age<17803)&...
            (NEEMCalc_Age>14705))));
        NEEMHS1_SH=mean(NEEMWAIS_Calc(find((NEEMCalc_Age<17803)&...
            (NEEMCalc_Age>14705))));
        NEEMHS1_SH_Std=std(NEEMWAIS_Calc(find((NEEMCalc_Age<17803)&...
            (NEEMCalc_Age>14705))));
        
        NEEMOD1_IPG_Abs(i,j)=NEEMHS1_NH-NEEMHS1_SH;
        NEEMOD1_IPG_Abs_Std(i,j)=sqrt(NEEMHS1_NH_Std^2+NEEMHS1_SH_Std^2);
        NEEMOD1_rIPG(i,j)=((NEEMHS1_NH/NEEMHS1_SH)-1)*100;
        
        NEEMBA2_NH=mean(NEEM_Calc(find((NEEMCalc_Age<14353)&...
            (NEEMCalc_Age>13735))));
        NEEMBA2_NH_Std=std(NEEM_Calc(find((NEEMCalc_Age<14353)&...
            (NEEMCalc_Age>13735))));
        NEEMBA2_SH=mean(NEEMWAIS_Calc(find((NEEMCalc_Age<14353)&...
            (NEEMCalc_Age>13735))));
        NEEMBA2_SH_Std=std(NEEMWAIS_Calc(find((NEEMCalc_Age<14353)&...
            (NEEMCalc_Age>13735))));

        NEEMBA2_IPG_Abs(i,j)=(NEEMBA2_NH-NEEMBA2_SH);
        NEEMBA2_IPG_Abs_Std(i,j)=sqrt(NEEMBA2_NH_Std^2+NEEMBA2_SH_Std^2);
        NEEMBA2_rIPG(i,j)=((NEEMBA2_NH/NEEMBA2_SH)-1)*100;
        
        NEEMBA1_NH=mean(NEEM_Calc(find((NEEMCalc_Age<13735)&...
            (NEEMCalc_Age>12943))));
        NEEMBA1_NH_Std=std(NEEM_Calc(find((NEEMCalc_Age<13735)&...
            (NEEMCalc_Age>12943))));
        NEEMBA1_SH=mean(NEEMWAIS_Calc(find((NEEMCalc_Age<13735)&...
            (NEEMCalc_Age>12943))));
        NEEMBA1_SH_Std=std(NEEMWAIS_Calc(find((NEEMCalc_Age<13735)&...
            (NEEMCalc_Age>12943))));

        NEEMBA1_IPG_Abs(i,j)=(NEEMBA1_NH-NEEMBA1_SH);
        NEEMBA1_IPG_Abs_Std(i,j)=sqrt(NEEMBA1_NH_Std^2+NEEMBA1_SH_Std^2);
        NEEMBA1_rIPG(i,j)=((NEEMBA1_NH/NEEMBA1_SH)-1)*100;
        
        NEEMYD_NH=mean(NEEM_Calc(find((NEEMCalc_Age<12449)&...
            (NEEMCalc_Age>11707))));
        NEEMYD_NH_Std=std(NEEM_Calc(find((NEEMCalc_Age<12449)&...
            (NEEMCalc_Age>11707))));
        NEEMYD_SH=mean(NEEMWAIS_Calc(find((NEEMCalc_Age<12449)&...
            (NEEMCalc_Age>11707))));
        NEEMYD_SH_Std=std(NEEMWAIS_Calc(find((NEEMCalc_Age<12449)&...
            (NEEMCalc_Age>11707))));

        NEEMYD_IPG_Abs(i,j)=(NEEMYD_NH-NEEMYD_SH);
        NEEMYD_IPG_Abs_Std(i,j)=sqrt(NEEMYD_NH_Std^2+NEEMYD_SH_Std^2);
        NEEMYD_rIPG(i,j)=((NEEMYD_NH/NEEMYD_SH)-1)*100;
        
    end
    
end
%%{
%Save best NEEM chronology:
NEEM_results=[NEEM_Depth NEEMBestAll NEEM_CH4];
save 'NEEM_Chronology.txt' NEEM_results -ascii

%Average over nt Monte Carlo iterations:
[nn,mm]=size(NEEMIPG);
for j=5:2:15    %loop through different degrees of smoothing
    for k=1:nn  %loop through time steps
        NEEMUpperBound(k,j)=max(NEEMIPG(k,:,j));
        NEEMLowerBound(k,j)=min(NEEMIPG(k,:,j));
        NEEMIPGMean(k,j)=mean(NEEMIPG(k,:,j));
        NEEMIPGStd(k,j)=std(NEEMIPG(k,:,j));
        NEEMIPGAbsMean(k,j)=nanmean(NEEMIPGAbs(k,:,j));
    end
end
%%

%Clear variables for next calculations:
clear NEEM_Calc WAIS_Calc

%Use best-fit chronology to save concentration and IPD results:

%Loop through different degrees of smoothing:
for j=5:2:15
    
    %Generate gaussian smoothing filter with width = w*dt (dt=50):
    w=j;
    window=gausswin(w);
    window=window/sum(window);
    
    %Interpolate both records onto evenly-spaced chronology:
    NEEM_Interp_CH4=interp1(NEEMBestAll,NEEM_CH4,NEEMInterp_Age);
    NEEM_Smooth_CH4=conv(NEEM_Interp_CH4,window,'same');

    NEEMWAIS_Interp_CH4=interp1(NEEMWAIS_Age,NEEMWAIS_CH4,NEEMInterp_Age);
    NEEMWAIS_Smooth_CH4=conv(NEEMWAIS_Interp_CH4,window,'same');
    
    %Save smoothed records:
    NEEM_Calc(:,j)=NEEM_Interp_CH4';
    NEEMWAIS_Calc(:,j)=NEEMWAIS_Interp_CH4';
    NEEMCalc_Age=NEEMInterp_Age;
end

%Choose amount of smoothing to save final results:
saveindex=11;
NEEMIPG_results=[NEEMCalc_Age NEEM_Calc(:,saveindex) NEEMWAIS_Calc(:,saveindex)...
    NEEMIPGAbsMean(:,saveindex) NEEMIPGMean(:,saveindex) NEEMIPGStd(:,saveindex)];
%save IPG_cont_results_onWAIS.txt IPG_results -ascii
%%
NEEM_Calc = NEEM_Calc(:,9);
NEEMWAIS_Calc = NEEMWAIS_Calc(:,9);
%Calculate interval-averages
%--------------------------------------------
NEEMHS1_NH=mean(NEEM_Calc(find((NEEMInterp_Age<17803)&(NEEMInterp_Age>14705))));
NEEMHS1_NH_Std=std(NEEM_Calc(find((NEEMInterp_Age<17803)&(NEEMInterp_Age>14705))));
NEEMHS1_SH=mean(NEEMWAIS_Calc(find((NEEMInterp_Age<17803)&(NEEMInterp_Age>14705))));
NEEMHS1_SH_Std=std(NEEMWAIS_Calc(find((NEEMInterp_Age<17803)&(NEEMInterp_Age>14705))));
        
NEEMHS1_Age=mean(NEEMInterp_Age(find((NEEMInterp_Age<17803)&(NEEMInterp_Age>14705))));
NEEMHS1_Dur=17803-14705;
NEEMHS1_IPG=mean(NEEMIPGMean(find((NEEMInterp_Age<17803)&(NEEMInterp_Age>14705)),:));
NEEMHS1_IPG_Std1=std(NEEMIPGMean(find((NEEMInterp_Age<17803)&(NEEMInterp_Age>14705)),:));
NEEMHS1_IPG_Std=mean(NEEMIPGStd(find((NEEMInterp_Age<17803)&(NEEMInterp_Age>14705)),:));
%OD1_IPG_Std = OD1_IPG_Std1 + OD1_IPG_Std2;
NEEMHS1_rIPG_Avg=mean(NEEMOD1_rIPG);
NEEMHS1_rIPG_Std=std(NEEMOD1_rIPG);
%--------------------------------------------
NEEMBA2_NH=mean(NEEM_Calc(find((NEEMInterp_Age<14353)&(NEEMInterp_Age>13735))));
NEEMBA2_NH_Std=std(NEEM_Calc(find((NEEMInterp_Age<14353)&(NEEMInterp_Age>13735))));
NEEMBA2_SH=mean(NEEMWAIS_Calc(find((NEEMInterp_Age<14353)&(NEEMInterp_Age>13735))));
NEEMBA2_SH_Std=std(NEEMWAIS_Calc(find((NEEMInterp_Age<14353)&(NEEMInterp_Age>13735))));
        
NEEMBA2_Age=mean(NEEMInterp_Age(find((NEEMInterp_Age<14353)&(NEEMInterp_Age>13735))));
NEEMBA2_Dur=14353-13735;
NEEMBA2_IPG=mean(NEEMIPGMean(find((NEEMInterp_Age<14353)&(NEEMInterp_Age>13735)),:));
NEEMBA2_IPG_Std1=std(NEEMIPGMean(find((NEEMInterp_Age<14353)&(NEEMInterp_Age>13735)),:));
NEEMBA2_IPG_Std=mean(NEEMIPGStd(find((NEEMInterp_Age<14353)&(NEEMInterp_Age>13735)),:));
%BA2_IPG_Std = BA2_IPG_Std1 + BA2_IPG_Std2;
NEEMBA2_rIPG_Avg=mean(NEEMBA2_rIPG);
NEEMBA2_rIPG_Std=std(NEEMBA2_rIPG);
%--------------------------------------------
NEEMBA1_NH=mean(NEEM_Calc(find((NEEMInterp_Age<13735)&(NEEMInterp_Age>12943))));
NEEMBA1_NH_Std=std(NEEM_Calc(find((NEEMInterp_Age<13735)&(NEEMInterp_Age>12943))));
NEEMBA1_SH=mean(NEEMWAIS_Calc(find((NEEMInterp_Age<13735)&(NEEMInterp_Age>12943))));
NEEMBA1_SH_Std=std(NEEMWAIS_Calc(find((NEEMInterp_Age<13735)&(NEEMInterp_Age>12943))));
        
NEEMBA1_Age=mean(NEEMInterp_Age(find((NEEMInterp_Age<13735)&(NEEMInterp_Age>12943))));
NEEMBA1_Dur=13735-12943;
NEEMBA1_IPG=mean(NEEMIPGMean(find((NEEMInterp_Age<13735)&(NEEMInterp_Age>12943)),:));
NEEMBA1_IPG_Std1=std(NEEMIPGMean(find((NEEMInterp_Age<13735)&(NEEMInterp_Age>12943)),:));
NEEMBA1_IPG_Std=mean(NEEMIPGStd(find((NEEMInterp_Age<13735)&(NEEMInterp_Age>12943)),:));
%BA1_IPG_Std = BA1_IPG_Std1 + BA1_IPG_Std2;
NEEMBA1_rIPG_Avg=mean(NEEMBA1_rIPG);
NEEMBA1_rIPG_Std=std(NEEMBA1_rIPG);
%--------------------------------------------
NEEMYD_NH=mean(NEEM_Calc(find((NEEMInterp_Age<12449)&(NEEMInterp_Age>11707))));
NEEMYD_NH_Std=std(NEEM_Calc(find((NEEMInterp_Age<12449)&(NEEMInterp_Age>11707))));
NEEMYD_SH=mean(NEEMWAIS_Calc(find((NEEMInterp_Age<12449)&(NEEMInterp_Age>11707))));
NEEMYD_SH_Std=std(NEEMWAIS_Calc(find((NEEMInterp_Age<12449)&(NEEMInterp_Age>11707))));

NEEMYD_Age=mean(NEEMInterp_Age(find((NEEMInterp_Age<12449)&(NEEMInterp_Age>11707))));
NEEMYD_Dur=12449-11707;
NEEMYD_IPG=mean(NEEMIPGMean(find((NEEMInterp_Age<12449)&(NEEMInterp_Age>11707)),:));
NEEMYD_IPG_Std1=std(NEEMIPGMean(find((NEEMInterp_Age<12449)&(NEEMInterp_Age>11707)),:));
NEEMYD_IPG_Std=mean(NEEMIPGStd(find((NEEMInterp_Age<12449)&(NEEMInterp_Age>11707)),:));
%YD_IPG_Std =YD_IPG_Std1 +YD_IPG_Std2;
NEEMYD_rIPG_Avg=mean(NEEMYD_rIPG);
NEEMYD_rIPG_Std=std(NEEMYD_rIPG);
%--------------------------------------------

%Save results:
NEEMAll_IPG=   [NEEMYD_Age  NEEMYD_Dur  NEEMYD_IPG  NEEMYD_IPG_Std;
            NEEMBA1_Age NEEMBA1_Dur NEEMBA1_IPG NEEMBA1_IPG_Std;
            NEEMBA2_Age NEEMBA2_Dur NEEMBA2_IPG NEEMBA2_IPG_Std;
            NEEMHS1_Age NEEMHS1_Dur NEEMHS1_IPG NEEMHS1_IPG_Std];
%             OD2_Age OD2_Dur OD2_IPG OD2_IPG_Std];
%  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now measure Corrected GISP2 record

%Number of Monte Carlo iterations:
nt=1000;

%Load WAIS Divide chronology (Buizert, 2015).
WAIS_Chronology = csvread('WAIS_Chronology_Buizert2015.csv',1,0);
WAISdepth_Chronology = WAIS_Chronology(:,1);
WAISiceage_Chronology = WAIS_Chronology(:,2);
WAISgasage_Chronology = WAIS_Chronology(:,3);

%Load data:
GISP2WAIS_Data=csvread('WAISfinal_GISP2_may22.csv',1,0);
%Average duplicates
for i = 1:(length(GISP2WAIS_Data)-1)
    if GISP2WAIS_Data((i+1),1) == GISP2WAIS_Data(i,1)
        for j = 1:length(GISP2WAIS_Data(1,:))
            GISP2WAIS_Data((i+1),j) = mean(GISP2WAIS_Data(i:(i+1),j));
        end
        GISP2WAIS_Data(i,:) = [0];
    else 
        GISP2WAIS_Data(i,:) = GISP2WAIS_Data(i,:);
    end 
end
GISP2WAIS_Data( ~any(GISP2WAIS_Data,2), : ) = [];

%Define variables:
GISP2WAIS_Depth=GISP2WAIS_Data(:,1);
GISP2WAIS_CH4=GISP2WAIS_Data(:,2);
GISP2WAIS_gasage = interp1(WAISdepth_Chronology,WAISgasage_Chronology,GISP2WAIS_Depth);

%Create evenly-spaced age scale for interpolation:
dt=50;
GISP2Interp_Age=[ceil(min(GISP2WAIS_gasage)):dt:floor(max(GISP2WAIS_gasage))]';

%Define initial tie points:
            %GISP2 depth %WAIS Depth
GISP2StartTies= [2065.924465, 2747.303781;
            2063.614009, 2744.011371;
            2039.394898, 2683.96069;
            2033.600000, 2668.300000;   %Tie replaced by me
            2024.446274, 2653.333375;
            2017.507493, 2635.717801;
            2014.661786, 2632.093108;
            2012.441478, 2626.456989;
            2005.957060, 2612.27433;
            1996.936978, 2598.559959;
            1994.943383, 2596.329822;
            1992.567955, 2592.671161;
            1988.972424, 2586.284441;
            %1974.938735, 2569.613681;  %removed due to unrealistic delt age and low data resolution
            %1972.000000, 2557.500000;  %(added by me originally) removed due to unrealistic delt age and low data resolution
            1946.630914, 2532.627352;
            1938.747025, 2519.074601;
            1917.000000, 2485.300000;   %Tie added by me
            1910.800000, 2475.400000;   %Tie added by me
            1894.938836, 2458.468788;
            1883.360557, 2442.527591;
            1876.935177, 2432.261851;
            1872.901412, 2421.943375;
            1867.410429, 2411.173673;
            1864.471867, 2404.914792;
            1858.010677, 2387.646424;
            1854.746433, 2380.661192;
            1854.379307, 2377.215276;
            1852.248637, 2370.955072;
            1844.323202, 2345.763621;
            1832.152722, 2300.893482;
            1827.141493, 2286.051225;
            1822.850091, 2272.367043];

%Define variables:
[n,m]=size(GISP2StartTies);
GISP2TieDepthGISP2=GISP2StartTies(:,1);
GISP2TieDepthWAIS=GISP2StartTies(:,2);
GISP2TieAge = interp1 (WAISdepth_Chronology,WAISgasage_Chronology,GISP2TieDepthWAIS);

%Load data:
GISP2_Data=csvread('GISP2final_gradient_nov22.csv',1,0);
%Define variables:
GISP2_Depth = GISP2_Data(:,1);
GISP2_CH4 = GISP2_Data(:,un);                    %un changes between 2 and 8 for uncorrected and corrected GISP2 data, respectively.
GISP2_Ca2 = GISP2_Data(:,7);
%Average duplicates 
     [C2,ia2,idx2] = unique(GISP2_Depth,'stable');
     GISP2_Depth = accumarray(idx2,GISP2_Depth,[],@mean); 
     GISP2_CH4 = accumarray(idx2,GISP2_CH4,[],@mean); 
     GISP2_Ca2 = accumarray(idx2,GISP2_Ca2,[],@mean); 

%Interpolate GISP2 data onto WAIS timescale for initial estimate:
GISP2_Age_Start=interp1(GISP2TieDepthGISP2,GISP2TieAge,GISP2_Depth);
GISP2_CH4_Start=interp1(GISP2_Age_Start,GISP2_CH4,GISP2WAIS_gasage);

%Calculate starting R-value:
GISP2Rstart=corrcoef(GISP2WAIS_CH4,GISP2_CH4_Start);

%Initialize solution arrays:
GISP2BestCorr=zeros(length(GISP2_CH4),1);
GISP2BestSmooth=zeros(length(GISP2_CH4),1);
GISP2R=zeros(nt,1);
GISP2Smoothness=zeros(nt,1);
GISP2BestIndex=zeros(nt,1);

%%
%Monte Carlo synchronization algorithm:
%Uncertainty in CH4xs correction (applied to each Monte Carlo Simulation,
%rather than each sample).
for i=1:nt   
    
    %Add analytical error to data (including error from CH4xs correction). 
    %.00573 subject to change depending on the std error of the CH4xs Ca2+
    %regression relationship (Currently based off of NEEM correction
    %sensitivity, though this shouldn't change).
    mTACg = .0945;
    %Calculation of GISP2 stdev 
    GISP2_stdev = 3;
    GISP2_CH4_Error = GISP2_CH4+normrnd(0,GISP2_stdev,[length(GISP2_CH4),1])+(GISP2_Ca2*R*Tstp*CH4xs_u(i)/(mCalcium*mTACg*Pstp));
    GISP2WAIS_CH4_Error = GISP2WAIS_CH4+normrnd(0,WAIS_stdev,[length(GISP2WAIS_CH4),1]);
    
    %Create evenly-spaced chronological anchors:   
    GISP2nties=60;%length(GISP2StartTies);
    GISP2TieAgeTemp=linspace(GISP2TieAge(1),GISP2TieAge(end),GISP2nties)';
    GISP2TieDepthTemp=interp1(GISP2TieAge,GISP2TieDepthGISP2,GISP2TieAgeTemp);
    
    %Define tie point age uncertainty:
    GISP2TieFix=25;
    GISP2TieError=[GISP2TieFix*randn(GISP2nties-2,1)];
    
    %Add uncertainty to all points except first and last:
    GISP2TieAgeError=[GISP2TieAgeTemp(1); GISP2TieAgeTemp(2:end-1)+GISP2TieError;...
        GISP2TieAgeTemp(end)];

    %Make sure tie points are still in decreasing order after adding age uncertainty:
    GISP2Test=GISP2TieAgeError(2:end)-GISP2TieAgeError(1:end-1);
    
    %If not, redo tie points until they all decrease:
    while max(GISP2Test)>0   
        GISP2TieError=[GISP2TieFix*randn(GISP2nties-2,1)];
        GISP2TieAgeError=[GISP2TieAgeTemp(1); GISP2TieAgeTemp(2:end-1)+GISP2TieError;...
        GISP2TieAgeTemp(end)];

        %Check again:
        GISP2Test=GISP2TieAgeError(2:end)-GISP2TieAgeError(1:end-1);
    end
    
    %Define tie points for this iteration:
    GISP2TieDepthTotal=GISP2TieDepthTemp;
    GISP2TieAgeTotal=GISP2TieAgeError;

    %Interpolate GISP2 data with error onto new chronology:
    GISP2_Age=interp1(GISP2TieDepthTotal,GISP2TieAgeTotal,GISP2_Depth);
    GISP2_Interp_CH4=interp1(GISP2_Age,GISP2_CH4_Error,GISP2Interp_Age);
    
    %Interpolate WAIS data onto evenly-spaced chronology:
    GISP2WAIS_Interp_CH4=interp1(GISP2WAIS_gasage,GISP2WAIS_CH4_Error,GISP2Interp_Age);

    %Calculate correlation between records:
    GISP2Rs=corrcoef(GISP2_Interp_CH4,GISP2WAIS_Interp_CH4);
    GISP2R(i)=GISP2Rs(2,1);

    %Calculate gas age annual layer thickness as a check:
    GISP2ALT=(GISP2_Depth(2:end)-GISP2_Depth(1:end-1))./...
        (GISP2_Age(2:end)-GISP2_Age(1:end-1));

    %Save GISP2 age scale that results in best fit:
    if GISP2R(i)>=max(GISP2R(1:i-1))
        GISP2BestAll=GISP2_Age;
        GISP2RBestAll=GISP2R(i);
        GISP2TieBestAll=[GISP2TieDepthTotal GISP2TieAgeTotal];
        GISP2ALTAll=GISP2ALT;
    end
    %%
    
    %Loop through different amounts of smoothing to test sensitivity:
    for j=5:2:15
    
        %Generate gaussian smoothing filter with width = w*dt (dt=50):
        w=j;
        window=gausswin(w);
        window=window/sum(window);
        
        %Smooth both sets of data with filter:
        GISP2_Smooth_CH4=conv(GISP2_Interp_CH4,window,'same');
        GISP2WAIS_Smooth_CH4=conv(GISP2WAIS_Interp_CH4,window,'same');
    
        %Interpolate both records onto WAIS age scale:
        GISP2_Smooth_onWAIS=interp1(GISP2Interp_Age,GISP2_Smooth_CH4,...
            GISP2WAIS_gasage);
        GISP2WAIS_Smooth_onWAIS=interp1(GISP2Interp_Age,GISP2WAIS_Smooth_CH4,...
            GISP2WAIS_gasage);
        
        %Transfer variable names for next calculations:
        GISP2_Calc=GISP2_Interp_CH4;
        GISP2WAIS_Calc=GISP2WAIS_Interp_CH4;
        GISP2Calc_Age=GISP2Interp_Age;
    
        %Calculate estimates of the IPD for this iteration and degree 
        %of smoothing:
        GISP2IPGAbs(:,i,j)=GISP2_Smooth_CH4-GISP2WAIS_Smooth_CH4;
        GISP2IPG(:,i,j)=((GISP2_Smooth_CH4./GISP2WAIS_Smooth_CH4)-1)*100;
%--------------------------------------------------- 
        %Calculate interval-average variables:
        GISP2st2_NH=mean(GISP2_Calc(find((GISP2Calc_Age<26956)&...
            (GISP2Calc_Age>24388))));
        GISP2st2_NH_Std=std(GISP2_Calc(find((GISP2Calc_Age<26956)&...
            (GISP2Calc_Age>24388))));
        GISP2st2_SH=mean(GISP2WAIS_Calc(find((GISP2Calc_Age<26956)&...
            (GISP2Calc_Age>24388))));
        GISP2st2_SH_Std=std(GISP2WAIS_Calc(find((GISP2Calc_Age<26956)&...
            (GISP2Calc_Age>24388))));

        GISP2st2_IPG_Abs(i,j)=GISP2st2_NH-GISP2st2_SH;
        GISP2st2_IPG_Abs_Std(i,j)=sqrt(GISP2st2_NH_Std^2+GISP2st2_SH_Std^2);
        GISP2st2_rIPG(i,j)=((GISP2st2_NH/GISP2st2_SH)-1)*100;
%--------------------------------------------------- 
        GISP2DO2_NH=mean(GISP2_Calc(find((GISP2Calc_Age<23097)&...
            (GISP2Calc_Age>22113))));
        GISP2DO2_NH_Std=std(GISP2_Calc(find((GISP2Calc_Age<23097)&...
            (GISP2Calc_Age>22113))));
        GISP2DO2_SH=mean(GISP2WAIS_Calc(find((GISP2Calc_Age<23097)&...
            (GISP2Calc_Age>22113))));
        GISP2DO2_SH_Std=std(GISP2WAIS_Calc(find((GISP2Calc_Age<23097)&...
            (GISP2Calc_Age>22113))));

        GISP2DO2_IPG_Abs(i,j)=GISP2DO2_NH-GISP2DO2_SH;
        GISP2DO2_IPG_Abs_Std(i,j)=sqrt(GISP2DO2_NH_Std^2+GISP2DO2_SH_Std^2);
        GISP2DO2_rIPG(i,j)=((GISP2DO2_NH/GISP2DO2_SH)-1)*100;
%---------------------------------------------------        
        GISP2LGM_NH=mean(GISP2_Calc(find((GISP2Calc_Age<21873)&...
            (GISP2Calc_Age>21217))));
        GISP2LGM_NH_Std=std(GISP2_Calc(find((GISP2Calc_Age<21873)&...
            (GISP2Calc_Age>21217))));
        GISP2LGM_SH=mean(GISP2WAIS_Calc(find((GISP2Calc_Age<21873)&...
            (GISP2Calc_Age>21217))));
        GISP2LGM_SH_Std=std(GISP2WAIS_Calc(find((GISP2Calc_Age<21873)&...
            (GISP2Calc_Age>21217))));

        GISP2LGM_IPG_Abs(i,j)=GISP2LGM_NH-GISP2LGM_SH;
        GISP2LGM_IPG_Abs_Std(i,j)=sqrt(GISP2LGM_NH_Std^2+GISP2LGM_SH_Std^2);
        GISP2LGM_rIPG(i,j)=((GISP2LGM_NH/GISP2LGM_SH)-1)*100;
%---------------------------------------------------         
        GISP2st1_NH=mean(GISP2_Calc(find((GISP2Calc_Age<20415)&...
            (GISP2Calc_Age>17803))));
        GISP2st1_NH_Std=std(GISP2_Calc(find((GISP2Calc_Age<20415)&...
            (GISP2Calc_Age>17803))));
        GISP2st1_SH=mean(GISP2WAIS_Calc(find((GISP2Calc_Age<20415)&...
            (GISP2Calc_Age>17803))));
        GISP2st1_SH_Std=std(GISP2WAIS_Calc(find((GISP2Calc_Age<20415)&...
            (GISP2Calc_Age>17803))));
             
        GISP2st1_IPG_Abs(i,j)=GISP2st1_NH-GISP2st1_SH;
        GISP2st1_IPG_Abs_Std(i,j)=sqrt(GISP2st1_NH_Std^2+GISP2st1_SH_Std^2);
        GISP2st1_rIPG(i,j)=((GISP2st1_NH/GISP2st1_SH)-1)*100;
%---------------------------------------------------         
        GISP2HS1_NH=mean(GISP2_Calc(find((GISP2Calc_Age<17803)&...
            (GISP2Calc_Age>14705))));
        GISP2HS1_NH_Std=std(GISP2_Calc(find((GISP2Calc_Age<17803)&...
            (GISP2Calc_Age>14705))));
        GISP2HS1_SH=mean(GISP2WAIS_Calc(find((GISP2Calc_Age<17803)&...
            (GISP2Calc_Age>14705))));
        GISP2HS1_SH_Std=std(GISP2WAIS_Calc(find((GISP2Calc_Age<17803)&...
            (GISP2Calc_Age>14705))));
        
        GISP2HS1_IPG_Abs(i,j)=GISP2HS1_NH-GISP2HS1_SH;
        GISP2HS1_IPG_Abs_Std(i,j)=sqrt(GISP2HS1_NH_Std^2+GISP2HS1_SH_Std^2);
        GISP2HS1_rIPG(i,j)=((GISP2HS1_NH/GISP2HS1_SH)-1)*100;
        
    end
    
end
%%{
%Save best GISP2 chronology:
GISP2_results=[GISP2_Depth GISP2BestAll GISP2_CH4];
save 'GISP2_Chronology.txt' GISP2_results -ascii

%Average over nt Monte Carlo iterations:
[nn,mm]=size(GISP2IPG);
for j=5:2:15    %loop through different degrees of smoothing
    for k=1:nn  %loop through time steps
        GISP2UpperBound(k,j)=max(GISP2IPG(k,:,j));
        GISP2LowerBound(k,j)=min(GISP2IPG(k,:,j));
        GISP2IPGMean(k,j)=mean(GISP2IPG(k,:,j));
        GISP2IPGStd(k,j)=std(GISP2IPG(k,:,j));
        GISP2IPGAbsMean(k,j)=nanmean(GISP2IPGAbs(k,:,j));
    end
end
%%

%Clear variables for next calculations:
clear GISP2_Calc WAIS_Calc

%Use best-fit chronology to save concentration and IPD results:

%Loop through different degrees of smoothing:
for j=5:2:15
    
    %Generate gaussian smoothing filter with width = w*dt (dt=50):
    w=j;
    window=gausswin(w);
    window=window/sum(window);
    
    %Interpolate both records onto evenly-spaced chronology:
    GISP2_Interp_CH4=interp1(GISP2BestAll,GISP2_CH4,GISP2Interp_Age);
    GISP2_Smooth_CH4=conv(GISP2_Interp_CH4,window,'same');

    GISP2WAIS_Interp_CH4=interp1(GISP2WAIS_gasage,GISP2WAIS_CH4,GISP2Interp_Age);
    GISP2WAIS_Smooth_CH4=conv(GISP2WAIS_Interp_CH4,window,'same');
    
%     Save smoothed records:
    GISP2_Calc(:,j)=GISP2_Smooth_CH4';
    GISP2WAIS_Calc(:,j)=GISP2WAIS_Smooth_CH4';
    GISP2Calc_Age=GISP2Interp_Age;
end

%Choose amount of smoothing to save final results:
saveindex=9;
GISP2IPG_results=[GISP2Calc_Age GISP2_Calc(:,saveindex) GISP2WAIS_Calc(:,saveindex)...
    GISP2IPGAbsMean(:,saveindex) GISP2IPGMean(:,saveindex) GISP2IPGStd(:,saveindex)];
%save IPG_cont_results_onWAIS.txt IPG_results -ascii
%%
GISP2_Calc = GISP2_Calc(:,9);
GISP2WAIS_Calc = GISP2WAIS_Calc(:,9);
%Calculate interval-averages
%--------------------------------------------
GISP2st2_NH=mean(GISP2_Calc(find((GISP2Interp_Age<26956)&(GISP2Interp_Age>24388))));
GISP2st2_NH_Std=std(GISP2_Calc(find((GISP2Interp_Age<26956)&(GISP2Interp_Age>24388))));
GISP2st2_SH=mean(GISP2WAIS_Calc(find((GISP2Interp_Age<26956)&(GISP2Interp_Age>24388))));
GISP2st2_SH_Std=std(GISP2WAIS_Calc(find((GISP2Interp_Age<26956)&(GISP2Interp_Age>24388))));
        
GISP2st2_Age=mean(GISP2Interp_Age(find((GISP2Interp_Age<26956)&(GISP2Interp_Age>24388))));
GISP2st2_Dur=26956-24388;
GISP2st2_IPG=mean(GISP2IPGMean(find((GISP2Interp_Age<26956)&(GISP2Interp_Age>24388)),:));
GISP2st2_IPG_Std1=std(GISP2IPGMean(find((GISP2Interp_Age<26956)&(GISP2Interp_Age>24388)),:));
GISP2st2_IPG_Std=mean(GISP2IPGStd(find((GISP2Interp_Age<26956)&(GISP2Interp_Age>24388)),:));
%BA2_IPG_Std = BA2_IPG_Std1 + BA2_IPG_Std2;
GISP2st2_rIPG_Avg=mean(GISP2st2_rIPG);
GISP2st2_rIPG_Std=std(GISP2st2_rIPG);
%--------------------------------------------
GISP2DO2_NH=mean(GISP2_Calc(find((GISP2Interp_Age<23097)&(GISP2Interp_Age>22113))));
GISP2DO2_NH_Std=std(GISP2_Calc(find((GISP2Interp_Age<23097)&(GISP2Interp_Age>22113))));
GISP2DO2_SH=mean(GISP2WAIS_Calc(find((GISP2Interp_Age<23097)&(GISP2Interp_Age>22113))));
GISP2DO2_SH_Std=std(GISP2WAIS_Calc(find((GISP2Interp_Age<23097)&(GISP2Interp_Age>22113))));
        
GISP2DO2_Age=mean(GISP2Interp_Age(find((GISP2Interp_Age<23097)&(GISP2Interp_Age>22113))));
GISP2DO2_Dur=23097-22113;
GISP2DO2_IPG=mean(GISP2IPGMean(find((GISP2Interp_Age<23097)&(GISP2Interp_Age>22113)),:));
GISP2DO2_IPG_Std1=std(GISP2IPGMean(find((GISP2Interp_Age<23097)&(GISP2Interp_Age>22113)),:));
GISP2DO2_IPG_Std=mean(GISP2IPGStd(find((GISP2Interp_Age<23097)&(GISP2Interp_Age>22113)),:));
%BA1_IPG_Std = BA1_IPG_Std1 + BA1_IPG_Std2;
GISP2DO2_rIPG_Avg=mean(GISP2DO2_rIPG);
GISP2DO2_rIPG_Std=std(GISP2DO2_rIPG);
%--------------------------------------
GISP2LGM_NH=mean(GISP2_Calc(find((GISP2Interp_Age<21873)&(GISP2Interp_Age>21217))));
GISP2LGM_NH_Std=std(GISP2_Calc(find((GISP2Interp_Age<21873)&(GISP2Interp_Age>21217))));
GISP2LGM_SH=mean(GISP2WAIS_Calc(find((GISP2Interp_Age<21873)&(GISP2Interp_Age>21217))));
GISP2LGM_SH_Std=std(GISP2WAIS_Calc(find((GISP2Interp_Age<21873)&(GISP2Interp_Age>21217))));

GISP2LGM_Age=mean(GISP2Interp_Age(find((GISP2Interp_Age<21873)&(GISP2Interp_Age>21217))));
GISP2LGM_Dur=21873-21217;
GISP2LGM_IPG=mean(GISP2IPGMean(find((GISP2Interp_Age<21873)&(GISP2Interp_Age>21217)),:));
GISP2LGM_IPG_Std1=std(GISP2IPGMean(find((GISP2Interp_Age<21873)&(GISP2Interp_Age>21217)),:));
GISP2LGM_IPG_Std=mean(GISP2IPGStd(find((GISP2Interp_Age<21873)&(GISP2Interp_Age>21217)),:));
%LGM_IPG_Std = LGM_IPG_Std1 + LGM_IPG_Std2;                                  %Need to include mean uncertainty in IPD
GISP2LGM_rIPG_Avg=mean(GISP2LGM_rIPG);
GISP2LGM_rIPG_Std=std(GISP2LGM_rIPG);
%--------------------------------------------
GISP2st1_NH=mean(GISP2_Calc(find((GISP2Interp_Age<20415)&(GISP2Interp_Age>17803))));
GISP2st1_NH_Std=std(GISP2_Calc(find((GISP2Interp_Age<20415)&(GISP2Interp_Age>17803))));
GISP2st1_SH=mean(GISP2WAIS_Calc(find((GISP2Interp_Age<20415)&(GISP2Interp_Age>17803))));
GISP2st1_SH_Std=std(GISP2WAIS_Calc(find((GISP2Interp_Age<20415)&(GISP2Interp_Age>17803))));
        
GISP2st1_Age=mean(GISP2Interp_Age(find((GISP2Interp_Age<20415)&(GISP2Interp_Age>17803))));
GISP2st1_Dur=20415-17803;
GISP2st1_IPG=mean(GISP2IPGMean(find((GISP2Interp_Age<20415)&(GISP2Interp_Age>17803)),:));
GISP2st1_IPG_Std1=std(GISP2IPGMean(find((GISP2Interp_Age<20415)&(GISP2Interp_Age>17803)),:));
GISP2st1_IPG_Std=mean(GISP2IPGStd(find((GISP2Interp_Age<20415)&(GISP2Interp_Age>17803)),:));
%OD2_IPG_Std = OD2_IPG_Std1 + OD2_IPG_Std2;
GISP2st1_rIPG_Avg=mean(GISP2st1_rIPG);
GISP2st1_rIPG_Std=std(GISP2st1_rIPG);
%--------------------------------------------
GISP2HS1_NH=mean(GISP2_Calc(find((GISP2Interp_Age<17803)&(GISP2Interp_Age>14705))));
GISP2HS1_NH_Std=std(GISP2_Calc(find((GISP2Interp_Age<17803)&(GISP2Interp_Age>14705))));
GISP2HS1_SH=mean(GISP2WAIS_Calc(find((GISP2Interp_Age<17803)&(GISP2Interp_Age>14705))));
GISP2HS1_SH_Std=std(GISP2WAIS_Calc(find((GISP2Interp_Age<17803)&(GISP2Interp_Age>14705))));
        
GISP2HS1_Age=mean(GISP2Interp_Age(find((GISP2Interp_Age<17803)&(GISP2Interp_Age>14705))));
GISP2HS1_Dur=17803-14705;
GISP2HS1_IPG=mean(GISP2IPGMean(find((GISP2Interp_Age<17803)&(GISP2Interp_Age>14705)),:));
GISP2HS1_IPG_Std1=std(GISP2IPGMean(find((GISP2Interp_Age<17803)&(GISP2Interp_Age>14705)),:));
GISP2HS1_IPG_Std=mean(GISP2IPGStd(find((GISP2Interp_Age<17803)&(GISP2Interp_Age>14705)),:));
%OD1_IPG_Std = OD1_IPG_Std1 + OD1_IPG_Std2;
GISP2HS1_rIPG_Avg=mean(GISP2HS1_rIPG);
GISP2HS1_rIPG_Std=std(GISP2HS1_rIPG);

% %--------------------------------------------

%Save results:
GISP2All_IPG=   [GISP2HS1_Age GISP2HS1_Dur GISP2HS1_IPG GISP2HS1_IPG_Std;
            GISP2st1_Age GISP2st1_Dur GISP2st1_IPG GISP2st1_IPG_Std;
            GISP2LGM_Age GISP2LGM_Dur GISP2LGM_IPG GISP2LGM_IPG_Std;
            GISP2DO2_Age GISP2DO2_Dur GISP2DO2_IPG GISP2DO2_IPG_Std;
            GISP2st2_Age GISP2st2_Dur GISP2st2_IPG GISP2st2_IPG_Std];

%% compile all data in cell array and save array
gradient(un,1) = {[NEEMBestAll,NEEM_CH4,NEEM_Depth]};
gradient(un,2) = {[NEEMWAIS_Age,NEEMWAIS_CH4]};
gradient(un,3) = {[NEEMTieAge]};
gradient(un,4) = {[NEEMInterp_Age,NEEMIPGMean,NEEMIPGStd,NEEMIPGAbsMean(:,saveindex)]};
gradient(un,5) = {[NEEMAll_IPG]};
gradient(un,6) = {[GISP2BestAll,GISP2_CH4,GISP2_Depth]};
gradient(un,7) = {[GISP2WAIS_gasage,GISP2WAIS_CH4]};
gradient(un,8) = {[GISP2TieAge]};
gradient(un,9) = {[GISP2Interp_Age,GISP2IPGMean,GISP2IPGStd,GISP2IPGAbsMean(:,saveindex)]};
gradient(un,10) = {[GISP2All_IPG]};

%Save interval averages 
if un == 8
NEEM_IA = [NEEMHS1_NH, NEEMHS1_SH; 
           NEEMBA2_NH, NEEMBA2_SH;
           NEEMBA1_NH, NEEMBA1_SH;
           NEEMYD_NH, NEEMYD_SH];
GISP2_IA = [GISP2st2_NH, GISP2st2_SH;
            GISP2DO2_NH, GISP2DO2_SH;
            GISP2LGM_NH, GISP2LGM_SH;
            GISP2st1_NH, GISP2st1_SH;
            GISP2HS1_NH, GISP2HS1_SH];
end 

%% compile 
end 
%Remove 0 rows and save all data 
gradient = [gradient(2,:);gradient(8,:)];
save('gradient');
save('NEEM_IA');
save('GISP2_IA');
