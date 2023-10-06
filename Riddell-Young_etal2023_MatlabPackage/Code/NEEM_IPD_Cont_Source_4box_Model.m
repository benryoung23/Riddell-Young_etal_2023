%=====================================================================
% Continuous 4-Box Source Model
%=====================================================================
%This code generates a continuous estimate of the distribution of
%methane sources across the last glacial termination using the output of
%GISP2_and_NEEMrIPDa.m. The output of this code is published in Riddell-Young et al., 2023.
% "Tropical Sources mainly controlled last glacial maximum and deglacial
% methane variability." Nature Geoscience.
% 
% The script uses a four-box model and monte-carlo style error propagation
% to calculate the methane source distribution across the last glacial
% maximum and deglaciation, as well as uncertainties. 
%---------------------------------------------------------------------
%%
clear all
%clf
format long eng
%Number of Monte Carlo iterations (was set to 10000):
nt=1000;
%Define parameters of 4-box model:
%Define exchange rates
TPexchange = 0.22; %Rasmussen et al., 1984 was 0.21 and 0.59. These are updated vallues based on new data.
TTexchange = 0.45; 
nP = TPexchange^-1;
nT = TTexchange^-1;
%Define lifetimes
Nlifetime = 15.6;  %Baumgartner 2012
Tlifetime = 6.8;
Slifetime = 22.4;
lambdaN = 1/Nlifetime;
lambdaT = 1/Tlifetime;
lambdaS = 1/Slifetime;
%Calculate removal matrix:
omega = [(lambdaN + nP) -nP 0 0;
    -nP (lambdaT + nP + nT) -nT 0;
    0 -nT (lambdaT + nP + nT) -nP;
    0 0 -nP (lambdaS + nP)];
%Model constants:
Mtrop=178E18; %mass of troposphere is moles of dry air from
%Trenberth 2005
MCH4=16; %g per mole of CH4

%Load WAIS chronology (Buizert, 2015).
WAIS_Chronology = csvread('WAIS_Chronology_Buizert2015.csv',1,0);
WAISdepth_Chronology = WAIS_Chronology(:,1);
WAISiceage_Chronology = WAIS_Chronology(:,2);
WAISgasage_Chronology = WAIS_Chronology(:,3);

%Load data: For NEEM data, use column 8 for CH4xs corrected data and column
%2 for uncorrected data.
NEEM_data=csvread('NEEMfinal_nov22.csv',1,0);
%Average duplicates
     [Depth,ia2,idx2] = unique(NEEM_data(:,1),'stable');
     val2 = accumarray(idx2,NEEM_data(:,8),[],@mean);           %un changes between 2 and 8 for uncorrected and corrected GISP2 data, respectively.
     val3 = accumarray(idx2,NEEM_data(:,7),[],@mean);
     val4 = accumarray(idx2,NEEM_data(:,5),[],@mean);
     NEEM_Data2 = [Depth val2 val3 val4];   

%Define variables:
NEEM_Depth = NEEM_Data2(:,1);
NEEM_Age = NEEM_Data2(:,4)-50;  %The NEEM gas age scale is reported on yb2k. WD and GISP2 are by1950, so 50 years was subtracted from the age to match GISP2 and WD.
NEEM_CH4 = NEEM_Data2(:,2); 
NEEM_Ca2 = NEEM_Data2(:,3);

WAIS_data=csvread('WAISfinal_NEEM_may22.csv',1,0);
WAIS_Depth=WAIS_data(:,1);
%WAIS_Age=WAIS_data(:,4);
WAIS_CH4=WAIS_data(:,2);
WAIS_Age = interp1(WAISdepth_Chronology,WAISgasage_Chronology,WAIS_Depth);
%Average duplicates 
     [C2,ia2,idx2] = unique(WAIS_Depth,'stable');
     WAIS_Depth = accumarray(idx2,WAIS_Depth,[],@mean); 
     WAIS_CH4 = accumarray(idx2,WAIS_CH4,[],@mean); 
     WAIS_Age = accumarray(idx2,WAIS_Age,[],@mean); 

nn=length(NEEM_Age);
mm=length(WAIS_Age);
%Initialize solution array:
F=zeros(nn,4,nt);
%Generate gaussian smoothing filter with width = w*dt (dt=50):
w=9; %use this for whole record
window=gausswin(w);
window=window/sum(window);
%Create evenly-spaced age scale for interpolation:
dt=50;
Interp_Age=[ceil(min(WAIS_Age)):dt:floor(max(WAIS_Age))];
nnn=length(Interp_Age);
%Define initial tie points:
            %NEEM depth %WAIS depth
StartTies= [1558.00 2480.375;          %From Julia Rosen, but changed to a more realistic value
            1540.88 2441.831;          %From Julia Rosen
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
%Define Southern Source strength
FS=[Interp_Age',ones(152,1)*15];
FS2 = [Interp_Age',ones(152,1)*1]; 
%Define variables:
[n,m]=size(StartTies);
TieDepthNEEM=StartTies(:,1);
TieDepthWAIS=StartTies(:,2);
TieAge = interp1 (WAISdepth_Chronology,WAISgasage_Chronology,TieDepthWAIS);

%Monte Carlo analysis:
%Uncertainty in CH4xs correction (applied to each Monte Carlo Simulation,
%rather than each sample).
CH4xs_u = normrnd(0,0.0012,[1000,1]);
for i=1:nt
%Initialize this for variable sink experiments:
%modeltime=Interp_Age(1);
                                                            
%calculate NEEM stdev
NEEM_stdev = 4;
NEEM_CH4_Error = NEEM_CH4+normrnd(0,NEEM_stdev,[length(NEEM_CH4),1])+(NEEM_Ca2*8.314*273*CH4xs_u(i)/(40.078*.092*101.325));
WAIS_CH4_Error = WAIS_CH4+normrnd(0,2,[length(WAIS_CH4),1]);
%Create evenly-spaced chronological anchors:
nties=50; %length(StartTies);
TieAgeTemp=linspace(TieAge(1),TieAge(end),nties)';
TieDepthTemp=interp1(TieAge,TieDepthNEEM,TieAgeTemp);
%Define tie point age uncertainty:
TieFix=20;
TieError=[TieFix*randn(nties-2,1)];
%Add uncertainty to all points except first and last:
TieAgeError=[TieAgeTemp(1); TieAgeTemp(2:end-1)+TieError;...
TieAgeTemp(end)];
%Make sure tie points are still in decreasing order after adding
%age uncertainty:
Test=TieAgeError(2:end)-TieAgeError(1:end-1);
%If not, redo tie points until they all decrease:
while max(Test)>0
TieError=TieFix*randn(nties-2,1);
TieAgeError=[TieAgeTemp(1); TieAgeTemp(2:end-1)+TieError;...
TieAgeTemp(end)];
%Check again:
Test=TieAgeError(2:end)-TieAgeError(1:end-1);
end
%Define tie points for this iteration:
TieDepthTotal=TieDepthTemp;
TieAgeTotal=TieAgeError;
%%
%Interpolate NEEM data with error onto new chronology:
NEEM_Age=interp1(TieDepthTotal,TieAgeTotal,NEEM_Depth);
NEEM_Interp_CH4=interp1(NEEM_Age,NEEM_CH4_Error,Interp_Age);
%Interpolate WAIS data onto evenly-spaced chronology:
WAIS_Interp_CH4=interp1(WAIS_Age,WAIS_CH4_Error,Interp_Age);
%Calculate correlation between records:
Rs=corrcoef(NEEM_Interp_CH4,WAIS_Interp_CH4);
R(i)=Rs(2,1);
%Calculate gas age annual layer thickness as a check:
ALT=(NEEM_Depth(2:end)-NEEM_Depth(1:end-1))./(NEEM_Age(2:end)...
-NEEM_Age(1:end-1));
%Save NEEM age scale that results in best fit:
if R(i)>=max(R(1:i-1))
BestAll=NEEM_Age;
RBestAll=R(i);
TieBestAll=[TieDepthTotal TieAgeTotal];
ALTAll=ALT;
end
%Smooth both sets of data with filter:
NEEM_Smooth_CH4=conv(NEEM_Interp_CH4,window,'same');
WAIS_Smooth_CH4=conv(WAIS_Interp_CH4,window,'same');
%Transfer variable names for next calculations:
NEEM_Calc=NEEM_Smooth_CH4;
WAIS_Calc=WAIS_Smooth_CH4;
Calc_Age=Interp_Age;
%Calculate rIPD for each iteration:
IPG(:,i)=((NEEM_Calc./WAIS_Calc)-1)*100;
%Calculate source fluxes:
for j=1:nnn

%%
%Reset solution vector:
C=zeros(4,1);
%Calculate burdens in N and S boxes (%need to adjust NH concentration by 7%
%of gradient (Brook et al, 2000). S box is homogeneous.:
C(1,1)=((NEEM_Calc(j)-0.07*(NEEM_Calc(j) - WAIS_Calc(j)))/(10^9))*...
Mtrop/4*MCH4/(10^12);
C(4,1)=(WAIS_Calc(j)/(10^9))*Mtrop/4*MCH4/(10^12);
FSB = FS2(j,2);

% Solve the box model
C(3,1) = (FSB - (lambdaS + nP)*C(4,1)) / -nP;
CTS(j,i) = C(3,1);
C(2,1) = (-(lambdaT + nT + nP)*C(3,1) + nP*C(4,1) - nP*C(1,1) - nT*C(3,1)) / (-nT - lambdaT - nP - nT);
CTN(j,i) = C(2,1);
Fn(j,i) = (lambdaN + nP)*C(1,1) -nP*C(2,1);
Ft(j,i) = -nP*C(1,1) + (lambdaT + nP + nT)*C(2,1) - nT*C(3,1);
Ft2(j,i) = -nT*C(2,1) + (lambdaT + nP + nT)*C(3,1) - nP*C(4,1);

end
end
%%
%Rearrange source estimates into NH, tropical, and SH arrays:
for i=1:nt
FN(:,i)=Fn(:,i); FT(:,i)=Ft(:,i)*2; FS(:,i)=ones(152,1)*15;
end
%Calculate max, min, mean, and std of box sources and IPD:
for i=1:nnn
FNMax(i,1)=max(FN(i,:)); FNMin(i,1)=min(FN(i,:));
FNMean(i,1)=mean(FN(i,:)); FNStd(i,1)=std(FN(i,:));
CTNMean(1,1)=mean(CTN(i,:));
CTSMean(1,1)=mean(CTS(i,:));
FTMax(i,1)=max(FT(i,:)); FTMin(i,1)=min(FT(i,:));
FTMean(i,1)=mean(FT(i,:)); FTStd(i,1)=std(FT(i,:));
FSMax(i,1)=max(FS(i,:)); FSMin(i,1)=min(FS(i,:));
FSMean(i,1)=mean(FS(i,:)); FSStd(i,1)=std(FS(i,:));
IPGMax(i,1)=max(IPG(i,:)); IPGMin(i,1)=min(IPG(i,:));
IPGMean(i,1)=mean(IPG(i,:)); IPGStd(i,1)=std(IPG(i,:));
end
%Calculate interval-average NH source estimates and standard
%deviations:
Avg_NSource(1,1)=mean(FNMean(find((Interp_Age<18500)&...
(Interp_Age>17800))));
Avg_NSource(1,2)=std(FNMean(find((Interp_Age<18500)&...
(Interp_Age>17800)))) + mean(FNMax(find((Interp_Age<18500)&...
(Interp_Age>17800)))) - mean(FNMean(find((Interp_Age<18500)&...
(Interp_Age>17800))));
Avg_NSource(2,1)=mean(FNMean(find((Interp_Age<17800)&...
(Interp_Age>16400))));
Avg_NSource(2,2)=std(FNMean(find((Interp_Age<17800)&...
(Interp_Age>16400)))) + mean(FNMax(find((Interp_Age<17800)&...
(Interp_Age>16400)))) - mean(FNMean(find((Interp_Age<17800)&...
(Interp_Age>16400))));
Avg_NSource(3,1)=mean(FNMean(find((Interp_Age<16000)&...
(Interp_Age>14900))));
Avg_NSource(3,2)=std(FNMean(find((Interp_Age<16000)&...
(Interp_Age>14900)))) + mean(FNMax(find((Interp_Age<16000)&...
(Interp_Age>14900)))) - mean(FNMean(find((Interp_Age<16000)&...
(Interp_Age>14900))));
Avg_NSource(4,1)=mean(FNMean(find((Interp_Age<14500)&...
(Interp_Age>14150))));
Avg_NSource(4,2)=std(FNMean(find((Interp_Age<14500)&...
(Interp_Age>14150)))) + mean(FNMax(find((Interp_Age<14500)&...
(Interp_Age>14150)))) - mean(FNMean(find((Interp_Age<14500)&...
(Interp_Age>14150))));
Avg_NSource(5,1)=mean(FNMean(find((Interp_Age<13600)&...
(Interp_Age>13400))));
Avg_NSource(5,2)=std(FNMean(find((Interp_Age<13600)&...
(Interp_Age>13400)))) + mean(FNMax(find((Interp_Age<13600)&...
(Interp_Age>13400)))) - mean(FNMean(find((Interp_Age<13600)&...
(Interp_Age>13400))));
Avg_NSource(6,1)=mean(FNMean(find((Interp_Age<12500)&...
(Interp_Age>11900))));
Avg_NSource(6,2)=std(FNMean(find((Interp_Age<12500)&...
(Interp_Age>11900)))) + mean(FNMax(find((Interp_Age<12500)&...
(Interp_Age>11900)))) - mean(FNMean(find((Interp_Age<12500)&...
(Interp_Age>11900))));
%Calculate interval-average tropical source estimates and standard
%deviations:
Avg_TSource(1,1)=mean(FTMean(find((Interp_Age<18500)&...
(Interp_Age>17800))));
Avg_TSource(1,2)=std(FTMean(find((Interp_Age<18500)&...
(Interp_Age>17800)))) + mean(FTMax(find((Interp_Age<18500)&...
(Interp_Age>17800)))) - mean(FTMean(find((Interp_Age<18500)&...
(Interp_Age>17800))));
Avg_TSource(2,1)=mean(FTMean(find((Interp_Age<17800)&...
(Interp_Age>16400))));
Avg_TSource(2,2)=std(FTMean(find((Interp_Age<17800)&...
(Interp_Age>16400)))) + mean(FTMax(find((Interp_Age<17800)&...
(Interp_Age>16400)))) - mean(FTMean(find((Interp_Age<17800)&...
(Interp_Age>16400))));
Avg_TSource(3,1)=mean(FTMean(find((Interp_Age<16000)&...
(Interp_Age>14900))));
Avg_TSource(3,2)=std(FTMean(find((Interp_Age<16000)&...
(Interp_Age>14900)))) + mean(FTMax(find((Interp_Age<16000)&...
(Interp_Age>14900)))) - mean(FTMean(find((Interp_Age<16000)&...
(Interp_Age>14900))));
Avg_TSource(4,1)=mean(FTMean(find((Interp_Age<14500)&...
(Interp_Age>14150))));
Avg_TSource(4,2)=std(FTMean(find((Interp_Age<14500)&...
(Interp_Age>14150)))) + mean(FTMax(find((Interp_Age<14500)&...
(Interp_Age>14150)))) - mean(FTMean(find((Interp_Age<14500)&...
(Interp_Age>14150))));
Avg_TSource(5,1)=mean(FTMean(find((Interp_Age<13600)&...
(Interp_Age>13400))));
Avg_TSource(5,2)=std(FTMean(find((Interp_Age<13600)&...
(Interp_Age>13400)))) + mean(FTMax(find((Interp_Age<13600)&...
(Interp_Age>13400)))) - mean(FTMean(find((Interp_Age<13600)&...
(Interp_Age>13400))));
Avg_TSource(6,1)=mean(FTMean(find((Interp_Age<12500)&...
(Interp_Age>11900))));
Avg_TSource(6,2)=std(FTMean(find((Interp_Age<12500)&...
(Interp_Age>11900)))) + mean(FTMax(find((Interp_Age<12500)&...
(Interp_Age>11900)))) - mean(FTMean(find((Interp_Age<12500)&...
(Interp_Age>11900))));
%Calculate interval-average tropical source estimates and standard
%deviations:
Avg_SSource(1,1)=mean(FSMean(find((Interp_Age<18500)&...
(Interp_Age>17800))));
Avg_SSource(1,2)=std(FSMean(find((Interp_Age<18500)&...
(Interp_Age>17800)))) + mean(FSMax(find((Interp_Age<18500)&...
(Interp_Age>17800)))) - mean(FSMean(find((Interp_Age<18500)&...
(Interp_Age>17800))));
Avg_SSource(2,1)=mean(FSMean(find((Interp_Age<17800)&...
(Interp_Age>16400))));
Avg_SSource(2,2)=std(FSMean(find((Interp_Age<17800)&...
(Interp_Age>16400)))) + mean(FSMax(find((Interp_Age<17800)&...
(Interp_Age>16400)))) - mean(FSMean(find((Interp_Age<17800)&...
(Interp_Age>16400))));
Avg_SSource(3,1)=mean(FSMean(find((Interp_Age<16000)&...
(Interp_Age>14900))));
Avg_SSource(3,2)=std(FSMean(find((Interp_Age<16000)&...
(Interp_Age>14900)))) + mean(FSMax(find((Interp_Age<16000)&...
(Interp_Age>14900)))) - mean(FSMean(find((Interp_Age<16000)&...
(Interp_Age>14900))));
Avg_SSource(4,1)=mean(FSMean(find((Interp_Age<14500)&...
(Interp_Age>14150))));
Avg_SSource(4,2)=std(FSMean(find((Interp_Age<14500)&...
(Interp_Age>14150)))) + mean(FSMax(find((Interp_Age<14500)&...
(Interp_Age>14150)))) - mean(FSMean(find((Interp_Age<14500)&...
(Interp_Age>14150))));
Avg_SSource(5,1)=mean(FSMean(find((Interp_Age<13600)&...
(Interp_Age>13400))));
Avg_SSource(5,2)=std(FSMean(find((Interp_Age<13600)&...
(Interp_Age>13400)))) + mean(FSMax(find((Interp_Age<13600)&...
(Interp_Age>13400)))) - mean(FSMean(find((Interp_Age<13600)&...
(Interp_Age>13400))));
Avg_SSource(6,1)=mean(FSMean(find((Interp_Age<12500)&...
(Interp_Age>11900))));
Avg_SSource(6,2)=std(FSMean(find((Interp_Age<12500)&...
(Interp_Age>11900)))) + mean(FSMax(find((Interp_Age<12500)&...
(Interp_Age>11900)))) - mean(FSMean(find((Interp_Age<12500)&...
(Interp_Age>11900))));

%Use best-fit chronology to plot concentrations:
NEEM_Interp_CH4=interp1(BestAll,NEEM_CH4,Interp_Age);
NEEM_Smooth_CH4=conv(NEEM_Interp_CH4,window,'same');
WAIS_Interp_CH4=interp1(WAIS_Age,WAIS_CH4,Interp_Age);
WAIS_Smooth_CH4=conv(WAIS_Interp_CH4,window,'same');
%Choose period to plot:
period=1; %1 = all, 2 = OD-BA, 3 = BA-YD
if period==1
starttime=11500; endtime=19000; CH4lower=300; CH4upper=800;
elseif period==2
starttime=14400; endtime=15000; CH4lower=450; CH4upper=720;
else
starttime=12400; endtime=13200; CH4lower=450; CH4upper=720;
end

%Load Baumgartner 3 box data
%              AgeM  Dur  N Nstd T Tstd S Sstd
Baumgartner = [19019 2611 39 2 58 4 12 0
               16254 3099 50 7 65 14 12 0
               14044 617 65 4 104 11 15 0
               13339 792 66 5 112 12 15 0
               12078 742 48 4 79 8 15 0];

%% Compile and export output

NEEM4box = cell(8,1);
NEEM4box(1,1) = {[Interp_Age',IPGMean,IPGStd]};
NEEM4box(2,1) = {[BestAll,NEEM_CH4]};
NEEM4box(3,1) = {[WAIS_Age,WAIS_CH4]};
NEEM4box(4,1) = {[Interp_Age',FNMean,FNMin,FNMax,FNStd]};
NEEM4box(5,1) = {[Interp_Age',FTMean,FTMin,FTMax,FTStd]};
NEEM4box(6,1) = {[Interp_Age',FS2(:,2)]};
NEEM4box(7,1) = {[Calc_Age',NEEM_Smooth_CH4']};
NEEM4box(8,1) = {[Calc_Age',WAIS_Smooth_CH4']};

save('NEEM4box');
