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
TPexchange = 0.22; %Rasmussen et al., 1984 (should probably update this)
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

%Load data: For GISP2 data, use column 8 for CH4xs corrected data and column
%2 for uncorrected data.
GISP2_data=csvread('GISP2final_gradient_nov22.csv',1,0);
GISP2_Depth=GISP2_data(:,1);
GISP2_Age=GISP2_data(:,5);
GISP2_CH4=(GISP2_data(:,8));
GISP2_Ca2 = GISP2_data(:,7);
%Average duplicates 
     [C2,ia2,idx2] = unique(GISP2_Depth,'stable');
     GISP2_Depth = accumarray(idx2,GISP2_Depth,[],@mean); 
     GISP2_CH4 = accumarray(idx2,GISP2_CH4,[],@mean); 
     GISP2_Ca2 = accumarray(idx2,GISP2_Ca2,[],@mean); 

WAIS_data=csvread('WAISfinal_GISP2_may22.csv',1,0);
%Average duplicates
for i = 1:(length(WAIS_data)-1)
    if WAIS_data((i+1),1) == WAIS_data(i,1)
        for j = 1:length(WAIS_data(1,:))
            WAIS_data((i+1),j) = mean(WAIS_data(i:(i+1),j));
        end
        WAIS_data(i,:) = [0];
    else 
        WAIS_data(i,:) = WAIS_data(i,:);
    end 
end
WAIS_data( ~any(WAIS_data,2), : ) = [];
WAIS_Depth=WAIS_data(:,1);
%WAIS_Age=WAIS_data(:,4);
WAIS_CH4=WAIS_data(:,2);
WAIS_Age = interp1(WAISdepth_Chronology,WAISgasage_Chronology,WAIS_Depth);


nn=length(GISP2_Age);
mm=length(WAIS_Age);
%Initialize solution array:
F=zeros(nn,4,nt);
%Generate gaussian smoothing filter with width = w*dt (dt=50):
w=9; %use this for whole record
%w=12; %use this for abrupt transitions
window=gausswin(w);
window=window/sum(window);
%Create evenly-spaced age scale for interpolation:
dt=50;
Interp_Age=[ceil(min(WAIS_Age)):dt:floor(max(WAIS_Age))];
nnn=length(Interp_Age);
%Define initial tie points:
            %GISP2 depth %WAIS depth
StartTies= [2065.924465, 2747.303781;
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
%Define Southern Source strength
FS=[Interp_Age',ones(254,1)*15];
FS2 = [Interp_Age',ones(254,1)*1]; 
%Define variables:
[n,m]=size(StartTies);
TieDepthGISP2=StartTies(:,1);
TieDepthWAIS=StartTies(:,2);
TieAge = interp1 (WAISdepth_Chronology,WAISgasage_Chronology,TieDepthWAIS);

%Monte Carlo analysis:
%Uncertainty in CH4xs correction (applied to each Monte Carlo Simulation,
%rather than each sample).
CH4xs_u = normrnd(0,0.0012,[1000,1]);
for i=1:nt
%Initialize this for variable sink experiments:
%modeltime=Interp_Age(1);
 
%calculate GISP2 stdev
GISP2_stdev = 3;
GISP2_CH4_Error = GISP2_CH4+normrnd(0,GISP2_stdev,[length(GISP2_CH4),1])+(GISP2_Ca2*8.314*273*CH4xs_u(i)/(40.078*.092*101.325));
WAIS_CH4_Error = WAIS_CH4+normrnd(0,2,[length(WAIS_CH4),1]);
%Create evenly-spaced chronological anchors:
nties=50; %length(StartTies);
TieAgeTemp=linspace(TieAge(1),TieAge(end),nties)';
TieDepthTemp=interp1(TieAge,TieDepthGISP2,TieAgeTemp);
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
%Interpolate GISP2 data with error onto new chronology:
GISP2_Age=interp1(TieDepthTotal,TieAgeTotal,GISP2_Depth);
GISP2_Interp_CH4=interp1(GISP2_Age,GISP2_CH4_Error,Interp_Age);
%Interpolate WAIS data onto evenly-spaced chronology:
WAIS_Interp_CH4=interp1(WAIS_Age,WAIS_CH4_Error,Interp_Age);
%Calculate correlation between records:
Rs=corrcoef(GISP2_Interp_CH4,WAIS_Interp_CH4);
R(i)=Rs(2,1);
%Calculate gas age annual layer thickness as a check:
ALT=(GISP2_Depth(2:end)-GISP2_Depth(1:end-1))./(GISP2_Age(2:end)...
-GISP2_Age(1:end-1));
%Save GISP2 age scale that results in best fit:
if R(i)>=max(R(1:i-1))
BestAll=GISP2_Age;
RBestAll=R(i);
TieBestAll=[TieDepthTotal TieAgeTotal];
ALTAll=ALT;
end
%Smooth both sets of data with filter:
GISP2_Smooth_CH4=conv(GISP2_Interp_CH4,window,'same');
WAIS_Smooth_CH4=conv(WAIS_Interp_CH4,window,'same');
%Transfer variable names for next calculations:
GISP2_Calc=GISP2_Smooth_CH4;
WAIS_Calc=WAIS_Smooth_CH4;
Calc_Age=Interp_Age;
%Calculate rIPD for each iteration:
IPG(:,i)=((GISP2_Calc./WAIS_Calc)-1)*100;
%Calculate source fluxes:
for j=1:nnn

%%
%Reset solution vector:
C=zeros(4,1);
%Calculate burdens in N and S boxes (%need to adjust NH concentration by 7%
%of gradient (Brook et al, 2000). S box is homogeneous.:
C(1,1)=((GISP2_Calc(j)-0.07*(GISP2_Calc(j) - WAIS_Calc(j)))/(10^9))*...
Mtrop/4*MCH4/(10^12);
C(4,1)=(WAIS_Calc(j)/(10^9))*Mtrop/4*MCH4/(10^12);
FSB = FS2(j,2);

% Solve the box model
C(3,1) = (FSB - (lambdaS + nP)*C(4,1)) / -nP;
C(2,1) = (-(lambdaT + nT + nP)*C(3,1) + nP*C(4,1) - nP*C(1,1) - nT*C(3,1)) / (-nT - lambdaT - nP - nT);
Fn(j,i) = (lambdaN + nP)*C(1,1) -nP*C(2,1);
Ft(j,i) = -nP*C(1,1) + (lambdaT + nP + nT)*C(2,1) - nT*C(3,1);
Ft2(j,i) = -nT*C(2,1) + (lambdaT + nP + nT)*C(3,1) - nP*C(4,1);

end
end
%%
%Rearrange source estimates into NH, tropical, and SH arrays:
for i=1:nt
FN(:,i)=Fn(:,i); FT(:,i)=Ft(:,i)*2; FS(:,i)=ones(254,1)*12;
end
%Calculate max, min, mean, and std of box sources and IPD:
for i=1:nnn
FNMax(i,1)=max(FN(i,:)); FNMin(i,1)=min(FN(i,:));
FNMean(i,1)=mean(FN(i,:)); FNStd(i,1)=std(FN(i,:));
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
GISP2_Interp_CH4=interp1(BestAll,GISP2_CH4,Interp_Age);
GISP2_Smooth_CH4=conv(GISP2_Interp_CH4,window,'same');
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

GISP24box = cell(8,1);
GISP24box(1,1) = {[Interp_Age',IPGMean,IPGStd]};
GISP24box(2,1) = {[BestAll,GISP2_CH4]};
GISP24box(3,1) = {[WAIS_Age,WAIS_CH4]};
GISP24box(4,1) = {[Interp_Age',FNMean,FNMin,FNMax,FNStd]};
GISP24box(5,1) = {[Interp_Age',FTMean,FTMin,FTMax,FTStd]};
GISP24box(6,1) = {[Interp_Age',FS2(:,2)]};
GISP24box(7,1) = {[Calc_Age',GISP2_Smooth_CH4']};
GISP24box(8,1) = {[Calc_Age',WAIS_Smooth_CH4']};

save('GISP24box');
