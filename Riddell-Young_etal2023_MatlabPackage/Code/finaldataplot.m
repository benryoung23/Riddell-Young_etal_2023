%This code plots all of the raw data for each core separated by measurement
%campaign. Colors correspond to the cores used, shapes correspond to the
%campaign (see legend). The data plotted here is corrected for
%solubility, blanks, gravitational fractionation, and campaign offsets.

close all
clear all

%% Load data

%NEEM first
NEEMdata = readtable('NEEM_Alldata_Nov22.xlsx');
NEEMdate = table2array(NEEMdata(:,1));
NEEMdepth = table2array(NEEMdata(:,5));
NEEMgasage = table2array(NEEMdata(:,7))-50; %50 years subtracted to report on the yb1950 rather than yb2k scale.
NEEMCH4u = table2array(NEEMdata(:,16));
NEEMCH4 = table2array(NEEMdata(:,20));

%Then GISP
GISP2data = readtable('GISP2_Alldata_Nov22.xlsx');
GISP2date = table2array(GISP2data(:,1));
GISP2depth = table2array(GISP2data(:,5));
GISP2gasage = table2array(GISP2data(:,7));
GISP2CH4u = table2array(GISP2data(:,16));
GISP2CH4 = table2array(GISP2data(:,20));

%Then WAIS
WAISolddata = readtable('WAIS_2010-13_Alldata.xlsx');
WAISolddate = table2array(WAISolddata(:,1));
WAISolddepth = table2array(WAISolddata(:,5));
WAISoldgasage = table2array(WAISolddata(:,7));
WAISoldCH4 = table2array(WAISolddata(:,19));

WAISnewdata = readtable('WAIS_2018-19_Alldata.xlsx');
WAISnewdate = table2array(WAISnewdata(:,1));
WAISnewdepth = table2array(WAISnewdata(:,5));
WAISnewgasage = table2array(WAISnewdata(:,7));
WAISnewCH4 = table2array(WAISnewdata(:,19));

%Change dates to year
     NEEMdate = year(NEEMdate);
     GISP2date = year(GISP2date); 
     WAISnewdate = year(WAISnewdate); 
     WAISolddate = year(WAISolddate); 

%% Average duplicates in late WAIS and GISP2 data '

%First GISP2
     [C2,ia2,idx2] = unique(GISP2depth,'stable');
     GISP2depthD = accumarray(idx2,GISP2depth,[],@mean); 
     GISP2dateD = accumarray(idx2,GISP2date,[],@mean);
     GISP2gasageD = accumarray(idx2,GISP2gasage,[],@mean); 
     GISP2CH4D = accumarray(idx2,GISP2CH4,[],@mean); 
     GISP2CH4uD = accumarray(idx2,GISP2CH4u,[],@mean); 

%Then WAIS 2018/19
     [C2,ia2,idx2] = unique(WAISnewdepth,'stable');
     WAISnewdepthD = accumarray(idx2,WAISnewdepth,[],@mean); 
     WAISnewdateD = accumarray(idx2,WAISnewdate,[],@mean);
     WAISnewgasageD = accumarray(idx2,WAISnewgasage,[],@mean); 
     WAISnewCH4D = accumarray(idx2,WAISnewCH4,[],@mean); 

%Then WAIS 2012/13
     [C2,ia2,idx2] = unique(WAISolddepth,'stable');
     WAISolddepthD = accumarray(idx2,WAISolddepth,[],@mean); 
     WAISolddateD = accumarray(idx2,WAISolddate,[],@mean);
     WAISoldgasageD = accumarray(idx2,WAISoldgasage,[],@mean); 
     WAISoldCH4D = accumarray(idx2,WAISoldCH4,[],@mean);


%% Separate by year

%Create 0 arrays for each year to compare
NEEM2011 = zeros(1,5);
NEEM2012 = zeros(1,5); 
NEEM2013 = zeros(1,5);
WAIS2012 = zeros(1,4);
WAIS2013 = zeros(1,4);
WAIS2018 = zeros(1,4);
WAIS2019 = zeros(1,4);
GISP2018 = zeros(1,5);
GISP2019 = zeros(1,5);

%Sort by year
%Start with NEEM
for i = 1:length(NEEMdate)
    if NEEMdate(i) == 2010 || NEEMdate(i) == 2011
       NEEM2011 = [NEEM2011; NEEMdate(i),NEEMdepth(i),NEEMgasage(i),NEEMCH4(i),NEEMCH4u(i)];
    end
    if NEEMdate(i) == 2012
       NEEM2012 = [NEEM2012; NEEMdate(i),NEEMdepth(i),NEEMgasage(i),NEEMCH4(i),NEEMCH4u(i)];
    end
    if  NEEMdate(i) == 2013
       NEEM2013 = [NEEM2013; NEEMdate(i),NEEMdepth(i),NEEMgasage(i),NEEMCH4(i),NEEMCH4u(i)];
    end 
end 
    
%Next GISP2 
for i = 1:length(GISP2dateD)
    if GISP2dateD(i) == 2018
       GISP2018 = [GISP2018; GISP2dateD(i),GISP2depthD(i),GISP2gasageD(i),GISP2CH4D(i),GISP2CH4uD(i)];
    end
    if GISP2dateD(i) == 2019
       GISP2019 = [GISP2019; GISP2dateD(i),GISP2depthD(i),GISP2gasageD(i),GISP2CH4D(i),GISP2CH4uD(i)];
    end
end 

%Now WAIS
for i = 1:length(WAISnewdateD)
    if WAISnewdateD(i) == 2018
       WAIS2018 = [WAIS2018; WAISnewdateD(i),WAISnewdepthD(i),WAISnewgasageD(i),WAISnewCH4D(i)];
    end
    if WAISnewdateD(i) == 2019
       WAIS2019 = [WAIS2019; WAISnewdateD(i),WAISnewdepthD(i),WAISnewgasageD(i),WAISnewCH4D(i)];
    end
end 
for i = 1:length(WAISolddateD)
    if WAISolddateD(i) == 2012
       WAIS2012 = [WAIS2012; WAISolddateD(i),WAISolddepthD(i),WAISoldgasageD(i),WAISoldCH4D(i)];
    end
    if  WAISolddateD(i) == 2013
       WAIS2013 = [WAIS2013; WAISolddateD(i),WAISolddepthD(i),WAISoldgasageD(i),WAISoldCH4D(i)];
    end 
end 

%remove first zero row (probably an easier way to do this
NEEM2011(1,:) = [];
NEEM2012(1,:) = []; 
NEEM2013(1,:) = [];
WAIS2012(1,:) = [];
WAIS2013(1,:) = [];
WAIS2018(1,:) = [];
WAIS2019(1,:) = [];
GISP2018(1,:) = [];
GISP2019(1,:) = [];

%% Sort by depth

NEEM2011 = sortrows(NEEM2011,2);
NEEM2012 = sortrows(NEEM2012,2);
NEEM2013 = sortrows(NEEM2013,2);
WAIS2012 = sortrows(WAIS2012,2);
WAIS2013 = sortrows(WAIS2013,2);
WAIS2018 = sortrows(WAIS2018,2);
WAIS2019 = sortrows(WAIS2019,2);
GISP2018 = sortrows(GISP2018,2);
GISP2019 = sortrows(GISP2019,2);

%% Create arrays of all data and sort by depth

%Compile 
NEEMALL = [NEEM2011; NEEM2012; NEEM2013];
GISP2ALL = [GISP2018; GISP2019];
WAISALL = [WAIS2012; WAIS2013; WAIS2018; WAIS2019];
%Sort by age
NEEMALL = sortrows(NEEMALL,2);
GISP2ALL = sortrows(GISP2ALL,2);
WAISALL = sortrows(WAISALL,2);

%% Plot all data 

f = figure(1);
f.Position = [100 100 900 750];
hold on
set(gcf,'color','w');
%Plot All NEEM
a = plot (NEEM2011(:,3)/1000,NEEM2011(:,4),'d','Color',[0.8500 0.3250 0.0980]);
hold on 
b = plot (NEEM2012(:,3)/1000,NEEM2012(:,4),'*','Color',[0.8500 0.3250 0.0980]);
hold on
c = plot (NEEM2013(:,3)/1000,NEEM2013(:,4),'o','Color',[0.8500 0.3250 0.0980]);
%Plot All GISP2
d = plot (GISP2018(:,3)/1000,GISP2018(:,4),'d','Color',[0.3010 0.7450 0.9330]);
hold on
e = plot (GISP2019(:,3)/1000,GISP2019(:,4),'*','Color',[0.3010 0.7450 0.9330]);
hold on
%Plot All WAIS
f = plot (WAIS2012(:,3)/1000,WAIS2012(:,4),'d','Color','b');
hold on
g = plot (WAIS2013(:,3)/1000,WAIS2013(:,4),'*','Color','b');
hold on
h = plot (WAIS2018(:,3)/1000,WAIS2018(:,4),'o','Color','b');
hold on
i = plot (WAIS2019(:,3)/1000,WAIS2019(:,4),'+','Color','b');
hold on
%Plot all data from each core connected
plot (NEEMALL(:,3)/1000,NEEMALL(:,4),'-','Color',[0.8500 0.3250 0.0980])
hold on
plot (GISP2ALL(:,3)/1000,GISP2ALL(:,4),'-','Color',[0.3010 0.7450 0.9330])
hold on
plot (WAISALL(:,3)/1000,WAISALL(:,4),'-','Color','b')
hold on

%Set labels
legend ([a b c d e f g h i],'NEEM 2011','NEEM 2012','NEEM 2013','GISP2 2018','GISP2 2019','WAIS 2012','WAIS 2013','WAIS 2018','WAIS 2019','Location','northwest','EdgeColor','w')
xlabel('Original Gas Age (kyr)','FontSize',18)
ylabel('CH_4 (ppb)','FontSize',18)

%modify graph appearance 
axis ([11420/1000 27070/1000 340 785])
ax = gca;
set (ax,'Xdir','reverse','FontSize',18)
%%
f = figure (2)
f.Position = [100 100 750 900];
hold on
set(gcf,'color','w');
subplot (3,1,1)
hold on
c = plot (NEEM2013(:,3)/1000,NEEM2013(:,4),'o','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerSize',1);
hold on
a = plot (NEEM2011(:,3)/1000,NEEM2011(:,4),'^','MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerSize',1);
hold on 
b = plot (NEEM2012(:,3)/1000,NEEM2012(:,4),'*','MarkerEdgeColor',[0.1660 0.8740 0.1880],'MarkerFaceColor',[0.1660 0.8740 0.1880],'MarkerSize',1);
hold on
% plot (NEEMALL(:,3)/1000,NEEMALL(:,4),'-','Color',[0.8500 0.3250 0.0980])
% hold on
legend ([a b c],'NEEM 2011','NEEM 2012','NEEM 2013','FontSize',13,'EdgeColor','none','Color','none','Location','NorthWest')
xlabel ('NEEM Age (kyr)','FontSize',15)
ylabel ('[CH_4] (ppb)','FontSize',15)
%xlim ([1260 1580])
ylim ([300 800])
set(gca,'FontSize',15)
axis ([11420/1000 27070/1000 340 785])
ax = gca;
set (ax,'Xdir','reverse','FontSize',18)

subplot(3,1,2)
hold on
plot (GISP2018(:,3)/1000,GISP2018(:,4),'d','MarkerEdgeColor',[0.3010 0.7450 0.9330],'MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerSize',1);
hold on
plot (GISP2019(:,3)/1000,GISP2019(:,4),'+','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840],'MarkerSize',1.5);
hold on
% plot (GISP2ALL(:,3)/1000,GISP2ALL(:,4),'-','Color',[0.3010 0.7450 0.9330])
% hold on
legend ('GISP2 2018','GISP2 2019','FontSize',13,'EdgeColor','none','Color','none','Location','NorthWest')
xlabel ('GISP2 Age (kyr)','FontSize',15)
ylabel ('[CH_4] (ppb)','FontSize',15)
set(gca,'FontSize',15)
axis ([11420/1000 27070/1000 340 785])
ax = gca;
set (ax,'Xdir','reverse','FontSize',18)

subplot(3,1,3)
hold on
plot (WAIS2012(:,3)/1000,WAIS2012(:,4),'*','MarkerEdgeColor',[0.1660 0.8740 0.1880],'MarkerFaceColor',[0.1660 0.8740 0.1880],'MarkerSize',1);
hold on
plot (WAIS2013(:,3)/1000,WAIS2013(:,4),'o','MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerSize',1);
hold on
plot (WAIS2018(:,3)/1000,WAIS2018(:,4),'d','MarkerEdgeColor',[0.3010 0.7450 0.9330],'MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerSize',1);
hold on
plot (WAIS2019(:,3)/1000,WAIS2019(:,4),'+','MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840],'MarkerSize',1.5);
hold on
% plot (WAISALL(:,3)/1000,WAISALL(:,4),'-','Color','b')
% hold on
legend ('WAIS 2012','WAIS 2013','WAIS 2018','WAIS 2019','FontSize',13,'EdgeColor','none','Color','none','Location','NorthWest')
xlabel ('WAIS Age (kyr)','FontSize',15)
ylabel ('[CH_4] (ppb)','FontSize',15)
set(gca,'FontSize',15)
axis ([11420/1000 27070/1000 340 785])
ax = gca;
set (ax,'Xdir','reverse','FontSize',18)
annotation('textbox',[.07 .86 .1 .1],'String','(a)','EdgeColor','none','FontSize',20,'FontWeight','bold','Color','k')
annotation('textbox',[.07 .56 .1 .1],'String','(b)','EdgeColor','none','FontSize',20,'FontWeight','bold','Color','k')
annotation('textbox',[.07 .26 .1 .1],'String','(c)','EdgeColor','none','FontSize',20,'FontWeight','bold','Color','k')
annotation('line',[.747 .747],[.12 .92],'LineStyle','-.')