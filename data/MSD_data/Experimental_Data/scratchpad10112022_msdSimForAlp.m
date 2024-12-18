%% To load data. This will take a while
x20=load('20nmAqLSFixedMSD.mat');
x40=load('40nmPFVFixedMSD.mat');
x50=load('50nmVuldiFixedMSD.mat');
x18= load('40nm-18FixedMSD.mat');
x15= load('40nm+15FixedMSD.mat');

%% Finding the fits on the data 
figure
sizeT=60
[a20,b20,c20,rsquare20,sizeV20]=fitLinMSD(x20,15)
[a40,b40,c40,rsquare40,sizeV]=fitLinMSD(x40,sizeT)
[a50,b50,c50,rsquare50,sizeV]=fitLinMSD(x50,sizeT)
[a18,b18,c18,rsquare18,sizeV]=fitLinMSD(x18,sizeT)
[a15,b15,c15,rsquare15,sizeV]=fitLinMSD(x15,sizeT)

%% Calculates the localization error.
sigma40 =mean(c40(5:10))
sigma50=mean(c50(5:10))
sigma18=mean(c18(5:10))
sigma15=mean(c15(5:10))
sigma20=mean(c20(10:13))

%% Constants for simulations
tStep =0.03; %time step, for every data set but 20s should be 0.03, for 20s 0.1
kCycle = tStep;
tsteps =1000; %Steps of random walk
sizeF =tsteps 
numTrajs =20; %number of trajectories to simulate    

%% This section simulates the random walks and calculates the MSD. 
% At the moment it's not really saving individual trajectories for simulated data, but you can add that if you like
    counter =1;
    %clear tempFull
    %tempX = zeros(sizeF,1,numTrajs);
    %tempY = zeros(sizeF,1,numTrajs);
    %tempZ = zeros(sizeF,1,numTrajs);
    %tempFull = zeros(sizeF,1,numTrajs);

    clear tempX
    clear tempY
    clear tempZ

    
for iiFile=1:numTrajs 
    iiFile
    %try
 
    
    scratchpad09292022_simulatingRandomWalkStickCaps %Here this calls the simulation script, change the name to whichever code you want to run
     %clear msdX
     %clear msdY
     %clear msdZ
     %clear msdFull
     %clear zeroed
  
    zeroed = [Xprime,Yprime,Zprime, (1:tsteps)'];
    lastStep = zeroed(end,4);
  
    for iiSpacing =1:sizeF
    
    %iiSpacing

            counterTemp =1;

            for iiPossible = 0:iiSpacing-1
                indexes = (iiPossible:iiSpacing:lastStep)';
                %indexes

                for iiIndex =2:size(indexes,1)
                    if(any(zeroed(:,4)==indexes(iiIndex))&& any(zeroed(:,4)==indexes(iiIndex-1)))
                        tempX(iiSpacing,counterTemp,iiFile) = (zeroed(zeroed(:,4)==indexes(iiIndex),1)-zeroed(zeroed(:,4)==indexes(iiIndex-1),1)).^2;
                        tempY(iiSpacing,counterTemp,iiFile) = (zeroed(zeroed(:,4)==indexes(iiIndex),2)-zeroed(zeroed(:,4)==indexes(iiIndex-1),2)).^2;
                        tempZ(iiSpacing,counterTemp,iiFile) = (zeroed(zeroed(:,4)==indexes(iiIndex),3)-zeroed(zeroed(:,4)==indexes(iiIndex-1),3)).^2;
                        tempFull(iiSpacing,counterTemp,iiFile) = tempX(iiSpacing,counterTemp,iiFile)+tempY(iiSpacing,counterTemp,iiFile)+tempZ(iiSpacing,counterTemp,iiFile);


                        counterTemp =counterTemp +1;
                    end 
                    
                end

            end         
    end

   tempX(:,:,iiFile);
     %catch
     %end
end


Mx=[];
My=[];
Mz=[];
Mfull=[];
for iiFile = 1:numTrajs

        Mx = [Mx,tempX(:,:,iiFile)];
        My = [My,tempY(:,:,iiFile)];
        Mz = [Mz,tempZ(:,:,iiFile)];
        Mfull = [Mfull,tempFull(:,:,iiFile)];

end



Mx(Mx==0)=NaN;
sx = nanstd(Mx,[],2); 
rowMeanX = nanmean(Mx,2); %MSDx
My(Mz==0)=NaN;
sy = nanstd(My,[],2); 
rowMeanY = nanmean(My,2); %MSDy
Mz(Mz==0)=NaN;
sz = nanstd(Mz,[],2); 
rowMeanZ = nanmean(Mz,2); %MSDz
Mfull(Mfull==0)=NaN;
sf = nanstd(Mfull,[],2); 
rowMeanFull = nanmean(Mfull,2); %MSDfull




%% Plot individual MSD for each trajectory
figure
for iiFile = 1: numTrajs
tempTraj = tempFull(:,:,iiFile);
tempTraj(tempTraj==0)=NaN;
trajMSD = nanmean(tempTraj,2); 
hold on
plot(trajMSD(1:500))
end

set(gca,'xscale','log')
set(gca,'yscale','log')




%% Figure averaged x,y,z,full
figure, plot(1:size(rowMeanX),rowMeanX)
hold on, plot(1:size(rowMeanY),rowMeanY)
plot(1:size(rowMeanZ),rowMeanZ)
 plot(1:size(rowMeanFull),rowMeanFull)

%% Figure average x, y ,z full all start at zero 
figure, plot(1:size(rowMeanX),rowMeanX-min(rowMeanX))
hold on, plot(1:size(rowMeanY),rowMeanY-min(rowMeanY))
plot(1:size(rowMeanZ),rowMeanZ-min(rowMeanZ))
 plot(1:size(rowMeanFull),rowMeanFull-min(rowMeanFull))
 
 %% MSD in micrometer instead of nanometers
figure, plot((1:size(rowMeanX))*kCycle,rowMeanX*10^-6)
hold on, plot((1:size(rowMeanX))*kCycle,rowMeanY*10^-6)
plot((1:size(rowMeanX))*kCycle,rowMeanZ*10^-6)
 plot((1:size(rowMeanX))*kCycle,rowMeanFull*10^-6)
        

%% Here what I did was run the above MSD simulation with different versions of the random walk, 
% saved the MSD result(rowMeanFull) in a varaible eg simRowMeanFull5ReDraw
% and then this lets me plot all of the different plots in one figure 
figure;hold on
kCycle =0.03 %% remember 0.03 for 40,50,18 and 15, and 0.1 for 20 

t4Fit=((1:250)*kCycle);
 xFit2 =(x15.rowMeanFull(1:250)-sigma15) %sigma is the localization error, it is calculate in the data section at the beggining
 plot(t4Fit,xFit2,'k.-')


t4Fit=((1:250)*kCycle);
xFit = (rowMeanFull(1:250))
plot(t4Fit,xFit,'m.-')
t4Fit=((1:250)*kCycle);
xFit = (simRowSave(1:250))
plot(t4Fit,xFit,'y.-')

  t4Fit=((1:250)*0.03)
  plot (t4Fit,simRowMeanFull15ReDraw(1:250),'r.-')
    t4Fit=((1:250)*0.03)
  plot (t4Fit,simRowMeanFull15Sticky(1:250),'g.-')
    t4Fit=((1:250)*0.03)
  plot (t4Fit,simRowMeanFull15FullCell(1:250),'b.-')

%%Plot in log scale 
%set(gca,'xscale','log') 
%set(gca,'yscale','log')
%% You can use this for plotting stuff, it's convenient to use with the curve fitting app in matlab
%To access the curve fitting app in matlab click on apps and then curve
%fitting you can select the data and then apply different fits 

figure;hold on
t4Fit=((1:15)*kCycle);
xFit = (simRowMeanFull50ReDraw(1:15));
plot(t4Fit,xFit,'r.--')

t4Fit=((1:15)*kCycle);
xFit2 =(x50.rowMeanFull(1:15));
plot(t4Fit,xFit2,'k.-')
%set(gca,'xscale','log')
%set(gca,'yscale','log')


%% I used this to find the diffusion coefficients we are using for the simulations, you will probably won't need it. 
% For the fit in the figure you can go to the top menu, select tools, and
% in the bottom it says basic fitting, that opens a nice little tool to do
% linear fitting and some easy stuff like that . 
figure;hold on

t4Fit=((3:8)*0.1);
xFit2 =(x20.rowMeanFull(3:8)-sigma20);
plot(t4Fit,xFit2,'r.-')

figure; hold on
t4Fit=((1:5)*0.1);
xFit2 =(x20.rowMeanFull(1:5)-sigma20);
plot(t4Fit,xFit2,'k.-')
%set(gca,'xscale','log')
%set(gca,'yscale','log')


%% bootstrap for error bars (not necessary to run the other stuff)
intervalsX =bootci(200,{@nanmean,Mx'})';
intervalsY =bootci(200,{@nanmean,My'})';
intervalsZ =bootci(200,{@nanmean,Mz'})'
intervalsF =bootci(200,{@nanmean,Mfull'})';

%% This plots stuff with error bars requires the bootstrap step above 
figure
colorX = [55,126,184]/255; %yz
colorY =[77,175,74]/255;%xz
colorZ =[228,26,28]/255;%xy

sizeT =	600;
% figure, plot(log(kCycle*(1:sizeT)),log(rowMeanX(1:sizeT)*10^-6),'o','MarkerSize',6,'LineWidth',2,'Color' , colorX)
% hold on 
% lsline
%  plot(log(kCycle*(1:sizeT)),log(rowMeanY(1:sizeT)*10^-6),'o','MarkerSize',6,'LineWidth',2,'Color' , colorY)
% lsline
%  plot(log(kCycle*(1:sizeT)),log(rowMeanZ(1:sizeT)*10^-6),'o','MarkerSize',6,'LineWidth',2, 'Color' , colorZ)
% lsline


plot(log(kCycle*(1:sizeT)),log(rowMeanFull(1:sizeT)*10^-6),'ok','MarkerSize',6,'LineWidth',2,'MarkerEdgeColor','b')
lsline
figure
%ylim([0,0.06]);
errorbar((1:sizeT)*kCycle,rowMeanFull(1:sizeT)*10^-6,intervalsF(1:sizeT,2)*10^-6-rowMeanFull(1:sizeT)*10^-6,rowMeanFull(1:sizeT)*10^-6-intervalsF(1:sizeT,1)*10^-6,...
    '.k-','LineWidth',2)
%title('MSD')
set(gca, 'FontSize', 20)%,'BinWidth', 50)
xlabel('\tau (s)')
ylabel('MSD (\mum^2)')
title('50 nm MSD')
set(gca,'xscale','log')
set(gca,'yscale','log')
%plot ((1:sizeT)*kCycle,(1:sizeT)*kCycle,'g--')



x =  log(kCycle*(1:sizeT))';
y =  log(rowMeanFull(1:sizeT)*10^-6);
set(gca, 'FontSize', 20)%,'BinWidth', 50)
xlabel('log (\tau)')
ylabel('log (MSD)')
p = polyfit(x,y,1)
yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal
title (strcat ('MSD 40nm = Dt^{\alpha}, \alpha= ',num2str(p(1))) )
p(2)

set(gca,'FontSize',24)
set(gca,'LineWidth',2.5)


