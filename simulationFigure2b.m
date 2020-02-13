% Matlab script used to generate Figure 2b of the article:
%
% "Low-Complexity Distributed XL-MIMO for Multiuser Detection"
%
% @author Victor Croisfelt Rodrigues
% @author Abolfazl Amiri
% @author Taufik Abrao
% @author Elisabeth de Carvalho
% @author Petar Popovski
%

%Initialization
close all;
clearvars;

%Fixing random seed
rng(0)

%% Loading simulation setup
load('20191217_CRD_Norm1_M100_S4_K10',...
    'M','S','Ms','K','diagNorm',...
    'p','noiseVariancedBmRange',...
    'updateSchedule','meanConv',...
    'numSetups','numRealizations')

%% Algorithm parameters

%Define number of iterations vector
numIterations = meanConv;

%Number of bounds
bounds = 2;

%% Communication parameters

%Define modulation order
mOrder = 4;

%Data fusion methods
dataFusion = ["DEDF","DLDF"];

%% Simulation

%Prepare to save simulation results
SER_RZF = zeros(length(noiseVariancedBmRange),length(dataFusion),numSetups);
SER_RKA = zeros(length(noiseVariancedBmRange),bounds,length(updateSchedule),length(dataFusion),numSetups);

%Go through all setups
for t = 1:numSetups

    setup = tic; % setting up a timer: setup
    warning('off','MATLAB:sqrtm:SingularMatrix') % omitting warning related to sqrtm function

    %Output setup progress
    disp([num2str(t) ' setups out of ' num2str(numSetups)]);

    %Prepare to save parfor simulation results
    par_SER_RZF = zeros(length(noiseVariancedBmRange),length(dataFusion));
    par_SER_RKA = zeros(length(noiseVariancedBmRange),bounds,length(updateSchedule),length(dataFusion));

    %Generate mobile communication setup
    [channelGaindB,Rdistrib] = functionExampleSetup(M,S,K,diagNorm);

    %Go through all different SNR values
    for r = 1:length(noiseVariancedBmRange)

        %Output SNR progress
        disp([num2str(r) ' SNRs out of ' num2str(length(noiseVariancedBmRange))]);

        %Extract current noise variance
        noiseVariancedBm = noiseVariancedBmRange(r);

        %Compute channel gain over noise
        channelGainOverNoise = zeros(Ms,K,S);
        channelGainOverNoise(channelGaindB ~= 0) = channelGaindB(channelGaindB ~= 0) - noiseVariancedBm;

        %Generate random signal sent by each one of the K users
        x = randi([0 mOrder-1],[K 1]);

        %Modulate
        xMod = qammod(x,mOrder,'gray','UnitAveragePower',true);

        %Tracking some values
        xHatMod_RZF = zeros(K,numRealizations,S);
        xHatMod_RKA = zeros(K,numRealizations,S,bounds,length(updateSchedule));

        %Go through each subarray
        for ss = 1:S

            %Channel realizations
            H = functionChannelRealizations(Ms,K,channelGainOverNoise(:,:,ss),Rdistrib(:,:,:,ss),numRealizations);

            %Received signal
            Y = functionReceivedSignal(Ms,K,p,numRealizations,xMod,H);

            %RZF receive combining matrix
            V_RZF = functionRZF(Ms,K,p,numRealizations,H);

            %Obtain RZF signal estimates
            xHatMod_RZF(:,:,ss) = functionSignalEstimate(Ms,K,numRealizations,Y,V_RZF);

            %Go through all different USs
            for us = 1:length(updateSchedule)

                %Go through all the bounds
                for b = 1:bounds

                    %Run RKA
                    V_RKA = functionRKA(Ms,K,p,numRealizations,numIterations(r,b,us),Rdistrib(:,:,:,ss),H,updateSchedule(us));

                    %Obtain RKA signal estimates
                    xHatMod_RKA(:,:,ss,b,us) = functionSignalEstimate(Ms,K,numRealizations,Y,V_RKA);

                end

            end

        end

        %Go through all different data fusion methods
        for d = 1:length(dataFusion)

            %RZF data fusioning
            xHatMod_RZF_fusioned = functionDataFusion(S,K,channelGainOverNoise,numRealizations,xHatMod_RZF,dataFusion(d));

            %Demodulating
            xHat_RZF = qamdemod(xHatMod_RZF_fusioned,mOrder,'gray','UnitAveragePower',true);

            %Storing results
            par_SER_RZF(r,d) = sum(sum(xHat_RZF ~= x,2))/numRealizations/K;

            %Go through all different USs
            for us = 1:length(updateSchedule)

                %Go through all the bounds
                for b = 1:bounds

                    %RKA data fusioning
                    xHatMod_RKA_fusioned = functionDataFusion(S,K,channelGainOverNoise,numRealizations,xHatMod_RKA(:,:,:,b,us),dataFusion(d));

                    %Demodulating
                    xHat_RKA = qamdemod(xHatMod_RKA_fusioned,mOrder,'gray','UnitAveragePower',true);

                    %Storing results
                    par_SER_RKA(r,b,us,d) = sum(sum(xHat_RKA ~= x,2))/numRealizations/K;

                end

            end

        end

    end

    %Save parfor simulation results
    SER_RZF(:,:,t) = par_SER_RZF;
    SER_RKA(:,:,:,:,t) = par_SER_RKA;

    toc(setup)

end

%% Data extraction
meanSER_RZF = mean(SER_RZF,3);
meanSER_RKA = mean(SER_RKA,5);

%% Save simulation results
save([date,'_SER_',diagNorm,'_M',num2str(M),'_S',num2str(S),'_K',num2str(K)])

%% Plotting simulation results

figure;
subplot(2,2,1)
hold on; box on;

plot(noiseVariancedBmRange,meanSER_RZF(:,1),'k','LineWidth',1)

plot(noiseVariancedBmRange,meanSER_RKA(:,1,1,1),'k*','LineWidth',1)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,2,1),'kd','LineWidth',1)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,3,1),'kh','LineWidth',1)

plot(noiseVariancedBmRange,meanSER_RKA(:,1,1,1),'*--','LineWidth',1)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,2,1),'d--','LineWidth',1)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,3,1),'h--','LineWidth',1)

xlabel('Noise variance [dBm]');
ylabel('SER');

set(gca,'YScale','log')

xlim([min(noiseVariancedBmRange) max(noiseVariancedBmRange)])

legend('RZF benchmark','Power-based update schedule','Uniform update schedule','Active-antenna-based updated schedule','Location','Northwest');

title('Losing 10\% of RZF performance: Distributed equal data fusion (DEDF)')

subplot(2,2,2)
hold on; box on;

plot(noiseVariancedBmRange,meanSER_RZF(:,2),'k','LineWidth',1)

plot(noiseVariancedBmRange,meanSER_RKA(:,1,1,2),'k*','LineWidth',1)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,2,2),'kd','LineWidth',1)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,3,2),'kh','LineWidth',1)

plot(noiseVariancedBmRange,meanSER_RKA(:,1,1,2),'*--','LineWidth',1)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,2,2),'d--','LineWidth',1)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,3,2),'h--','LineWidth',1)

xlabel('Noise variance [dBm]');
ylabel('SER');

set(gca,'YScale','log')

xlim([min(noiseVariancedBmRange) max(noiseVariancedBmRange)])

legend('RZF benchmark','Power-based update schedule','Uniform update schedule','Active-antenna-based updated schedule','Location','Northwest');

title('Losing 10\% of RZF performance: Distributed linear data fusion (DLDF)')

subplot(2,2,3)
hold on; box on;

plot(noiseVariancedBmRange,meanSER_RZF(:,1),'k','LineWidth',1.5)

plot(noiseVariancedBmRange,meanSER_RKA(:,2,1,1),'k*','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,2,2,1),'kd','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,2,3,1),'kh','LineWidth',1.5)

plot(noiseVariancedBmRange,meanSER_RKA(:,2,1,1),'*-.','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,2,2,1),'d-.','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,2,3,1),'h-.','LineWidth',1.5)

xlabel('Noise variance [dBm]');
ylabel('SER');

set(gca,'YScale','log')

xlim([min(noiseVariancedBmRange) max(noiseVariancedBmRange)])

legend('RZF benchmark','Power-based update schedule','Uniform update schedule','Active-antenna-based updated schedule','Location','Northwest');

title('Losing 1\% of RZF performance: Distributed equal data fusion (DEDF)')

subplot(2,2,4)
hold on; box on;

plot(noiseVariancedBmRange,meanSER_RZF(:,2),'k','LineWidth',1.5)

plot(noiseVariancedBmRange,meanSER_RKA(:,2,1,2),'k*','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,2,2,2),'kd','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,2,3,2),'kh','LineWidth',1.5)

plot(noiseVariancedBmRange,meanSER_RKA(:,2,1,2),'*-.','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,2,2,2),'d-.','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,2,3,2),'h-.','LineWidth',1.5)

xlabel('Noise variance [dBm]');
ylabel('SER');

set(gca,'YScale','log')

xlim([min(noiseVariancedBmRange) max(noiseVariancedBmRange)])

legend('RZF benchmark','Power-based update schedule','Uniform update schedule','Active-antenna-based updated schedule','Location','Northwest');

title('Losing 1\% of RZF performance: Distributed linear data fusion (DLDF)')

%%

figure(1);
hold on; box on; ax = gca;

plot(noiseVariancedBmRange,meanSER_RZF(:,2),'k','LineWidth',2)

plot(noiseVariancedBmRange,meanSER_RKA(:,1,1,2),'k*','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,2,2),'kd','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,3,2),'kh','LineWidth',1.5)

ax.ColorOrderIndex = 1;

plot(noiseVariancedBmRange,meanSER_RKA(:,1,1,2),'*--','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,2,2),'d--','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,1,3,2),'h--','LineWidth',1.5)

%xlabel('Noise variance [dBm]');
%ylabel('Average SER');

%set(gca,'YScale','log')

%ylim([1e-4 1])
%xlim([min(noiseVariancedBmRange) max(noiseVariancedBmRange)])

%legend('RZF benchmark','Power-based update schedule','Uniform update schedule','Active-antenna-based updated schedule','Location','Northwest');

%title('Norm. 2 \& Losing 10\%')

%subplot(1,2,2)
%hold on; box on; ax = gca;

%plot(noiseVariancedBmRange,meanSER_RZF(:,2),'k','LineWidth',1.5)

%plot(noiseVariancedBmRange,meanSER_RKA(:,2,1,2),'k*','LineWidth',1.5)
%plot(noiseVariancedBmRange,meanSER_RKA(:,2,2,2),'kd','LineWidth',1.5)
%plot(noiseVariancedBmRange,meanSER_RKA(:,2,3,2),'kh','LineWidth',1.5)

ax.ColorOrderIndex = 1;

plot(noiseVariancedBmRange,meanSER_RKA(:,2,1,2),'*-.','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,2,2,2),'d-.','LineWidth',1.5)
plot(noiseVariancedBmRange,meanSER_RKA(:,2,3,2),'h-.','LineWidth',1.5)

xlabel('Noise variance [dBm]');
ylabel('Average SER');

set(gca,'YScale','log')

ylim([1e-3 .3])
xlim([min(noiseVariancedBmRange) max(noiseVariancedBmRange)])

legend('RZF benchmark','Power','Uniform','Active-antennas','Location','best');

title('Normalization 1 of $\mathbf{D}_{k}$')
