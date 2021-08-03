% setup
clear all; close all
%% Set Parameters
% Time parameters
tStart = 0;
tStop = 20000;
dt = 1;

% EtOH Dose Parameters
doseStep = 0.5;
doseMin = 1;
doseMax = 5;
EtOH = doseMin:doseStep:doseMax;

% Noise
includeNoise = 0; % 1 or 0 Depending on Desire to Include Noise. 1 = Noise. 0 = No Noise.

% Run paramscript containing parameters
% Simulations with noise on require a different parameter set than
% simulations without noise.
% Sims with Noise on: Theta = 4.8, Beta = 0.7
% Sims with Noise off: Theta = -1, Beta = 7
if includeNoise == 0
    paramscript_withoutNoise
elseif includeNoise == 1
    paramscript_withNoise
else
    paramscript_withoutNoise
end

%% Simulate over Doses of EtOH
for i = 1:length(EtOH)
    TSpan = tStart:dt:tStop;
    noise = createNoise(params.sigN, 0, (tStop), (tStop)/dt, 0, 1); % Noise is created here. Set to 0 if Noise = 0. createNoise is a custom function
    if includeNoise == 0
        noise = noise*0;
    elseif includeNoise ~= 1
        noise = noise*0;
    end
    params.Adesp = params.Adesp / EtOH(i);
    intBE = griddedInterpolant(TSpan,noise(1,:));   % Interpolate noise to fit within RK4
    intBI = griddedInterpolant(TSpan,noise(2,:));   
    
    params.EtOH = EtOH(i); %Iteratively set EtOH Parameter

    
    odefun = @(T, Y) [(-Y(1) + params.gE * binary_stepE((params.Jee * Y(1) - params.Jei * Y(2) - Y(3) + intBE(T) + params.thetaE),params.Edesp))/(params.tauE); ...
        (-Y(2) + params.gI * binary_stepI((params.Jie * Y(1) - params.Jii * Y(2) + intBI(T) + params.thetaI),params.Idesp))/(params.tauI); ...
        (-Y(3) * (1/params.EtOH) + (params.betaE) * Y(1))/(params.tauA)];
    
    % ODEFUN is the set of equations from Jercog 2017.
    
    y0 = [0;0;2];   %Initial conditions. rE is index 1, rI is index 2, rA is index 3
                    %Adaptation current is set to 2 
                    
    [Y] = RK4(odefun, TSpan, y0);   %Integrate using custom RK4 function 
    % Compare results
    res{i} = Y;

    figure; hold on;
    plot(TSpan, Y, '-', 'linewidth', 2);
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')
    xlabel('Time (ms)')
    ylabel('Rates (Hz)')
    title(['Model Rates at EtOH Concentration ' num2str(EtOH(i))])
    if includeNoise == 0
        filename = ['figures/withoutNoise/RatesOverTime' num2str(i)];
        saveas(gcf, filename, 'fig')
        saveas(gcf, filename, 'svg')
        saveas(gcf, filename, 'png')
    elseif includeNoise == 1
        filename = ['figures/withNoise/RatesOverTime' num2str(i)];
        saveas(gcf, filename, 'fig')
        saveas(gcf, filename, 'svg')
        saveas(gcf, filename, 'png')
    else 
        disp('No valid flag provided, figures unsaved')
    end

    %% Calculate UDS
    % Logical statement is required here due to noise creating up or down
    % states that are less than 50 ms in duration. Jercog 2017 adds up/down
    % states less than 50 ms to neighboring state durations. Here I am
    % smoothing the data in the noise condition to accomplish the same thing. 
    
    if includeNoise == 0
        up_idx = (Y(:,1)) > 1;
        down_idx = (Y(:,1)) < 1;
    else
        up_idx = smooth(Y(:,1),25) > 1;
        down_idx = smooth(Y(:,1),25) < 1;
    end
    
    labeledMatrix1 = bwlabel(up_idx);
    measurements1 = regionprops(labeledMatrix1, 'Area');
    upDurations = [measurements1.Area];

    labeledMatrix2 = bwlabel(down_idx);
    measurements2 = regionprops(labeledMatrix2, 'Area');
    downDurations = [measurements2.Area];

    minDuration = 1;


    allupDurations{i} = upDurations(upDurations > minDuration);
    alldownDurations{i} = downDurations(downDurations > minDuration);

    CV_U{i} = sqrt(var(upDurations))/length(upDurations);
    CV_D{i} = sqrt(var(downDurations))/length(downDurations);

end
%%
if includeNoise == 0
    filename = ['figures/withoutNoise/results.mat'];
    save(filename,'res')
elseif includeNoise == 1
    filename = ['figures/withNoise/results.mat'];
    save(filename,'res')
else
    disp('No valid flag provided')
end
close all
%% Analyze
% Calculate UDS Statistics Over EtOH "Dose"
% In the for loop below the first and last down state or up state are
% removed in case they are cut off by tStart or tStop times.

meanalldownDurations = zeros(length(alldownDurations),1);
meanallupDurations = zeros(length(allupDurations),1);

for i = 1:length(alldownDurations)
    meanalldownDurations(i) = median(alldownDurations{i}(1+1:end-1));
    meanallupDurations(i) = median(allupDurations{i}(1+1:end-1));
    semalldownDurations(i) = std(alldownDurations{i}(1+1:end-1))/sqrt(length(alldownDurations{i}(1+1:end-1)));
    semallupDurations(i) = std(allupDurations{i}(1+1:end-1))/sqrt(length(allupDurations{i}(1+1:end-1)));
end

figure(1)
plot(EtOH,meanalldownDurations,'bo','MarkerSize',8)
hold on
plot(EtOH,meanallupDurations,'ro','MarkerSize',8)
errorbar(EtOH,meanalldownDurations,semalldownDurations)
errorbar(EtOH,meanallupDurations,semallupDurations)
xlim([0.8 5.2])
xlabel('EtOH Concentration Variable')
ylabel('Duration (ms)')
title(['Mean Up and Down State Durations' 'Noise = ' num2str(includeNoise)])
% legend('Down State Durations', 'Up State Durations')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')

if includeNoise == 0
    filename = ['figures/withoutNoise/UDSDurPerDose'];
    print('-dpng', filename);
    saveas(gcf, filename, 'fig')
    saveas(gcf, filename, 'svg')
elseif includeNoise == 1
    filename = ['figures/withNoise/UDSDurPerDose'];
    print('-dpng', filename);
    saveas(gcf, filename, 'fig')
    saveas(gcf, filename, 'svg')
else
    disp('No valid flag provided')
end

%% Plot Adaptation Current Over EtOH Dose
if includeNoise == 0
%     evalDoseWindow = [1 1.5 2 2.5 3 3.5 4 4.5 5];
%     refData = res{evalDoseWindow(end)};
%     [~,maxpeakAll] = findpeaks(refData(:,3));
%     maxWindow = maxpeakAll(2) - maxpeakAll(1);
%     lenMax = length(1:maxWindow);
%     [~,refpeakAll] = findpeaks(res{1}(:,3));
%     refPeak = refpeakAll(1);
% 
%     figure(2); hold on
% 
%     for i = evalDoseWindow
%         presRes = res{i};
%         [~,presPeaks] = findpeaks(presRes(1:maxWindow,3));
%         presRefPeak = presPeaks(1);
%         winStart = presRefPeak - refPeak + 1;
%         winEnd = winStart + lenMax - 1;
%         plot(1:maxWindow, presRes(winStart:winEnd,3),'linewidth', 2)
%     end
%     legend(num2str(evalDoseWindow'))
%     xlabel('Time (ms)')
%     ylabel('Adaptation Current')
%     set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')
% 
%     % Plot Adaptation Current Over EtOH Dose - First Peak Only
% 
% 
%     figure(3); hold on
% 
%     for i = evalDoseWindow
%         presRes = res{i};
%         [~,presPeaks] = findpeaks(presRes(:,3));
%         presRefPeak = presPeaks(1);
%         presWin = presPeaks(2) - presPeaks(1);
%         plot(1:presWin, presRes(1:presWin,3),'linewidth', 2)
%     end
%     legend(num2str(evalDoseWindow'))
%     xlabel('Time (ms)')
%     ylabel('Adaptation Current')
%     set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')


    % Plot Adaptation Current Over EtOH Dose - First Peak Only - Aligned


    evalDoseWindow = [1 2 3 4 5];
    refData = res{evalDoseWindow(end)};
    [~,maxpeakAll] = findpeaks(refData(:,3));
    maxWindow = maxpeakAll(2) - maxpeakAll(1);
    lenMax = length(1:maxWindow);
    [~,refpeakAll] = findpeaks(res{1}(:,3));
    refPeak = refpeakAll(1);

    figure(4); hold on

    for i = evalDoseWindow
        presRes = res{i};
        [~,presPeaks] = findpeaks(presRes(:,3));
        presRefPeak = presPeaks(1);
        peakWindow = presPeaks(2) - presPeaks(1);
        presRes(peakWindow:end,3) = NaN;
        winStart = presRefPeak - refPeak + 1;
        winEnd = winStart + lenMax - 1;
        plot(1:maxWindow, presRes(winStart:winEnd,3),'linewidth', 2)
    end
    legend(num2str(evalDoseWindow'))
    xlabel('Time (ms)')
    ylabel('Adaptation Current')
    set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')

    if includeNoise == 0
        filename = ['figures/withoutNoise/AdaptPeakAligned'];
        saveas(gcf, filename, 'png')
        saveas(gcf, filename, 'svg')
        saveas(gcf, filename, 'fig')
    else
        filename = ['figures/withNoise/AdaptPeakAligned'];
        saveas(gcf, filename, 'png')
        saveas(gcf, filename, 'svg')
        saveas(gcf, filename, 'fig')
    end

end
%% Phase Diagrams

% These can be run independent of any results
% Require some parameter changes 
odefun = @(T, Y) [(-Y(1) + params.gE * binary_stepE((params.Jee * Y(1) - params.Jei * Y(2) - Y(3) + intBE(T) + params.thetaE),params.Edesp))/(params.tauE); ...
    (-Y(2) + params.gI * binary_stepI((params.Jie * Y(1) - params.Jii * Y(2) + intBI(T) + params.thetaI),params.Idesp))/(params.tauI); ...
    (-Y(3) * (1/params.EtOH) + (params.betaE * Y(1)))/(params.tauA)];

% Figure 4, Panel B
evalWindow = 40;
re = 0:1:evalWindow-1;
ri = 0:1:evalWindow-1;
ra = 0:1:evalWindow-1;
ra = ra*0;

E = zeros(evalWindow,1)';
I = zeros(evalWindow,1)';
params.thetaE = 0;
params.thetaI = 0;
for z = 1:evalWindow
   E(z) = params.gE * binary_stepE(z + params.thetaE, params.Edesp);
   I(z) = params.gI * binary_stepI(z + params.thetaI, params.Idesp);
end
 figure(5) 
 plot(E,'r-','LineWidth',3); hold on; plot(I,'b-','LineWidth',3)
 xlabel('Input z (a.u.)')
 ylabel('Rate ?_x(z)')
 set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')

 % Figure 4, Panel C
params.thetaE = 0;
params.thetaI = 0;
evalWindow = 20;
dt = 0.01;
re = -5:dt:evalWindow;
ri = -5:dt:evalWindow;
ra = 1;


nullE = params.gE * binary_stepE((params.Jee * re - params.Jei * ri - ra + params.thetaE),params.Edesp);
nullI = params.gI * binary_stepI((params.Jie * re - params.Jii * ri + params.thetaI),params.Idesp);

% quiver(re,ri,nullE,nullI)
nullE(re < 0) = nullE(re < 0) * -1;
re(re < 0) = 0;
% nullE(nullE < 0) = nullE(nullE < 0) * -1;
nullI(nullI < 0) = 0;
figure(6)
plot(re,nullE,'r-','LineWidth',3); hold on; plot(ri,nullI,'b-','LineWidth',3)
xlabel('Rate r_e (Hz)')
ylabel('Rate r_i (Hz)')
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold')

ylim([-6 15])
xlim([-0.3 4.3])

