% EENG 311- Final Project
% Group 3
% December 15, 2021
% This code creates a survival curve from randomly generated data using
% both the Kaplan-Meier Estimator and the Nelson-Aalen Estimator. This code
% compares how the two behave when there is both censored and non-censored
% data.
%% Non-censored Estimations
% define number of participants
n = 100;
% define for start and end times
data = NaT(n+1, 2);
% define for id #, time, and status
data2 = zeros(n+1, 2);
atRisk = zeros(1, n+1);
survive = zeros(1, n+1);
% index of participant
for i = 2:n+1
    data2(i, 1) = i;
end

% loop that assigns diagnosis and death dates
for i = 1:n
    % decide random diagnosed date
    day1 = randi([1 31],1,1);
    month1 = randi([1 12],1,1);
    year1 = randi([2019 2021],1,1);
    data(i+1, 1) = datetime(year1, month1, day1);
    
    % define days and months above if statements
    month2 = 0;
    day2 = 0;
    % find death year first
    year2 = randi([year1 2021],1,1);
    if year2 == year1
        % if years equal, month must be same or later
        month2 = randi([month1 12],1,1);
        % if years and months equal, day must be same or later
        if month2 == month1
            day2 = randi([day1 31],1,1);
        end
    else
        month2 = randi([1 12],1,1);
        day2 = randi([1 31],1,1);
    end
    % random death date
    data(i+1, 2) =  datetime(year2, month2, day2);
end

% find amount of time between events
for i = 2:n+1
    data2(i, 2) = days(data(i, 2) - data(i, 1));
end

% sort the data based on time between events
sortedData = sortrows(data2, 2);

% sort the event dates relative to time between events
sortedDates = NaT(n+1, 2);
for i = 2:n+1
   x = sortedData(i, 1);
   sortedDates(i, 1) = data(x, 1);
   sortedDates(i, 2) = data(x, 2);
end

% set at risk to be everyone, and P(survive) = 1 to start
atRisk(1) = n;
survive(1) = 1;

% use Kaplan-Meier to estimate survival probability
for i = 2:n
    j = i;
    numDeath = 0;
    % if multiple events happened in the same amount of time
    if sortedData(j, 2) == sortedData(j+1, 2)
        % count how many events occured in same time
        while sortedData(j, 2) == sortedData(j+1, 2)
            % if there was a death
            numDeath = numDeath + 1;
            atRisk(j) = atRisk(j - 1) - 1;
            if j <= 1000
                j = j + 1;
            end
        end
        survive(i) = survive(i - 1) * (1 - (numDeath / atRisk(i)));
    else
        numDeath = numDeath + 1;
        atRisk(i) = atRisk(i - 1) - 1;
        survive(i) = survive(i - 1) * (1 - (numDeath / atRisk(i)));
    end
    i = j;
    numDeath = 0;
end

% Estimate hazard function with Nelson-Aalen Estimator
hazard = zeros(1, n+1);
atRisk(1) = n;
for i = 2:n    
    j = i;
    numDeath = 0;
    % if multiple events happened in the same amount of time
    if sortedData(j, 2) == sortedData(j+1, 2)
        % count how many events occured in same time
        while sortedData(j, 2) == sortedData(j+1, 2)
            % if there was a death
            numDeath = numDeath + 1;
            atRisk(j) = atRisk(j - 1) - 1;
            if j <= n
                j = j + 1;
            end
        end
        % estimate hazard probability
        hazard(i) = hazard(i-1) + (numDeath / atRisk(i));
    else
        % check for censoring
        numDeath = numDeath + 1;
        atRisk(i) = atRisk(i - 1) - 1;            
        % estimate hazard probability
        hazard(i) = hazard(i-1) + (numDeath / atRisk(i));
    end
    i = j;
    numDeath = 0;
    % check for the last element
    if i == n
        numDeath = numDeath + 1;
        atRisk(i + 1) = atRisk(i) - 1;            
        hazard(i + 1) = hazard(i) + (numDeath / atRisk(i + 1));
    end
end

% Estimate a survival curve from hazard function
sFromH = zeros(1, n+1);
for i = 1:n+1
    sFromH(i) = exp(-1 * hazard(i));
end

% compare two methods of finding survival curve
figure(1)
hold on
pl4 = plot(sortedData(:, 2), sFromH);
pl3 = plot(sortedData(:, 2), survive);
title('$\hat{S}$(t) from Kaplan-Meier vs $\hat{S}$(t) from Nelson-Aalen', 'Interpreter','latex');
ylabel('P(t)');
xlabel('Time (days)');
h = [pl4 pl3];
legend(h, '$\hat{S}$(t) from Kaplan-Meier', '$\hat{S}$(t) from Nelson-Aalen', 'Interpreter','latex');
hold off

%% Censored Estimations
% define for start and end times
dataCensor = NaT(n+1, 2);
% define for id #, time, and status
dataWithCensoring2 = zeros(n+1, 3);
atRiskCensor = zeros(1, n+1);
surviveCensor = zeros(1, n+1);
% define for indexes
for i = 2:n+1
    dataWithCensoring2(i, 1) = i;
end

diagnosed = NaT(1, n);
death = NaT(1, n);
% loop for setting diagnosis and death dates
for i = 1:n
    % use same diagnosed date as before
    dataCensor(i+1, 1) = data(i+1, 1);
    
    % use same death date as before
    dataCensor(i+1, 2) =  data(i+1, 2);
end

% find amount of time between events
for i = 2:n+1
    dataWithCensoring2(i, 2) = data2(i, 2);
end

% randomly show if death or censored using normal distribution
for i = 2:n+1
    r = normrnd(0.6, 0.1);
    if r >= 0.5
        dataWithCensoring2(i, 3) = 1;
    else
        dataWithCensoring2(i, 3) = 0;
    end
end

% sort the data based on time between events
sortedDataCensor = sortrows(dataWithCensoring2, 2);

% sort the event dates relative to time between events
sortedDatesCensored = NaT(n+1, 2);
for i = 2:n+1
   x = sortedDataCensor(i, 1);
   sortedDatesCensored(i, 1) = dataCensor(x, 1);
   sortedDatesCensored(i, 2) = dataCensor(x, 2);
end

% set at risk to be everyone, and P(survive) = 1 to start
atRiskCensor(1) = n;
surviveCensor(1) = 1;

% use Kaplan-Meier to estimate survival probability
for i = 2:n
    j = i;
    numDeath = 0;
    % if multiple events happened in the same amount of time
    if sortedDataCensor(j, 2) == sortedDataCensor(j+1, 2)
        % count how many events occured in same time
        while sortedDataCensor(j, 2) == sortedDataCensor(j+1, 2)
            % if there was a death
            if sortedDataCensor(j, 3) == 1
                numDeath = numDeath + 1;
                atRiskCensor(j) = atRiskCensor(j - 1) - 1;
            % if there was censored data
            else
                atRiskCensor(j) = atRiskCensor(j - 1) - 1;
            end
            if j <= 1000
                j = j + 1;
            end
        end
        surviveCensor(i) = surviveCensor(i - 1) * (1 - (numDeath / atRiskCensor(i)));
    else
        % if there was a death
        if sortedDataCensor(i, 3) == 1
            numDeath = numDeath + 1;
            atRiskCensor(i) = atRiskCensor(i - 1) - 1;
            surviveCensor(i) = surviveCensor(i - 1) * (1 - (numDeath / atRiskCensor(i)));
        % if there was censored data
        else
            atRiskCensor(i) = atRiskCensor(i - 1) - 1;
            surviveCensor(i) = surviveCensor(i - 1) * 1;
        end
    end
    i = j;
    numDeath = 0;
end

% Estimate hazard function with Nelson-Aalen Estimator
hazardCensor = zeros(1, n+1);
atRiskCensor(1) = n;
for i = 2:n    
    j = i;
    numDeath = 0;
    % if multiple events happened in the same amount of time
    if sortedDataCensor(j, 2) == sortedDataCensor(j+1, 2)
        % count how many events occured in same time
        while sortedDataCensor(j, 2) == sortedDataCensor(j+1, 2)
            % if there was a death
            if sortedDataCensor(j, 3) == 1
                numDeath = numDeath + 1;
                atRiskCensor(j) = atRiskCensor(j - 1) - 1;
            % if there was censored data
            else
                atRiskCensor(j) = atRiskCensor(j - 1) - 1;
            end
            if j <= n
                j = j + 1;
            end
        end
        % estimate hazard probability
        hazardCensor(i) = hazardCensor(i-1) + (numDeath / atRiskCensor(i));
    else
        % check for censoring
        if sortedDataCensor(i, 3) == 1
            numDeath = numDeath + 1;
            atRiskCensor(i) = atRiskCensor(i - 1) - 1;            
        else
            atRiskCensor(i) = atRiskCensor(i - 1) - 1;
        end
        % estimate hazard probability
        hazardCensor(i) = hazardCensor(i-1) + (numDeath / atRiskCensor(i));
    end
    i = j;
    numDeath = 0;
    % check for the last element
    if i == n
        % if there was a death
        if sortedDataCensor(i + 1, 3) == 1
            numDeath = numDeath + 1;
            atRiskCensor(i + 1) = atRiskCensor(i) - 1;            
        % if there was censored data
        else
            atRiskCensor(i + 1) = atRiskCensor(i) - 1;
        end
        hazardCensor(i + 1) = hazardCensor(i) + (numDeath / atRiskCensor(i + 1));
    end
end

% Plot survival Function
figure(2)
p1Censor = plot(sortedDataCensor(:, 2), surviveCensor);
hold on
p1 = plot(sortedData(:, 2), survive);
hold off
title({'Estimated Survival Curve $\hat{S}$(t)' 'Censored vs Non-Censored'}, 'Interpreter','latex');
ylabel('P(t)');
xlabel('Time (days)');
legA = [p1 p1Censor];
legend(legA, '$\hat{S}$(t)', '$\hat{S}$(t) Censored', 'Interpreter','latex');

% Plot the hazard function
figure(3)
p2Censor = plot(sortedDataCensor(:, 2), hazardCensor);
hold on
p2 = plot(sortedData(:, 2), hazard);
hold off
title({'Estimated Hazard Curve $\hat{H}$(t)' 'Censored vs Non-Censored'}, 'Interpreter','latex');
ylabel('Hazard Rate');
xlabel('Time (days)');
legA = [p2 p2Censor];
legend(legA, '$\hat{H}$(t)', '$\hat{H}$(t) Censored', 'Interpreter','latex');

% Estimate a survival curve from hazard function
sFromHCensor = zeros(1, n+1);
for i = 1:n+1
    sFromHCensor(i) = exp(-1 * hazardCensor(i));
end

% compare two methods of finding survival curve
figure(4)
hold on
pl4Censor = plot(sortedDataCensor(:, 2), sFromHCensor);
pl3Censor = plot(sortedDataCensor(:, 2), surviveCensor);
title('$\hat{S}$(t) from Kaplan-Meier vs $\hat{S}$(t) from Nelson-Aalen', 'Interpreter','latex');
ylabel('P(t)');
xlabel('Time (days)');
h = [pl4Censor pl3Censor];
legend(h, '$\hat{S}$(t) from Kaplan-Meier', '$\hat{S}$(t) from Nelson-Aalen', 'Interpreter','latex');
hold off

% calculate percent difference between two censored survival curves
percDiffKap = zeros(1, n);
for i = 1:n-1
    percDiffKap(i) = ((abs(sFromHCensor(i) - surviveCensor(i))) / ((sFromHCensor(i) + surviveCensor(i)) / 2)) * 100;
    %percDiff(i) = (abs(sFromH(i) - survive(i)) / survive(i)) * 100;
end

% plot percent difference between two censored survival curves
figure(5)
bar(percDiffKap, 'FaceColor', [0.4940 0.1840 0.5560]);
title('Percent Difference Between Kaplan-Meier and Nelson-Aalen $\hat{S}$(t), Censored', 'Interpreter','latex');
xlabel('Number of Event Occurences');
ylabel('Percent Difference');

% plot all of the survival curves, censored and not censored
figure(6)
hold on
new1 = plot(sortedData(:, 2), sFromH);
new2 = plot(sortedData(:, 2), survive);
new3 = plot(sortedDataCensor(:, 2), sFromHCensor);
new4 = plot(sortedDataCensor(:, 2), surviveCensor);
title('$\hat{S}$(t) from Kaplan-Meier vs $\hat{S}$(t) from Nelson-Aalen', 'Interpreter','latex');
ylabel('P(t)');
xlabel('Time (days)');
legNew = [new1 new2 new3 new4];
legend(legNew,'$\hat{S}$(t) from Kaplan-Meier', '$\hat{S}$(t) from Nelson-Aalen','$\hat{S}$(t) from Kaplan-Meier Censored', '$\hat{S}$(t) from Nelson-Aalen Censored', 'Interpreter','latex');
hold off


% calculate percent difference Nelson
percDiffNel = zeros(1, n);
for i = 1:n-1
    percDiffNel(i) = ((abs(sFromH(i) - sFromHCensor(i))) / sFromHCensor(i)) * 100;
    %percDiff(i) = (abs(sFromH(i) - survive(i)) / survive(i)) * 100;
end

% calculate percent difference Kaplan
percDiffKap = zeros(1, n);
for i = 1:n-1
    percDiffKap(i) = ((abs(survive(i) - surviveCensor(i))) / surviveCensor(i)) * 100;
    %percDiff(i) = (abs(sFromH(i) - survive(i)) / survive(i)) * 100;
end

% plot the percent difference between the censored/non-censored survival
% curves
figure(8)
% vals = [percDiffKap; percDiffNel];
% bar((1:100), vals);
% legend('Kaplan-Meier', 'Nelson-Aalen', 'Location','best');
hold on
bar1 = bar(percDiffKap);
bar2 = bar(percDiffNel);
title({'Percent Difference Between Kaplan-Meier and Nelson-Aalen $\hat{S}$(t),' 'Censored vs Non-Censored'}, 'Interpreter','latex');
xlabel('Number of Event Occurences');
ylabel('Percent Difference');
legBar = [bar1 bar2];
legend(legBar,'Kaplan-Meier', 'Nelson-Aalen', 'Interpreter','latex', 'Location','best');
hold off