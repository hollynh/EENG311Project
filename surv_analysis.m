% EENG 311- Final Project
% Group 3
% December 15, 2021
% This code analyzes how different models of censoring affect the estimated
% survival function from the Kaplan-Meier estimator
function [expec,std_dev,variance] = surv_analysis(varargin)
Defaults = {100};
Defaults(1:nargin) = varargin;

expec = zeros(10,1);
std_dev = zeros(10,1);
variance = zeros(10,1);

load("baseline.mat");
load("surv_data.mat");
hold on;

plot(base(:,1),base(:,4),'LineWidth',2);

y_raw = sort(data(:,4));

expec(1) = mean(y_raw);
std_dev(1) = std(y_raw);
variance(1) = var(y_raw);

y = unique(y_raw);
freq = zeros(length(y),1);
for i = 1:length(y)
    freq(i,1) = length(find(y_raw == y(i)));
end

%figure();
ecdf(y,'Frequency',freq,'Function','survivor');
title('No censoring (raw data)');
legend('true survivorship','KM - no censoring');

l_y = length(y_raw);
cens = zeros(l_y/100,3);

uni = (1:100:l_y)';
cens(:,1) = y_raw(uni);

r_norm = sort(randi(l_y,l_y/100,1));
cens(:,2) = y_raw(r_norm);

r_uni = sort(unidrnd(l_y,l_y/100,1));
cens(:,3) = y_raw(r_uni);


avg_decens = zeros(l_y/100,3);
clos_decens = zeros(l_y/100,3);

avg = zeros(l_y/100,6);

avg(:,1) = uni - 50;
avg(:,2) = uni + 50;

avg(:,3) = r_norm - 50;
avg(:,4) = r_norm + 50;

avg(:,5) = r_uni - 50;
avg(:,6) = r_uni + 50;

avg(avg<1)=1;
avg(avg>l_y)=l_y;

odds = 1:2:5;
for i = 1:length(avg)
    for j = 1:3
        avg_decens(i,j) = mean(y_raw(avg(i,odds(j)):avg(i,odds(j)+1)));
        [~,closest] = min(abs(y_raw(avg(i,odds(j)):avg(i,odds(j)+1)) - avg_decens(i,j)));
        clos_decens(i,j) = y_raw(avg(i,odds(j)) - 1 +  closest);
    end
end
clear odds;

expec(2) = mean(cens(:,1));
std_dev(2) = std(cens(:,1));
variance(2) = var(cens(:,1));
expec(3) = mean(cens(:,2));
std_dev(3) = std(cens(:,2));
variance(3) = var(cens(:,2));
expec(4) = mean(cens(:,3));
std_dev(4) = std(cens(:,3));
variance(4) = var(cens(:,3));
expec(5) = mean(avg_decens(:,1));
std_dev(5) = std(avg_decens(:,1));
variance(5) = var(avg_decens(:,1));
expec(6) = mean(avg_decens(:,2));
std_dev(6) = std(avg_decens(:,2));
variance(6) = var(avg_decens(:,2));
expec(7) = mean(avg_decens(:,3));
std_dev(7) = std(avg_decens(:,3));
variance(7) = var(avg_decens(:,3));
expec(8) = mean(clos_decens(:,1));
std_dev(8) = std(clos_decens(:,1));
variance(8) = var(clos_decens(:,1));
expec(9) = mean(clos_decens(:,2));
std_dev(9) = std(clos_decens(:,2));
variance(9) = var(clos_decens(:,2));
expec(10) = mean(clos_decens(:,3));
std_dev(10) = std(clos_decens(:,3));
variance(10) = var(clos_decens(:,3));

for i = 1:3
    y_cen = unique(cens(:,i));
    y_avg = unique(avg_decens(:,i));
    y_clos = unique(clos_decens(:,i));

    freq_cen = zeros(length(y_cen),1);
    freq_avg = zeros(length(y_avg),1);
    freq_clos = zeros(length(y_clos),1);

    for j = 1:length(y_cen)
        freq_cen(j,1) = length(find(cens(:,i) == y_cen(j)));
    end
    
    for j = 1:length(y_avg)
        freq_avg(j,1) = length(find(avg_decens(:,i) == y_avg(j)));
    end
    
    
    for j = 1:length(y_clos)
        freq_clos(j,1) = length(find(clos_decens(:,i) == y_clos(j)));
    end

    figure();
    hold on;
    plot(base(:,1),base(:,4),'LineWidth',2);
    ecdf(y_cen,'Frequency',freq_cen,'Function','survivor');
    ecdf(y_avg,'Frequency',freq_avg,'Function','survivor');
    ecdf(y_clos,'Frequency',freq_clos,'Function','survivor');

    switch i
        case 1
            title('Uniform censoring');
            legend('true survivorship','KM - uniform censoring','KM - average decensoring','KM - closest decensoring');
        case 2
            title('Normal censoring');
            legend('true survivorship','KM - normal censoring','KM - average decensoring','KM - closest decensoring');
        case 3
            title('Uniform random censoring');
            legend('true survivorship','KM - uniform random censoring','KM - average decensoring','KM - closest decensoring');
    end
end