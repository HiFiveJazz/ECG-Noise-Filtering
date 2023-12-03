load('data/train_data.mat')
ecg=train_data{1}.ecg;
fs=train_data{1}.fs;
ecg1 = ecg(:,1);
ecg2 = ecg(:,2);
ecg3 = ecg(:,3);
x = linspace(0, 152800, 152801);
%plot(x, ecg1)

average_ecg = mean(ecg, 2);
plot(x, average_ecg)
hold on
plot(x, ecg1)
plot(x, ecg2)
plot(x, ecg3)
legend("avg", "1", "2", "3");

