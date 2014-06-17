%% the following is the code for the first simulation, that is one fixed source, one moving target
% close all; clear all;
% tp = [];
% targets = [23; 20; 32; 6; 34; 10; 13; 8; 41];
% length = 0:(length(targets) - 1); % length = # cells between source and target
% 
% for T = targets' % have to use transpose or it will set T equal to the vector targets, instead of iterating through each element of targets
%     close all;
%     temp = cell_flows(T);
%     tp = [tp temp];
% end
% 
% close all;
% plot(length, tp);


%% the following is the code for the second simulation, that is two fixed source, and two moving targets with path overlapping
close all; clear all;
tp1 = [];
tp2 = [];
temp = [];
targets = [71 70; 67 57; 58 8; 10 65; 14 66; 45 60; 53 21; 55 16; 54 19; 50 17];
length = 0:(length(targets) - 1); % length = # cells between source and target

for i = 1 : 10
    close all;
    temp = cell_flows(targets(i, :));
    tp1 = [tp1 temp(1)];
    tp2 = [tp2 temp(2)];
end

close all;
plot(length, tp1);
figure;
plot(length, tp2);
