figure 
% Example data
x = [640.4, 132.7, 217.6, 557.2, 47.4,114.9,190.6,41.4];
y = [0.97002, 0.96856, 0.96808, 0.96948, 0.92944,0.94910, 0.96324, 0.82043];
labels = {'1A', '1B', '1C', '1D','1E','1F','1G','1H'};
% Plot scatter plot
scatter(x, y, 'filled');
hold on;
% Set the offset for labels
label_offset = 10;
% Add labels for each point
for i = 1:length(x)
    text(x(i) + label_offset, y(i), labels{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
end
% Add axis labels
xlabel('Model Duration (seconds)');
ylabel('Correlation Value');
% Add a title
title('Simulation 1');


figure
x = [695.2,143.8,239.9,617.6,53.8,140.5,210.6,45.08];
y = [0.83489,0.29763,0.82426,0.57523,0.47554,0.39766,0.51021,0.61210];
labels = {'2A', '2B', '2C', '2D','2E','2F','2G','2H'};
% Plot scatter plot
scatter(x, y, 'filled');
hold on;
% Set the offset for labels
label_offset = 10;
% Add labels for each point
for i = 1:length(x)
    text(x(i) + label_offset, y(i), labels{i}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
end
% Add axis labels
xlabel('Model Duration (seconds)');
ylabel('Correlation Value');
% Add a title
title('Simulation 2');