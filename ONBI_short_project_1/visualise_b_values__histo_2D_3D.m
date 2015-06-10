%% visualise b values - histograms

% for b = 1:10
%     b
    figure
    subplot 221
    nbins = 30;
    hold on
    histogram(trainingdata(1:(size(names,1)/2) ,b(1)), nbins)
    histogram(trainingdata((size(names,1)/2)+1 :(size(names,1)),b(1)), nbins)
    hold off
    ylabel 'frequency'
    xlabel ([num2str(b(1)),' , ',trainingdataLabels(b(1),:)])
    legend ('DETERMINE', 'MESA', 'Location', 'best')
%     pause
% end

%% visualise b values - 2D plot

subplot 222
xlabel ([num2str(b(1)),' , ',trainingdataLabels(b(1),:)])
ylabel ([num2str(b(2)),' , ',trainingdataLabels(b(2),:)])

hold on
plot(trainingdata( 1:(size(names,1)/2) ,b(1)), trainingdata( 1:(size(names,1)/2) , b(2)), 'o')
plot(trainingdata( (size(names,1)/2)+1 :(size(names,1)),b(1)), trainingdata( (size(names,1)/2)+1:size(names,1) ,b(2)), 'o')
hold off
legend ('DETERMINE', 'MESA' ,'Location' ,'best')

%% visualise b values - 3D plot

subplot 223
hold on
plot3D(trainingdata( 1:(size(names,1)/2) ,b),'o')
plot3D(trainingdata((size(names,1)/2)+1:(size(names,1)),b),'o')
hold off
legend ('DETERMINE', 'MESA', 'Location' ,'best')
xlabel ([num2str(b(1)),' , ',trainingdataLabels(b(1),:)])
ylabel ([num2str(b(2)),' , ',trainingdataLabels(b(2),:)])
zlabel ([num2str(b(3)),' , ',trainingdataLabels(b(3),:)])
view(3)
