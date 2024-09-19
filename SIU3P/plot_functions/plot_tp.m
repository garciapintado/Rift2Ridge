% Plot track points

interval = 1;

X = zeros(length(TRACKP),length(TRACKP{1}(1,1:interval:end)));
Y = zeros(length(TRACKP),length(TRACKP{1}(2,1:interval:end)));
for n = 1:length(TRACKP)
    X(n,:) = TRACKP{n}(1,1:interval:end);
    Y(n,:) = TRACKP{n}(2,1:interval:end);
end

plot(X/1000,Y/1000)

hold on
%plot_box
plot_vx