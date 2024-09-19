function plot_tracers(TRACER,n,choice)

% PLOT_TRACER (TRACER,N,CHOICE) Plots the tracers evolution taken from the 
% array TRACER, in N steps. CHOICE indicate how you want to plot it, 1 is
% for plotting the horizontal tracers along time and the lines that match
% equivalent points along time, 2 is for plotting the horizontal tracers
% along time, 3 is for plotting the lines that match equivalent points 
% along time and 4 is for plotting equivalent points with the same color.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 02-12-2011. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

array_size = size (TRACER);
step = floor(array_size(4)/n);
mycolor = rand (array_size(3),3);

switch (choice)
    case {1}
   for j = 1:array_size(3)
        for i = 1:step:array_size(4)
            reordered_x_tracer = sort(TRACER(1,:,j,i));
            order = zeros(1,array_size(2));
            for u = 1:array_size(2)
                order(u) = find (reordered_x_tracer(u)==TRACER(1,:,j,i));
            end
            reordered_y_tracer = TRACER(2,order,j,i);
            plot (reordered_x_tracer,reordered_y_tracer...
                ,'-o','color',mycolor(j,:))
            hold on
        end
    end

    % Ploting lines that match the same point along time.

    % Search for zeros
    counter = 0;
    for k = 1:array_size(4)
        if TRACER(2,1,1,k) < 0
            counter = counter+1;
        else
            break
        end
    end
    mov_lines = zeros(2,counter);

    for j = 1:array_size(3)
        for i = 1:array_size(2)
            for k = 1:counter
                mov_lines(1,k) = TRACER(1,i,j,k);
                mov_lines(2,k) = TRACER(2,i,j,k);
            end
            plot (mov_lines(1,:),mov_lines(2,:),'--','color',mycolor(j,:))
        end
    end
    hold off
    
    case {2}
    for j = 1:array_size(3)
        for i = 1:step:array_size(4)
            reordered_x_tracer = sort(TRACER(1,:,j,i));
            order = zeros(1,array_size(2));
            for u = 1:array_size(2)
                order(u) = find (reordered_x_tracer(u)==TRACER(1,:,j,i));
            end
            reordered_y_tracer = TRACER(2,order,j,i);
            plot (reordered_x_tracer,reordered_y_tracer...
                ,'.','color',mycolor(j,:))
            hold on
        end
    end
    hold off
    
    case {3}
    % Ploting lines that match the same point along time.

    % Search for zeros
    counter = 0;
    for k = 1:array_size(4)
        if TRACER(2,1,1,k) < 0
            counter = counter+1;
        else
            break
        end
    end
    mov_lines = zeros(2,counter);

    for j = 1:array_size(3)
        for i = 1:array_size(2)
            for k = 1:counter
                mov_lines(1,k) = TRACER(1,i,j,k);
                mov_lines(2,k) = TRACER(2,i,j,k);
            end
            plot (mov_lines(1,:),mov_lines(2,:),'-','color',mycolor(j,:))
        end
        hold on
    end
    hold off
    
    case {4}
    mycolor = rand (array_size(2),3);
    mov_lines = zeros(2,array_size(4));
    for j = 1:array_size(3)
        for i = 1:array_size(2)
            for k = 1:step:array_size(4)
                mov_lines(1,k) = TRACER(1,i,j,k);
                mov_lines(2,k) = TRACER(2,i,j,k);
            end
            plot (mov_lines(1,:),mov_lines(2,:),'o','color',mycolor(i,:))
        end
        hold on
    end
    hold off

    otherwise
        disp('Input of the 3rd input not valid. You need to enter numbers between 1 and 3.')
end
end