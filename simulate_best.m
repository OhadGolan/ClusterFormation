%%%Simulation for calculating the number of white and colored cells in the
%%%cluster of the best cluster life time and recruitment time

function [sorted_final_cell]=simulate_best(i_best, t_final_best)
dt=1;                    %the time step between each interval in hours
% t_division=6;                         %average division  time

%cell cycle times
Time_start_divide=10; %the probability to divide starts increasing at 10 h - after that the probability to divide increases
Time_end_divide=13;   %if cells haven't divided until this point, the probability to divide is 1

precent_of_colored_cells_tissue=0.25;   %precent of colored cells in the tissue
precent_of_white_cells_tissue=0.75;     %precent of white cells in tissue


 T_final=t_final_best;
 i=i_best;
 
        p_recruitment= 1/i;
        p_recruitment_colored=p_recruitment*precent_of_colored_cells_tissue;        %the probability
        p_recruitment_white=p_recruitment*precent_of_white_cells_tissue;            %the probability
        
        cells_final=zeros(3000,3);               %declare the final cells count array
        sorted_final_cell=zeros(3000,3);         %declare the final sorted cells count array
        final_cell_array=zeros(3000,2);         %declare the 
        
        for run=1:90
            cell_array=zeros(100, 4); % the array of cells. first collum if there is a cell: 1 there is a cell, 0 means there is no cell. second collum shows if the cell is white or colored: white is 0, colored is 1. third colum shows the time of last devision of the cell.
            % fourth collum: probability of division. The probability of division goes according to a distribution that accumulates over time:
            %probability of division= (t-t_last_division)/((t-t_last_division) +
            %k_division) . This means that the probability to divide increases
            %constantly over time until some point where the division is almost
            %certain.
            
            start_of_cluster=0;  %a marker that shows if the cluster has started already or not. 0 no cluster yet, 1 there is already a cluster
            
            for t=0:dt:T_final                      %time step over all the processess from 0 to the final time
                
                %The first recruitment of the initial cells in the cluster happens
                %randomly through the whole time of the simulation
                if start_of_cluster==0                   % if the cluster hasn't started yet, check if the cluster is to start at this time point
                    start_of_cluster_rand=rand();        % generate a random number to see if the cluster is to be started
                    if start_of_cluster_rand<(1/T_final)  %if the random number for the start of the cluster if lower than the probability for the start of the cluster - start the cluster
                        cell_array(1,1)=1;                      %start the cluster with two cells
                        cell_array(2,1)=1;
                        rand_initial=rand();                    %generate a random variable to select what color the first two cells of the cluster
                        if rand_initial<precent_of_colored_cells_tissue         %set the first  two cells to be colored or white. if the randon number is lower than the precent of colored cells in the cluster the first two cells are colored
                            cell_array(1,2)=1;
                            cell_array(2,2)=1;
                        else
                            cell_array(1,2)=0;
                            cell_array(2,2)=0;
                        end
                        cell_array(1,3)=t;
                        cell_array(2,3)=t;
                        start_of_cluster=1;
                        
                    end
                end
                
                A=sum(cell_array);
                number_of_cells = A(1,1);           %find the total number of cells in the cluster.
                
                if start_of_cluster==1      %if the cluster has already started continue with proliferation and recruitment
                    
                    if number_of_cells >0
                        cell_array(:,4) = cell_array(:,1).* DivisionProb(Time_start_divide, Time_end_divide, t, cell_array(:,3));  %calculate the probability of division for each cell. the first term makes only the possible cells a probability to divide
                        
                        random_division_num=rand(number_of_cells, 1);                                                     % generate a random number for each cell
                        
                        cells_that_divide=find(ceil(cell_array(1:number_of_cells,4) - random_division_num));        %find the indices of which cells divide. all the cell with probability to decide lower than the random numnber will not divide.
                        cell_array([cells_that_divide], 3)= t;                                      % set the division time of the dividing cells to t
                        number_of_colored_cells_dividing = sum(cell_array([cells_that_divide],2));  %calculate the number of colored cells that divide
                        
                        cell_array(number_of_cells+1:(number_of_cells + length(cells_that_divide)), 1)=1;              %add the number of cells that divided to the cell array
                        cell_array(number_of_cells+1:(number_of_cells + number_of_colored_cells_dividing), 2)= 1;      % allocate a color to all the new cells -
                        cell_array(number_of_cells+1:(number_of_cells + length(cells_that_divide)), 3)=t;              % allocate the time of division to all the new cells
                    end
                    
                    A=sum(cell_array);
                    number_of_cells = A(1,1);          %find the total number of cells in the cluster again
                    
                    random_recruit_num=rand(2,1);      %generate two random number for the possible recruitment of a colored or a white cell
                    if random_recruit_num(1,1)<p_recruitment_white      %if the first random number is lower than the probability to recruit a white cell add a white cell with last division time of t
                        cell_array(number_of_cells + 1, 1)=1;
                        cell_array(number_of_cells + 1, 2)=0;
                        cell_array(number_of_cells + 1, 3)=t;
                    end
                    
                    A=sum(cell_array);
                    number_of_cells = A(1,1);                               %find the total number of cells in the cluster again
                    if random_recruit_num(1,1)> (1-p_recruitment_colored)        %if the first random number is lower than the probability to recruit a colored cell add a white cell with last division time of t
                        cell_array(number_of_cells + 1, 1)=1;
                        cell_array(number_of_cells + 1, 2)=1;
                        cell_array(number_of_cells + 1, 3)=t;
                    end
                    A=sum(cell_array);
                    number_of_cells = A(1,1);
                end
            end
            number_of_colored_cells=sum(cell_array(:,2)); %calculae rhe number of colored cells
            final_cell_array(run,1)=number_of_cells;
            final_cell_array(run,2)=number_of_colored_cells;
        end
        cells_final(:,1)=final_cell_array(:,1); % the total number of cells in each experiment
        cells_final(:,2)=final_cell_array(:,1)-final_cell_array(:,2); % the number of white cells in each experiment
        cells_final(:,3)=final_cell_array(:,2); % the number of colored cells in each experiment
        [sorted_final_cell, I]= sort(cells_final(:,1)); %sort the experiments by the number of cells
        sorted_final_cell(:,2)=cells_final(I,2);        %fit the number of white cells accordingly
        sorted_final_cell(:,3)=cells_final(I,3);        %fit the number of colored cells accordingly
        sorted_final_cell(find(sorted_final_cell(:,3)==0), :)=[];   %remove the clusters with no colored cells
           

%plot figure
figure
bar(sorted_final_cell(:,2:3), 'stacked');
title('Simulation of clusters');
xlabel('Cluster');
ylabel('Number fo cells');
axis([-inf inf 0 20])
legend('cKit^+ cells', 'CKit^- cells');

end