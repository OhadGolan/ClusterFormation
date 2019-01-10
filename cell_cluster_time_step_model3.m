%%%Simulation for calculating the number of white and colored cells in the
%%%ellepthilial tissue cluster using a Markov chain

clear all
dt=1;                    %the time step between each interval in hours

%raw cluster data - cluster size and ratio of colored cells
raw_data_IgG=[2, 2, 2, 2, 2, 3, 5, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10; 1 , 1, 1, 1, 1, 1, 0.2, 0.2, 0.4, 4./6, 3./7, 1./7, 3./8, 1./8, 2./9, 3./10, 1./10]';
raw_data_dll4=[2, 2, 2, 2, 2, 3, 2, 3, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9, 10, 10, 12, 14, 16, 30, 32; 1, 1, 1, 1, 1, 2./3, 1, 1./3, 1./4, 1, 1./5, 5./5, 1./6, 2./6, 1./6, 1./6, 1./7, 1./7, 2./7, 3./7, 4./8, 1./8, 1./9, 2./9, 2./9, 1./10, 2./10, 1./12, 2./14, 6./16, 2./30, 2./32]';


precent_of_colored_cells_tissue=0.25;   %precent of colored cells in the tissue
precent_of_white_cells_tissue=0.75;     %precent of white cells in tissue

t_dll4_end=5; %the length of time of cluster growth with dll4 inhibitor

%cell cycle times
Time_start_divide=10; %the probability to divide starts increasing at 10 h - after that the probability to divide increases
Time_end_divide=13;   %if cells haven't divided until this point, the probability to divide is 1

%intialize parameters
PVAL_best=0;            %P-value fit for IgG treated cells
PVAL_best_dll4=0;       % P-value fit for adll4 treated cells

j=0;                    %initialize loop parameters
k=0;

%the parameter field to check the cluster life time for (T_final)
T_final_start=1;        
T_final_end=30;

%parameter field to check recruitment time for (i)
i_start=1;
i_end=20;

%start siulation
for T_final=T_final_start:1:T_final_end   %go over all cluster life time parameter field
    k=k+1;
    j=0;
    for i=i_start:1:i_end                  %go over all possible recruitment times
        p_recruitment= 1/i;                %recruitment probability is 1/probability of recruitment
        
        p_recruitment_colored=p_recruitment*precent_of_colored_cells_tissue;        %the probability for recruitment of a colored cell
        p_recruitment_white=p_recruitment*precent_of_white_cells_tissue;            %the probability for recruitment of a white cell
        j=j+1;
        
        cells_final=zeros(10000,3);               %declare the final cells count array
        sorted_final_cell=zeros(10000,3);         %declare the final sorted cells count array
        simulation_data_IgG=zeros(10000,3);          %declare the simulation data array IgG
        final_cell_array=zeros(10000,2);         %declare the final cell array
        
        for run=1:2000      %run the simulation 2000 times
            cell_array=zeros(100, 4); % the array of cells. first collum is if there is a cell: 1 there is a cell, 0 means there is no cell. 
            %second collum shows if the cell is white or colored: white is 0, colored is 1. 
            %third colum shows the time of last devision of the cell.
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
                        cell_array(1,3)=t;                      %set the time of cells last division as t
                        cell_array(2,3)=t;
                        start_of_cluster=1;                     %set the cluster flag to 1
                        
                    end
                end
                
                A=sum(cell_array);
                number_of_cells = A(1,1);           %find the total number of cells in the cluster.
                
                if start_of_cluster==1      %if the cluster has already started continue with proliferation and recruitment
                    
                    if number_of_cells >0
                        cell_array(:,4) = cell_array(:,1).* DivisionProb(Time_start_divide, Time_end_divide, t, cell_array(:,3));  %calculate the probability of division for each cell. the first term makes only the possible cells a probability to divide
                        
                        random_division_num=rand(number_of_cells, 1);                                                     % generate a random number for each cell
                        
                        cells_that_divide=find(ceil(cell_array(1:number_of_cells,4) - random_division_num));        %find the indices of which cells divide. all the cell with probability to devide lower than the random numnber will not divide.
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
        
        simulation_data_IgG=[sorted_final_cell(:,1), sorted_final_cell(:,3)./sorted_final_cell(:,1)]; %arrange the data to be used in the hotelling t2 test
        
        [PVAL(k,j), T_test(k,j)]=hotell2(raw_data_IgG,  simulation_data_IgG);  %ran hosteliing t-test
   
        if abs(PVAL(k,j))>abs(PVAL_best)     %checking if this is the best score, if it is set as best value
            PVAL_best=PVAL(k,j);
            i_best=i;
            t_final_best=T_final;
        end
    end
end

%plot the simulation
T_final_plot=T_final_start:1:T_final_end;
recruitment_time=i_start:1:i_end;

figure
imagesc(recruitment_time, T_final_plot, PVAL);
title('p-value ');
xlabel('Average recruitment time (h)');
ylabel('Cluster life time (h)');
zlabel('p-value');

simulate_best(i_best, t_final_best);  %simulate the best fit for the IgG treated cells




%% simulating the experiment with anti-dll4

dt=1;                    %the time step between each interval in hours
j=0;


T_final=t_final_best;   %set T_final as the best T_final from the IgG expeiment
i=i_best;               %set the recruitment time as best recruitment time found in IgG
p_recruitment= 1/i;

p_recruitment_colored=p_recruitment*precent_of_colored_cells_tissue;        %the probability
p_recruitment_white=p_recruitment*precent_of_white_cells_tissue;            %the probability

cells_final=zeros(3000,3);               %declare the final cells count array
sorted_final_cell=zeros(3000,3);         %declare the final sorted cells count array
simulation_data_dll4=zeros(3000,3);         %declare the final sorted cells count array
final_cell_array=zeros(3000,2);         %declare the final cell array

for m=1:40  %going over all the possible recruitment times after initiation of alpha-dll4  treatment
    for run=1:2000
        cell_array=zeros(100, 4); % the array of cells. first collum if there is a cell: 1 there is a cell, 0 means there is no cell. second collum shows if the cell is white or colored: white is 0, colored is 1. third colum shows the time of last devision of the cell.
        % fourth collum: probability of division. The probability of division goes according to a distribution that accumulates over time:
        %probability of division= (t-t_last_division)/((t-t_last_division) +
        %k_division) . This means that the probability to divide increases
        %constantly over time until some point where the division is almost
        %certain.
        
        start_of_cluster=0;  %a marker that shows if the cluster has started already or not. 0 no cluster yet, 1 there is already a cluster
        
        for t=0:dt:(T_final - t_dll4_end)                      %time step over all the processess from 0 to the final time minus the 5h of adll4 tretment. In the last 5 h recruitment time will change.
            
            %The first recruitment of the initial cells in the cluster happens
            %randomly through the whole time of the simulation
            if start_of_cluster==0                   % if the cluster hasn't started yet, check if the cluster is to start at this time point
                start_of_cluster_rand=rand();        % generate a random number to see if the cluster is to be started
                if start_of_cluster_rand<(1/(T_final))  %if the random number for the start of the cluster if lower than the probability for the start of the cluster - start the cluster
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
                    cell_array(1,3)=t;                                      %set the last division time of the new cells to t
                    cell_array(2,3)=t;
                    start_of_cluster=1;                                     %set the flag of cluster start to 1
                    
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
        
        
        p_recruitment= 1/m;                                                         %set the new recruitment time effected by adll4 treatment
        p_recruitment_colored=p_recruitment*precent_of_colored_cells_tissue;        %the probability
        p_recruitment_white=p_recruitment*precent_of_white_cells_tissue;            %the probability
        
        for t= (T_final - t_dll4_end) :dt: (T_final)  % run the second part of the simulation after the dll4 inhibitor has been added
            
            %The first recruitment of the initial cells in the cluster happens
            %randomly through the whole time of the simulation
            if start_of_cluster==0                   % if the cluster hasn't started yet, check if the cluster is to start at this time point
                start_of_cluster_rand=rand();        % generate a random number to see if the cluster is to be started
                if start_of_cluster_rand<(1/(t_dll4_end))  %if the random number for the start of the cluster if lower than the probability for the start of the cluster - start the cluster
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
    
    simulation_data_dll4=[sorted_final_cell(:,1), sorted_final_cell(:,3)./sorted_final_cell(:,1)]; %arrange the data to be used in the hotelling t2 test
    [PVAL_dll4(m), T_test_dll4(m)]=hotell2(raw_data_dll4,  simulation_data_dll4);                   %ran hoteeling t-test

%    
%    
    if abs(PVAL_dll4(m))>abs(PVAL_best_dll4)     %checking if this is the best score, if best set as best parameter
        PVAL_best_dll4=PVAL_dll4(m);
        m_best=m;
    end
end

%plot simulation
m=1:40;
j=1:40;

figure
plot(m, PVAL_dll4, 'linewidth', 3);
title('p-value ');
xlabel('Average recruitment time (h)');
ylabel('p-value');
i_best_dll4=m_best;
simulate_best_dll4(i_best, t_final_best,i_best_dll4);
