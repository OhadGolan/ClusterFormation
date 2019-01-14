# ClusterFormation
Notch ligand Dll4 impairs the recruitment of hemogenic cells into intra-aortic clusters and limits Hematopoietic stem cell production.

The code simulates the formation of hematopoietic stem cell clusters in the aortic cavity. This is a probabilistic markov chain model.
The code runs the simulation for various possible parameters of the system and fits the simulation to the measured data.

The parameters of the simulation are as follows:

%raw cluster data - cluster size and ratio of colored cells:
raw_data_IgG
raw_data_dll4


precent_of_colored_cells_tissue   %precent of colored cells in the tissue
precent_of_white_cells_tissue    %precent of white cells in tissue

t_dll4_end %the length of time of cluster growth with dll4 inhibitor

%cell cycle times
Time_start_divide %the probability to divide starts increasing at 10 h - after that the probability to divide increases
Time_end_divide   %if cells haven't divided until this point, the probability to divide is 1

%intialize parameters
PVAL_best            %P-value fit for IgG treated cells
PVAL_best_dll4       % P-value fit for adll4 treated cells


%the parameter field to check the cluster life time for (T_final)
T_final_start      
T_final_end

%parameter field to check recruitment time for (i)
i_start
i_end
