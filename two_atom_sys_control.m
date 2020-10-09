% two_atom_sys_control is a short script to do looping through seeds which
% are input into the two_atom_system function, and to output the useful
% results to a file

min_seed_number = 30;
max_seed_number = 100;
n = 40;             % number of nodes in the rbatom networks
k = 2;              % number of inputs to each node
nm = 4;             % number of vectors produced from each rbatom network
threshold = 8;     % number of iterations atoms are required to stick for to log the seeds
master_seed = 130;  % one seed to bind them (used for atom generation)
temp = 0.5;         % value between 0 and 1, affects rotation speed
y = 3;              % number of largest vectors to check for object-object bonding (<=numorbs)
rng(master_seed);
results = [];       % [n, k, nm, seed1, seed2, threshold, temp, y]
filedump = 1000;    % number of data points to dump in file at a time
output_file_name = 'two_atom_data_0609';
fileID = fopen(output_file_name,'w');
fprintf(fileID,'[n, k, nm, seed1, seed2, threshold, temp, y]');
fclose(fileID);

for seed1 = min_seed_number:max_seed_number
    for seed2 = min_seed_number:max_seed_number
        if two_atom_system(n, k, nm, seed1, seed2, threshold, temp, y) == 1
            data_vector = [round(n),round(k),roudn(nm),round(seed1),round(seed2),round(threshold),temp,round(y)];
            results = vertcat(results,data_vector);
        end
        [num_data, ~] = size(results);
        
        %{
        if num_data>= filedump
            fileID = fopen(output_file_name,'w');
            %dlmwrite(output_file_name,results,'-append');
            fprintf(fileID,'%\t',results);
            fclose(fileID);
            results = [];
        end
        %}
    end
end
disp('End of two_atom_sys_control')
