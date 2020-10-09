% control sets parameters for the system (e.g. temp, number of largest connections
% to check when bonding objects) as well as setting up a number of atoms.

% set parameters
temp = 0.5;             % value between 0 and 1
y = 3;                  % number of largest vectors to check for object-object bonding (<=numorbs)
starting_number = 2;   % number of atoms to begin with
master_seed = 130;      % one seed to bind them (used for atom generation)
total_events = 100;     % number of vecworld/implementer loops to do (each loop corresponds to one event (or sometimes zero events))

% set variables for statistics gathered
natoms = zeros(1,total_events);     % each element is the number of atoms at that timestep
nmols = zeros(1,total_events);
numsb = 0;
nummb = 0;
numub = 0;

% generate atoms (currently generate a load of atoms with the same n,k, numorbs but
% different seeds)
n = 20;
k = 2;
numorbs = 4;
rng(master_seed);
%seeds = randi(master_seed,1,starting_number);
seeds = [1,1];
objects{starting_number} = [];       % pre-allocate object cell array

for j = 1:starting_number
    objects{j} = rbatom(n, k, numorbs, seeds(j));    
end %set-up atoms


% start vecworld/implementer loop
for iteration = 1:total_events
   % vectorize objects
   v = zeros(1,length(objects));
   for n = 1:length(objects)
       v(n) = vectorizer(objects{n});
   end
   objects = objects(v==1);
   %disp(['Using objects ',num2str(v==1),' on iteration ',num2str(iteration)]);
    
    
   [con_vec, mode, which_affected] = vecworld(objects, temp, y);
   not_affected = setdiff([1:length(objects)], which_affected);
   if ~isempty(con_vec) && ~isempty(mode)
       disp(['Event ',num2str(mode),' occurs between ',num2str(con_vec{1,1}),' and ',num2str(con_vec{2,1}),' on iteration ',num2str(iteration)])
       new_stuff = implementer(objects(which_affected), con_vec, mode);     % update the affected objects
       objects = horzcat(objects(not_affected),new_stuff);                  % remove affected objects and stick on new
   else
       %disp(['No event on iteration ',num2str(iteration),' between objects ',num2str(which_affected)])
   end    
   
   nm = 0;
   na = 0;
   % perform a count of atoms and molecules
    for j = 1:length(objects)
        if objects{j}.isatm ==0
            nm = nm+1;
        elseif objects{j}.isatm==1
            na = na+1;
        else
            disp('Error in reading type of objects')
        end
    end
   natoms(iteration) = na;
   nmols(iteration) = nm;
   
   % log the event
   if ~isempty(mode)
   switch mode
       case 'mb'
           nummb = nummb + 1;
       case 'sb'
           numsb = numsb + 1;
       case 'ub'
           numub = numub + 1;
   end
   end
end
disp('End of run')


% display stats
figure
plot([1:total_events],natoms,'k.')
hold on
plot([1:total_events],nmols,'bx')
xlabel(gca,'time step')
ylabel(gca,'number')
ttltxt = ['Plot of number of atoms and molecules against time']; 
title(ttltxt, 'FontSize', 16, 'FontWeight', 'bold')