function tas = two_atom_system(n,k,nm,seed1,seed2, threshold, temp, y)
% function makes two rbatoms with n k and nm (numorbs) seeded by 1 and 2,
% and if the atoms stick for >threshold iterations at any point in time,
% tas returns 1



total_events = 100;     % number of vecworld/implementer loops to do (each loop corresponds to one event (or sometimes zero events))

% set variables for statistics gathered
natoms = zeros(1,total_events);     % each element is the number of atoms at that timestep
nmols = zeros(1,total_events);
numsb = 0;
nummb = 0;
numub = 0;

% generate atoms (currently generate a load of atoms with the same n,k, numorbs but
% different seeds)

objects{2} = [];       % pre-allocate object cell array

objects{2} = rbatom(n, k, nm, seed2);    
objects{1} = rbatom(n, k, nm, seed1);


% start vecworld/implementer loop
for iteration = 1:total_events
   % vectorize objects
   v = zeros(1,length(objects));
   for n = 1:length(objects)
       v(n) = vectorizer(objects{n});
   end
   objects = objects(v==1);
   if length(objects)<2
        tas = -1;
        return
   end
   
   %disp(['Using objects ',num2str(v==1),' on iteration ',num2str(iteration)]);
    
    
   [con_vec, mode, which_affected] = vecworld(objects, temp, y);
   not_affected = setdiff([1:length(objects)], which_affected);
   if ~isempty(con_vec) && ~isempty(mode)
       %disp(['Event ',num2str(mode),' occurs between ',num2str(con_vec{1,1}),' and ',num2str(con_vec{2,1}),' on iteration ',num2str(iteration)])
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


% look through natoms, if 0 for >threshold iterations, return 1
tick =zeros(1,total_events);
count = 1;
for t = 2:total_events
    if natoms(t) == natoms(t-1) && natoms(t)==0
        tick(count) = tick(count)+1;
    else
        count = count+1;
    end
end

max_stuck_time = max(tick);
if max_stuck_time >=threshold
    tas = 1;
else
    tas = 0;
end
disp('End of run')
end