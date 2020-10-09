function [con_vec, mode, which_affected] =  vecworld(object_array, temp, y)
% takes atoms and molecules in a cell array (along with temp for rotation
% scaling) and outputs an event type as string (one of 'sb', 'ub' or 'mb')
% and a connection vector (2x1 cell array)

% set parameters, initialize arrays
    iters = 100;                    % number of time steps to run for, this could be an input parameter
    c_rate = 1;                   % collisions per unit time
    vis = 1;                        % set to 1 to plot, 0 to not
    nm = 0;                         % number of molecules
    na = 0;                         % number of atoms
    y;                          % number of largest mag orbital to test
    which_affected = [];
    
    % count numbers of each object type
    for j = 1:length(object_array)
        if object_array{j}.isatm ==0
            nm = nm+1;
        elseif object_array{j}.isatm==1
            na = na+1;
        else
            disp('Error in reading type of object_array')
        end
    end
    
if ~isempty(object_array)
counter = 0;
% randomly select two species every 1/rate seconds
for c = 0:floor(1/c_rate):iters     % every 1/rate seconds
    % randomly select two species
    numthings = length(object_array);
    sel1 = randi(numthings);
    sel2 = randi(numthings);
    
    if numthings>1
    while sel1 == sel2
        sel2 = randi(numthings);
        counter = counter + 1;
    end
    end
    
    thing1 = object_array{sel1};
    thing2 = object_array{sel2};
    
    %USE ROTATION PROPERTY?
    
    %randomly rotate
    for t=(c-floor(1/c_rate)):c
            randrotate(thing1,temp);
            randrotate(thing2,temp);
    end
    
    % start testing in a pre-determined order
    
    % first the self-bonding
    % test self-bonding
    sc1 = sb(thing1, 1);
    sc2 = sb(thing2, 1);        % SECOND ARG 1 OR 2?
    if isempty(sc1)==0 || isempty(sc2)==0     % if either thing self bonds, return the connection
        % since choice of which is thing1/thing2 is random
        % we choose arb order for this part
        if isempty(sc1) == 0
            con_vec = sc1;  % return the con vector
            mode = 'sb';
            which_affected = [sel1];
            return
        end
        if isempty(sc2) == 0
            con_vec = sc2;
            mode = 'sb';
            which_affected = [sel2];
            return
        end
    end
    
    % test decomposition
    if isempty(thing1.selfcon) ==0 || isempty(thing2.selfcon) ==0
        ub1 = ub(thing1);
        ub2 = ub(thing2);
        if isempty(ub1) ==0
            con_vec = ub1;
            mode = 'ub';
            which_affected = [sel1];
            return
        end
        if isempty(ub2) ==0
            con_vec = ub1;
            mode = 'ub';
            which_affected = [sel2];
            return
        end
    end
    
    % test bonding between thing1 and thing2
    mc = mb(thing1,thing2,y);
    if isempty(mc) ==0
        con_vec = mc;
        mode = 'mb';
        which_affected = [sel1,sel2];
        return
    end
    
    % in the case that nothing happens, return empties
    con_vec = [];
    mode = [];
    which_affected = [sel1,sel2];
    
end %every 1/rate secs
    
else % if no objects passed
    con_vec = [];
    mode = [];
    which_affected = [];
end % length(objarr)>0

end %vecworld