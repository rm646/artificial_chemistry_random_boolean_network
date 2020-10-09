function un_bond = ub(object)
% tests all connections of a molecule and returns a connection_vector
% indicating the connection to be broken, in terms of the molecule's
% selfcon array. The testing is done on the molecules (if a molecule) con 
% array, which holds the (roughly) same info.

orbs = object.orbvec;
nmorbs= object.numorbs;
otols = object.orbtol;
bonded = object.bonded;    % 1 = bonded, 0 = not bonded
ubcon{2,1} = [];    %initialize matrix to hold connection info

if object.isatm ==0    % molecule case
    con = object.con;
    % run through con to extract orbital indices
    [~,ncons] = size(con);
    O = zeros(2,ncons);
    for j = 1:ncons
       element1 = con{1,j};
       element2 = con{2,j};
       O(1,j) = element1(1,2);
       O(2,j) = element2(1,2);
    end
    
    % using the orbital lookups, test the connections
    still_bonded = -ones(1,ncons);
    for j = 1:ncons
        still_bonded(j) = ob(orbs(:,O(1,j)),orbs(:,O(2,j)),otols(O(1,j)),otols(O(2,j)));
    end
    
    % randomly choose one of the possible breaks if they exist, return
    % empty mat if not.
    broken = still_bonded == 0;
    if ~isempty(broken)
        choice = broken(randi(length(broken)));
        selfcon = object.selfcon;
        un_bond = selfcon(:,choice);
    else
        un_bond = [];
    end
    
else    % atom case
    selfcon = object.selfcon;
    [~,ncons] = size(selfcon);
    O = zeros(2,ncons);
    for j = 1:ncons
       element1 = selfcon{1,j};
       element2 = selfcon{2,j};
       O(1,j) = element1(1,2);
       O(2,j) = element2(1,2);
    end
    
    % using the orbital lookups, test the connections
    still_bonded = -ones(1,ncons);
    for j = 1:ncons
        still_bonded(j) = ob(orbs(:,O(1,j)),orbs(:,O(2,j)),otols(O(1,j)),otols(O(2,j)));
    end
    
    % randomly choose one of the possible breaks
    broken = still_bonded == 0;
    if ~isempty(broken)
        choice = broken(randi(length(broken)));
        un_bond = selfcon(:,choice);
    else
        un_bond = [];
    end
end 


end %ub