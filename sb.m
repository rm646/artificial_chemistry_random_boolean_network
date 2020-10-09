function self_bond = sb(species, specind)
% function takes an atom or molecule and its index as an argument and tests its orbitals
% for self-bonding (all combinations). If a bond test is positive, the
% connection should be logged. The function returns a random connection of
% all possible connections.
orbs = species.orbvec;
nmorbs= species.numorbs;
otols = species.orbtol;
bonded = species.bonded;    % 1 = bonded, 0 = not bonded
sbcon{2,1} = [];    %initialize matrix to hold connection info

% search all pairings exhaustively
for j = 1:nmorbs
    for k = 1:nmorbs
        if j~=k && ob(orbs(:,j),orbs(:,k),otols(j),otols(k)) ==1 && (bonded(j) + bonded(k)) == 0
            con_vec{1,1} = [specind,j];
            con_vec{2,1} = [specind,k];
            sbcon = horzcat(sbcon,con_vec);
        end
    end
end

trimmed = sbcon(:,2:end);
[~,nsb] = size(trimmed);

if nsb>0
    if species.isatm ==1    % all's well if it is an atom
    self_bond = trimmed(:,randi(nsb));
    else
    selected_self_bond = trimmed(:,randi(nsb));
    % lookup the translation from molecular index to AO
    lookup = species.E;     % E = [MOind, Aind, AOind]
    element1 = selected_self_bond{1,1};
    element2 = selected_self_bond{2,1};
    MOind1 = element1(1,2);
    MOind2 = element2(1,2);
    [row1,~] = find(lookup(:,1)==MOind1);
    [row2,~] = find(lookup(:,1)==MOind2);
    Aind1 = lookup(row1,2);
    Aind2 = lookup(row2,2);
    AOind1 = lookup(row1,3);
    AOind2 = lookup(row2,3);
    
    % update the elements of the self_connection
    element1 = [Aind1,AOind1];
    element2 = [Aind2,AOind2];
    
    % return the tranformed connection
    self_bond{2,1} = element2;
    self_bond{1,1} = element1;
    
    
    end
    
else
self_bond = [];
end


end % fun