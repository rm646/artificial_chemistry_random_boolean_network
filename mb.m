function molecule_bond = mb(species1, species2, y)
% function takes two species and (excluding any bonded orbs) uses the
% largest y (in magnitude) on each and tests all possible combinations of those before
% finally selecting one of those possibilities randomly as the one to bond, returning the connection
% vector.

orbs1 = species1.orbvec;
nmorbs1 = species1.numorbs;
otols1 = species1.orbtol;
bonded1 = species1.bonded;    % 1 = bonded, 0 = not bonded

orbs2 = species2.orbvec;
nmorbs2 = species2.numorbs;
otols2 = species2.orbtol;
bonded2 = species2.bonded;

% extract magnitudes of orbital vectors
mags1 = zeros(1,nmorbs1);
mags2 = zeros(1,nmorbs2);
for j = 1:nmorbs1
    mags1(j) = norm(orbs1(:,j));
end
for j = 1:nmorbs2
    mags2(j) = norm(orbs2(:,j));
end

% sort the vectors in descending order of magnitude
[~,ind1] = sort(mags1, 'descend');
[~,ind2] = sort(mags2, 'descend');
%sorted_orbs1 = orbs1(:,ind1);
%sorted_orbs2 = orbs2(:,ind2);
sorted_bonded1 = bonded1(ind1);
sorted_bonded2 = bonded2(ind2);

% select the first y of the sorted vectors that are not bonded and track
% which index they had in the original orbs1 / orbs2 matrices
counter1 = 0;
for k=1:nmorbs1 
    if sorted_bonded1(k) ==0 && counter1 < y
    counter1 = counter1 + 1;    
    %selected_orbs1(:,counter1) = sorted_orbs1(:,k);
    lookup_ind1(counter1) = ind1(k);
    end
end %for

counter2 = 0;
for k=1:nmorbs2 
    if sorted_bonded2(k) ==0 && counter2 < y
    counter2 = counter2 + 1;
    %selected_orbs2(:,counter2) = sorted_orbs2(:,k);
    lookup_ind2(counter2) = ind2(k);
    end
end %for

% exhaustively test the combinations of these <=y largest orbitals for bonding
bondcount = 0;
for j1 = 1:counter1
    for j2 = 1:counter2
        isbond = ob(orbs1(:,lookup_ind1(j1)),orbs2(:,lookup_ind2(j2)),otols1(lookup_ind1(j1)),otols2(lookup_ind2(j2)));
        if isbond ==1
            bondcount = bondcount+1;
            O(1,bondcount) = lookup_ind1(j1);
            O(2,bondcount) = lookup_ind2(j2);
        end
    end
end

% randomly select one of the possibilities (if there is one)
if bondcount >0
    O = O(:,randi(bondcount));
    % format as usual connection_vector
    con_vec{2,1} = [2,O(2,1)];
    con_vec{1,1} = [1,O(1,1)];
    % return the connection
    molecule_bond = con_vec;
    
else
    molecule_bond = [];
end

end