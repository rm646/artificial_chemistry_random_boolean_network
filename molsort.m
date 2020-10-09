function [ms, list] = molsort(A)
% A is input as a 2 x N cell array whose elements are vectors
% containing an atom index and an orbital index.
% Atom index has to run from 1 --> natoms for the program to work.
% output is a row vector, with the element in the ith column indicating
% which group the ith atom is in (group number is the number of the lowest
% atom in that group). Also outputs a list of atoms in each group, with
% each row in atomlist corresponding to a group in ms.
sizea = size(A);
N0 = sizea(1,2);
C = zeros(2,N0);     % initialize

% extract matrix of atom indices
for j = 1:N0
    for i = 1:2
        mat = cell2mat(A(i,j));
        C(i,j) = mat(1,1);       % select atom index and put in matrix
    end
end

N1 = max(max(C));
group = [1:N1];
sizec = size(C);

for i = 1:sizec(1,2)
    % iterate over all the edges
    source = C(1,i);
    target = C(2,i);
    
    for k = 1:N1
        % find all targets and replace by source
        if k == target 
            if group(k) > group(source)
                group(k) = group(source);
            else
                group(source) = group(k);
            end
        end
    end
    
    
end

% now create cell array of cell arrays
% i.e. something to hold the various molecules' connection matrices
ug = unique(group);
nsg = length(unique(group));   % gives the number of groups/molecules/subgraphs
B{nsg,1} = [];
disp(C)
list = zeros(nsg,N1);
for i = 1:nsg;
    gi = ug(i);     % group index
    members = find(gi == group);    % gives the members of the gi group by atomic index
    list(i,:) = list(i,:) + members;
    cnew{2,1} = [];      % cell array to hold new connectivity matrix
    % search C columns from left to right (assumed equiv to chronological)
    % and search each column for any member of the group specified by gi
    % any column containing a member causes the corresponding A column to be appended to a group-specific
    % connectivity matrix/array to preserve ordering.
    for j = 1:sizec(1,2)        % run through columns of C
        for p = 1:2
            if isempty(find(C(p,j) == members,1)) == 0    % if a member of the group is in the jth column of C
                yes = 1;
            end
        end
        if yes ==1
            cnew = horzcat(cnew,A(:,j));    % then append the corresponding column of A to the group's new con matrix
        end
        yes = 0;
    end

    % we should now have the (nearly) finished product of a connectivity list, so
    % assign this to an array of such (nearly) finished products. Just need
    % to chop first column
    cnew = cnew(:,2:end);
    B{i,1} = cnew;
    clear cnew
end



% return the array of finished products
ms = B;
end 

%zeros(1, N)