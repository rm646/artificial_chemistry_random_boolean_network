function dc = disconnect(object, connection_vector)
% function takes an object and a connection_vector to remove as inputs,
% then effects the disconnection, returning the resultant objects.

% if object is an atom, simply rebuild it with an adjusted connection
% matrix

if object.isatm == 1
    con_mat = object.selfcon;
    [~,numcons] = size(con_mat);
    di = zeros(1,numcons);
    for j = 1:numcons
        di(j) = isequal(connection_vector, con_mat(:,j));
    end
    ind = di == 0;
    con_mat = con_mat(:,ind);   % remove the offending connection vector
    bonded = object.bonded;     % update bonding vector
    bonded(di==1) = 0;
    object.bonded = bonded;
    object.selfcon = [];        % clear history of self-connections
    [~,numcons] = size(con_mat);
    if numcons>0                % so long as no connections left
    for j = 1:numcons
        object = self_connect(object, con_mat(:,j));     % perform self_connections
    end
    end
    dc{1} = object;     % return resulting atom
end %ifatom

% if not, need to sort things out
if object.isatm == 0
    % retrieve con_matrix
    con_mat = object.selfcon;
    % update it (drop the removed connection)
    [~,numcons] = size(con_mat);
    di = zeros(1,numcons);
    for j = 1:numcons
        di(j) = isequal(connection_vector, con_mat(:,j));
    end
    ind = di == 0;
    con_mat = con_mat(:,ind);   % remove the connection vector
    
    % create another list of connections which has all atom-self
    % connections removed to allow easier sorting
    [~,numcons] = size(con_mat);
    count = 0;
    included = zeros(1,numcons);    % vector holding the column indices of the con_matrix to be included (atom-self connections excluded)
    % extract the atomic indices connection matrix (ignore orb indices)
    % leaving out the connection when atom-self
    for j = 1:numcons
        ele1 = con_mat{1,j};
        ele2 = con_mat{2,j};
        if ele1(1,1) ~= ele2(1,1)
            count = count + 1;
            included(count) = j;
        end
    end
    included = included(included~=0);   % clear any unset elements
    
    % create the selfless connection matrix
    selfless_con = con_mat(:,included);
    
    % sort the new list of connections, see how many groups there are
    [~,snumcons] = size(selfless_con);
    if snumcons>0
        [~,atommat] = molsort(selfless_con);        % ith row of atommat is a list of the atoms in the ith group
        [numgroups, maxmems] = size(atommat);
        if numgroups ==1
            archive = object.archive;
            dc{1} = rbmol(archive, con_mat);
        elseif numgroups ==2
            % sort the connection matrices
            % first find indices of atoms in each
            alist1 = atommat(1,:);
            alist1 = alist1(alist1~=0);         % remove zeroes
            alist2 = atommat(2,:);
            alist2 = alist2(alist2~=0);
            
            
            % create a list of atoms from con_mat to simplify searching
            A = zeros(size(con_mat));
            for j = 1:numcons
                for i = 1:2
                    element = con_mat{i,j};
                    A(i,j) = element(1,1);
                end
            end

            % exhaustively search this A matrix for the atoms in each alist
            % to provide indices for the columns relevant to each group
            grp_ind1 = zeros(1,numcons);
            grp_ind2 = zeros(1,numcons);
            for j = 1:numcons
               for i = 1:2
                    grp_ind1(j) = grp_ind1(j) + max(A(i,j) == alist1);          % does either end of this connection belong anywhere in alist1, if so, +1 to grp ind
                    grp_ind2(j) = grp_ind2(j) + max(A(i,j) == alist2);
               end
            end
            con1 = con_mat(:,grp_ind1>0);
            con2 = con_mat(:,grp_ind2>0);
            allatomsinvolved = object.archive;
            atoms1 = allatomsinvolved(alist1);   % select those atoms in the list
            atoms2 = allatomsinvolved(alist2);

            % depending on whether one of the groups is a single atom or not,
            % create rbmols/atoms/mix and return the results in a cell array
            if length(alist1) ==1
                if length(alist2)==1
                    dc{1} = atoms1{1};
                    dc{2} = atoms2{1};
                else
                    dc{1} = atoms1{1};
                    dc{2} = rbmol(atoms2, con2);
                end
            elseif length(alist1) >1
                if length(alist2)==1
                    dc{1} = rbmol(atoms1,con1);
                    dc{2} = atoms2{1};
                else
                    dc{1} = rbmol(atoms1,con1);
                    dc{2} = rbmol(atoms2, con2);
                end
            end
        
            
        
        
        end %numgrps
    else    %if no connections left
        allatomsinvolved = object.archive;  %ASSUME IF 0 CONNECTIONS, THEN TWO ATOMS MUST BE THE RESULT
        dc{1} = allatomsinvolved{1};
        dc{2} = allatomsinvolved{2};
    end
end %ifmol
end %function
