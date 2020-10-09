function sc = self_connect(object, connection_vector)
% function takes specified inputs, effects the connection and updates all
% properties, returning the resulting object.

% shorten variable names for clearer code
K = object.K;

if object.isatm ==1

% call prioritize function to get a priority list for each orb
ind = prioritize_flash(object);

% retrieve orb number from the connection vector
element1 = connection_vector{1,1};
element2 = connection_vector{2,1};
ref1 = element1(1,2);
ref2 = element2(1,2);

% effect the rewiring
K(1,ind(ref1, 1)) = ind(ref2, 1);
K(1,ind(ref2, 1)) = ind(ref1, 1);

% update properties
bonding_vector = object.bonded;
bonding_vector(ref1) = bonding_vector(ref1) + 1;
bonding_vector(ref2) = bonding_vector(ref2) + 1;
object.bonded = bonding_vector;
object.selfcon = horzcat(object.selfcon,connection_vector);
sg(object);
ff(object);
ft(object);

% return object
sc{1} = object;

end % if an atom

if object.isatm ==0
    % if a molecule, create a new molecule with the additional
    % self-connection
    object = rbmol(object.archive, horzcat(object.selfcon, connection_vector)); %!@!!!!!!
    % return the object
    sc{1} = object;
end % if a molecule

end % function