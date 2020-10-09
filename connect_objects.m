function co = connect_objects(object1, object2, o1, o2)
% inputs are two objects to connect and orbital references for each
% function's role is mainly re-indexing
% first copy each object to avoid unintentional edits
thing1 = object1.cp;
thing2 = object2.cp;

% extract the properties to re-index
K2 = thing2.K;
n2 = thing2.n;
no2 = thing2.numorbs;
if thing2.isatm==0
na2 = length(thing2.archive);
end

K1 = thing1.K;
n1 = thing1.n;
no1 = thing1.numorbs;
if thing1.isatm==0
na1 = length(thing1.archive);
end

%find the point referred to on objects, and create the connection vector
if thing1.isatm ==0 
    E1 = thing1.E;
    archived_atom_ind1 = E1(o1,2);   % E matrix should give the atom index
    arc_at_orb_ind1 = E1(o1,3);
    con{1,1} = [archived_atom_ind1, arc_at_orb_ind1];
else
    con{1,1} = [1,o1];    
end

if thing2.isatm ==0
    E2 = thing1.E;
    archived_atom_ind2 = E2(o1,2);   % E matrix should give the atom index
    arc_at_orb_ind2 = E2(o1,3);
    con{2,1} = [archived_atom_ind2 + na1, arc_at_orb_ind2];
else
    con{2,1} = [2,o2];
end

% create the new list of atoms from both archives
if thing1.isatm==0
    atom_array = thing1.archive;
else
    atom_array{1} = thing1;
end

if thing2.isatm==0
    atom_array = horzcat(atom_array, thing2.archive);
else
    thing2in_cell{1} = thing2;
    atom_array = horzcat(atom_array, thing2in_cell);
end

% create the new connection matrix from both old and the new vector
first_half = thing1.selfcon;

second_half = thing2.selfcon;   %needs re-indexing

final_bit = con;

%re-index second_half
[~,ncs2] = size(second_half);
for j = 1:ncs2
    for i = 1:2
        element = second_half{i,j};
        element(1,1) = element(1,1) + na1;
        second_half{i,j} = element;
    end
end

% stick it all together
con_mat = horzcat(first_half, second_half);
con_mat = horzcat(con_mat, final_bit);

% create the resultant molecule
co{1} = rbmol(atom_array, con_mat);


end % function