function ca = connect_atoms(atom1,atom2, ref1, ref2)
    % takes two atoms and references to an orbital on each (as a number
    % representing an index if one connection, row vector otherwise).
    % CHANGES ATOM1 TO HOLD PROPERTIES OF THE CONNECTED ATOM
    
    % shorten variable names for clearer code later
    K1 = atom1.K;
    n1 = atom1.n;
    no1 = atom1.numorbs;
    %f1 = atom1.flashes;
    %nn1 = floor(n1/no1);


    K2 = atom2.K;
    n2 = atom2.n;
    no2 = atom2.numorbs;
    %f2 = atom2.flashes;
    %nn2 = floor(n2/no2);
    
    x = 1;      % number of connections (later will be based on position)

    % re-index K of atom2 so that indices are unique in the end product
    K2 = K2 + n1;

    % call the prioritize function to get priority lists
    ind1 = prioritize_flash(atom1);
    ind2 = prioritize_flash(atom2);
    
    % now have a priority list of indices for both (hopefully) which we can
    % use to lookup which nodes to use based on the orbitals bonding.
    
    % re-wire (yay!) by dropping the first input to the chosen node and
    % replacing it with it's bonding partner's index
    % loop is for connecting multiple nodes
    for c = 1:x
        K1(1,ind1(ref1(c), c)) = ind2(ref2(c), c);
        K2(1,ind2(ref1(c), c)) = ind1(ref2(c), c);
    end
    
    % concatenate the two K matrices to get final result, and update
    % atom1's properties. Atom2 remains untouched.
    atom1.K = horzcat(K1,K2);
    atom1.boolmat = horzcat(atom1.boolmat,atom2.boolmat);
    atom1.initialstate = horzcat(atom1.initialstate,atom2.initialstate);
    atom1.D = vertcat(atom1.D,(atom2.D+atom1.n));
    atom1.n = atom1.n + atom2.n;
    
    % update bonding vector
    bvec = horzcat(atom1.bonded, atom2.bonded);
    bvec(ref1) = 1;
    bvec(ref2+atom1.numorbs) = 1;
    atom1.bonded = bvec;
    atom1.numorbs = atom1.numorbs + atom2.numorbs;
    sg(atom1);
    ff(atom1);  % update flashes, since we use it
    ca = 1;


end %function