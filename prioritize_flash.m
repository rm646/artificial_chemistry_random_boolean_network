function pr = prioritize_flash(rbatom)
% takes atom as input and based on its flash properties the program returns a 
% numorbsxsomething list of nodes, with each row corresponding to an
% orbital.

D = rbatom.D;   % lists orbital boundaries
%n = rbatom.n;
flashes = rbatom.flashes;
no = rbatom.numorbs;
orbflashes = -ones(no,(D(end,2) - D(end,1)));       % matrix to sort flashes into, uses fact that last orbital may be larger
for i = 1:no        % runs through orbitals
    for j = D(i,1):D(i,2)       % for each, runs through the relevant nodes
        orbflashes(i,(j-D(i,1)+1)) = flashes(1,j);
    end    
end
% orbflashes now has rows (each row corresponding to an orbital) of flashes property, which we will next sort and
% keep the indices of

[sof,soi] = sort(orbflashes, 2, 'descend');

% re-index each element in soi to be consistent with its index in the atom.
for i = 2:no
    soi(i,:) = soi(i,:) + D((i-1),2); 
end    

% set indices of unset to zero
zeroi = sof == -1;
soi(zeroi) = 0;

% return the priority list
pr = soi;

end %func