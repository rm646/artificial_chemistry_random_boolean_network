classdef rbatom < handle
    %stores RBN information and vectors generated
    
    properties
        % vecworld related
        orbvec;                         % 2 x numorbs real positive matrix, columns are the orbital vectors
        numorbs;                        % user set positive integer
        orbtol;                         % 1 x numorbs real positive vector, giving angular tolerance in radians for each orbvec
        bonded;                         % 1 x numorbs binary vector, 1 indicates the corresponding column in orbvec is a 'bonded' vector
        D;                              % numorbs x 2 integer matrix, detailing the divisions of nodes into orbital chunks
        isatm;                          % object property, =1 always for rbatoms and =0 always for rbmols
        rot;                            % real valued number detailing the total angle of rotation applied to the atom since creation
        selfcon;                        % 2 x number of connections cell array detailing connections between this atom's orbitals
        
        % RBN world related
        cyclemat;                       % T x n binary matrix, representing one cycle of the RBN dynamics for the given parameters
        initialstate;                   % 1 x n binary vector, representing t=0 states of nodes
        flashes;                        % 1 x n half-integer or integer vector, representing the slashes property for each node
        totals;                         % 1 x n integer vector, representing the sums of each node's values over a cycle
        boolmat;                        % binary matrix, a lookup table for node inputs
        seed;                           % positive integer seed for rng, possibly used to seed multiple times in the constructor method
        n;                              % total number of nodes
        k;                              % number of inputs to each node
        K;                              % 2 x n positive integer matrix, the jth column represents the indices of nodes to be used as inputs for the jth node
        iters;                          % total number of iterations to seek a cycle for before giving up, currently fixed but could be made to relate to k
        T;                              % length of cycle, set to iters if no cycle is found
        
        % just in case
        isskynet = 0;
        prob_sentience = 0;
        
    end % properties
    
    methods
        
        % vecworld related
        
        function dgen = dg(rbatom)
           % generates a numorbsx2 matrix of nodes indices. The ith row
           % indicates a start and stop point for the ith orbital for the
           % atom.
           no = rbatom.numorbs;
           nn = floor(rbatom.n/no);
           DD =zeros(no,2);
           DD(1,1) = 1;
           DD(1,2) = nn;
           for i= 2:no
               DD(i,1) = DD(i-1,1)+ (nn);
           end
           for i= 2:no-1
               DD(i,2) = DD(i-1,2)+ (nn);
           end
           DD(no,2) = rbatom.n;
           rbatom.D = DD;
           dgen = 1;
        end % generates a matrix to track division into orbitals
        
        function randrotate = rr(rbatom, temp)     % rotates all vectors described by orbvec by a random angle scaled by temp and stores that angle (mod 2pi) in the rot property
            dir = (-1)^randi(2);        % positive for anti-clockwise
            x = 2*pi*temp*dir*rand;     %the angle
            o = rbatom.orbvec;
            [~,max] = size(o);
            R = [cos(x),-sin(x);sin(x),cos(x)];
            % rotate each vector
            for j = 1:max
                o(:,j) = R*o(:,j);
            end
            % update rot property
            rbatom.rot = mod(rbatom.rot + x, 2*pi);
            randrotate = 1;
        end %rr
        
        function ur = updaterotate(rbatom)
            x = rbatom.rot;
            if ~isempty(x)
            o = rbatom.orbvec;
            [~,max] = size(o);
            R = [cos(x),-sin(x);sin(x),cos(x)];
            % rotate each vector
            for j = 1:max
                o(:,j) = R*o(:,j);
            end

            ur = 1;
            else
            ur = 0;
            end
        end
        
        % RBN related
        
        function boolgen = bg(rbatom)       % using seed property, generates the lookup table           
            rng(rbatom.seed)
            m = rbatom.n;
            kk = rbatom.k;      % rename
            rbatom.boolmat = dec2bin(floor(rand(1,m)*2^(2^kk)))'-48;
            boolgen = 1;
        end
        
        function istategen = ig(rbatom)       % generates a random initial state
            rng(rbatom.seed)                  % beware, maybe seeding every time you generate will define a relationship between parts you thought were independant
            rbatom.initialstate = randi([0,1],1,rbatom.n);
            istategen = 1;
        end
        
        function kgen = kg(rbatom)      % generates random connection matrix
            rng(rbatom.seed)
            m = rbatom.n;
            kk = rbatom.k;
            [~,con] = sort(rand(m,m));
            con = con(1:kk,:);
            rbatom.K = con;
            kgen = 1;
        end
        
        function stategen = sg(rbatom)      % uses initial conditions and bool matrix to set cyclemat to a matrix of one cycle
            kk = rbatom.k;
            m = rbatom.n;
            T0 = rbatom.iters;
            KK = rbatom.K;
            B = rbatom.boolmat;
            Pow = 2.^sum(triu(ones(kk),1));
            state = zeros(T0,m);             % pre-allocate matrix to store time-evolution of system
            X0 = rbatom.initialstate;
            X = X0(1,:);
            
            %find time-evolution over iters time steps
            
            for t = 1:T0
                V = (Pow * X(KK))+1;         % input value (k bits as number)
                X = diag(B(V,:));           % lookup boolean function, set X to the next state
                state(t,:) = X;             % store time step in state matrix
            end

            % search for a cycle, if found, set cyclemat to it
            [dupesout, ind] = unique(state,'rows');
            sized = size(dupesout);
            sizes = size(state);
            cyclefound = sized(1,1) < sizes(1,1);
            if cyclefound == 1      % search using tortoise and hare algorithm
                not_enough_iters = 0;
                tortoise = 1;
                hare = 2;
                while isequal(state(tortoise,:),state(hare,:)) == 0
                    tortoise = tortoise +1;
                    hare = hare + 2;
                    if hare == T0        % hare waits if tortoise gets to the end and a flag is set
                        hare = hare-2;
                        not_enough_iters = 1;
                    end
                end
                
                if not_enough_iters ==1
                disp('Cycle extraction failed')
                rbatom.cyclemat = state;
                stategen = 1;
                return
                end
                
                mu = 1;
                tortoise = 1;
                while isequal(state(tortoise,:),state(hare+1,:)) == 0       % add +1 to hare in argument because 0 index does not exist --> tortoise cannot be set to 0
                    tortoise = tortoise + 1;
                    hare = hare + 1;
                    mu = mu + 1;
                    if hare == T0        % hare waits if it gets to the end and a flag is set
                        hare = hare-1;
                        not_enough_iters = 1;
                    end
                end
                
                if not_enough_iters ==1
                disp('Cycle extraction failed')
                rbatom.cyclemat = state;
                stategen = 1;
                return
                end
                
                lam = 1;
                hare = tortoise + 1;
                while isequal(state(tortoise,:),state(hare,:)) == 0
                    hare = hare + 1;
                    lam = lam + 1;
                    if hare == T0        % hare waits if it gets to the end and a flag is set
                        hare = hare-1;
                        not_enough_iters = 1;
                    end
                end
                
                if not_enough_iters ==1
                disp('Cycle extraction failed')
                rbatom.cyclemat = state;
                stategen = 1;
                return
                end
                
                cstart = mu;
                cend = lam + mu;
                
                rbatom.cyclemat = state(cstart:(cend-1),:);
                stategen = 1;
                rbatom.T = lam;
                
            else
                disp('Cycle not found')
                rbatom.cyclemat = state;
                stategen = 1;
            end
            %}
                     
        end

        function dispcycle = dc(rbatom)      % displays the found cyclemat as an image sorted by the totals property (calls ft)
            unsorted = rbatom.cyclemat;
            ft(rbatom);
            ttls = rbatom.totals;
            [~,i] = sort(ttls);
            sorted = unsorted(:,i);
            C = imresize(sorted,5,'nearest');
            imshow(C,'Border','tight');
            ttltxt = ['RBN cycle for n = ' num2str(rbatom.n) ' and k = ' num2str(rbatom.k)]; 
            title(ttltxt, 'FontSize', 16, 'FontWeight', 'bold')
            dispcycle = 1;
        end
        
        function findflashes = ff(rbatom)     % uses cyclemat to set flashes property
            state = rbatom.cyclemat;
            [T0,m] = size(state);
            flash = zeros(1,m);       % initialize flashes vector
            
            for i = 1:m
                for t=2:T0
                    ischange = state(t,i) ~= state(t-1,i);
                    flash(1,i) = flash(1,i) + ischange;
                end
                flash(1,i) = flash(1,i)/2;      % by its definition
            end
            rbatom.flashes = flash;
            findflashes = 1;
        end %ff
        
        function findtotals = ft(rbatom)     % uses cyclemat to set totals property
            state = rbatom.cyclemat;
            ssize = size(state);
            m = ssize(1,2);
            T = ssize(1,1);
            tots = zeros(1,m);
            for t = 1:T
                tots = tots + state(t,:);
            end
            rbatom.totals = tots;
            findtotals = 1;
        end %ft
        
        % constructor method [important]
        % constructor only sets RBN related properites, the rest are left
        % to the vectorizer
        function A = rbatom(n, k, numorbs, seed, istate, K, boolmat)    % n, k and numorbs are required, the rest may be left to be randomly generated by not specifying them
            A.n = n;
            A.k = k;
            A.numorbs = numorbs;
            A.iters = 1500;
            A.isatm = 1;
            A.bonded = zeros(1,numorbs);
            switch(nargin)      % to allow optional arguments
                case 4
                    A.seed = seed;
                    ig(A);
                    kg(A);
                    bg(A);
                    
                case 5
                    A.seed = seed;
                    A.initialstate = istate;
                    kg(A);
                    bg(A);
                
                case 6
                    A.seed = seed;
                    A.initialstate = istate;
                    A.K = K;
                    bg(A);
                 
                case 7
                    A.seed = seed;
                    A.initialstate = istate;
                    A.K = K;
                    A.boolmat = boolmat;
                    
            end % switch	
            sg(A);
            ff(A);
            ft(A);
            dg(A);
        end %constructor        
        
        function cpya = cp(A)
            cpya = rbatom(A.n, A.k, A.numorbs, A.seed, A.initialstate, A.K, A.boolmat);
            vectorizer(cpya);
      
        end %copyatm
    end % methods
end % rbatom