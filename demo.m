function d = demo(s)
%Generates an RBN, lets it find an attractor, measures its properties, then
%passes the information to the vectorizer.

s;                   % random seed, for reproducibility  % 38 good %39 v stable (n=120) both 
rng(s,'twister');       %set random number generator mode
k = 2;                  % number of connections, K
n = 120;                 % number of nodes, N
T = 100;                % number of timesteps
sortmode = 'totals';   %select property to sort according to one of 'flashes', 'totals'

% the k random connections to each node n
% generate n random permutations of length n;
% select the first k items of each
[ignore,K] = sort(rand(n,n));
K = K(1:k,:);
disp(K)
% the random boolean function at each node n
B = dec2bin(floor(rand(1,n)*2^(2^k)))'-48;
disp(B)
Pow = 2.^sum(triu(ones(k),1)); % conversion function
disp(Pow)
%================================================
% initial conditions: 50% random initialisation
X0 = rand(1,n); X0(1,:) = X0(1,:) < 0.5;

% alternative initial conditions: all on
X0 = ones(1,n);

%================================================
state = zeros(T,n);             % pre-allocate matrix to store time-evolution of system 
X = X0(1,:);
flashes = ones(1,n);
totals = zeros(1,n);

for t = 1:T
    V = (Pow * X(K))+1;         % input value (k bits as number)
    X = diag(B(V,:));           % lookup boolean function, set X to the next state
    state(t,:) = X;
    totals = totals + X';
    if t>1
    for i = 1:n
        ischange = state(t-1, i) == state(t, i);
        flashes(1,i) = flashes(1,i) + ischange;
    end
    end    
end

%================================================
% sorting
% create row vector of sorted indices based on chosen sorting property
switch sortmode
    case 'flashes'      % how many times a node changes state over the total time
        [sorflash, sorind] = sort(flashes);
    
    case 'totals'
        [sortot, sorind] = sort(totals);
        
    otherwise
        sorind = [1:1:n];
end

sorstate = ones(size(state));   %initialize new, sorted state matrix
for i = 1:n
   sorstate(:,i) = state(:,sorind(1,i)); 
end

%================================================
% visualisation
partial = sorstate(1,:);
B = imresize(partial,5,'nearest');
imshow(B,'Border','tight');
pause(5)
for i = 2:T
    pause(0.1)
    partial = sorstate(1:i,:);
    B = imresize(partial,5,'nearest');
    imshow(B,'Border','tight');
end

end

