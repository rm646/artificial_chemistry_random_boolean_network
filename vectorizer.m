function v = vectorizer(object)
% takes an object, uses its cycle state matrix, and produces 
% a number of two dimensional column vectors from it. The number of vectors
% produced depends on the numorbs property of the rbatom. The vectors
% produced are stored in the rbatom's properties.

% set-up variables
state = object.cyclemat;
[T,n] = size(state);

numorbs = object.numorbs;
if object.isatm == 1
    dg(object);
end
D = object.D;
if isempty(D)
    disp('Empty D')
end

% crush the matrix in the vertical direction to become 2xn.
% convert first half of each column (which is equivalent to a bitstring of
% length floor(T/2)) into a number, then normalize by 2^(floor(T/2)).
% if not using base 2, but base a say, normalize with a^(...)
% Update: use last bit (if >2 bits) for sign


firstfew = floor(T/2);
lastfew = T-floor(T/2);
crush1 = zeros(2,n);            % initialize first compressed matrix
storage1 = zeros(1,firstfew);   % initialize vector to hold things before conversion to decimal
storage2 = zeros(1,lastfew);
signs = -ones(2,n);

if T<2
    disp('Cyclelength is <2, no vectors produced')
    v = 0;
    return
end

if T>3
for i = 1:n
    for t = 1:firstfew
        storage1(1,t) = state(t,i);
    end
        signs(1,i) = storage1(end);
        crush1(1,i)= (2^(1-firstfew))*bi2de(storage1(1:end-1))*((-1)^signs(1,i));        % factor of 2^ to normalise, keep last bit for sign
        
    for t = (firstfew+1):T
        storage2(1,t-firstfew) = state(t,i);
    end
        signs(2,i) = storage2(end);
        crush1(2,i) = (2^(1-lastfew))*bi2de(storage2(1:end-1))*((-1)^signs(2,i));
        
end %i

else
for i = 1:n
    for t = 1:firstfew
        storage1(1,t) = state(t,i);
    end
        crush1(1,i)= (2^(-firstfew))*bi2de(storage1);        % factor of 2^ to normalise, keep last bit for sign
    for t = (firstfew+1):T
        storage2(1,t-firstfew) = state(t,i);
    end
        crush1(2,i) = (2^(-lastfew))*bi2de(storage2);
end %i    
end %if T>3

% find average flashes property for each of the numorbs chunks
avflash = zeros(1,numorbs);
for k=1:numorbs
    flash = object.flashes;
    noderange = D(k,1):D(k,2);
    chunkofflash = flash(1,noderange);
    avflash(1,k) = sum(chunkofflash)/length(chunkofflash);
end %k

% use this to set orbtol
object.orbtol = avflash./T;


%now that we have a chunkified matrix crushed in the T direction,
%crush each chunk in the n direction too, by averaging.

orbvec = zeros(2,numorbs);      % initialize matrix to hold the final crush, columns are the vectors

for k =1:numorbs
    noderange = D(k,1):D(k,2);
    orbvec(1,k) = sum(crush1(1,noderange))/length(noderange);
    orbvec(2,k) = sum(crush1(2,noderange))/length(noderange);
end %k

% set orbvec
object.orbvec = orbvec;

% rotate according to what's been done already (i.e. history)
ur = updaterotate(object);

v = 1;
end %function