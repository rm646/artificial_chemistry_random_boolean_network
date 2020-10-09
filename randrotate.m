function v2 = randrotate(object, temp)
%Given a matrix of 2D column vectors to rotate, and temp to affect amount of rotation
%Function rotates v1 randomly to output v2

ri = randi(2);
o1 = object.orbvec;
o2 = zeros(size(o1));
x = 2*pi*temp*rand*(-1)^ri;      %random angle, temp should scale from 0-->1, BEWARE LINEAR SCALE
object.rot = object.rot + x;
R = [cos(x),-sin(x);sin(x),cos(x)];
[~,max] = size(o1);

for j = 1:max
    o2(:,j) = R*o1(:,j);
end

object.orbvec = o2;
v2 = x;
end

