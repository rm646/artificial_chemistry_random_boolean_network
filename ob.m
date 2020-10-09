function orb_bond_test = ob(orb1,orb2,t1,t2)
% function takes orb1 and orb2 as column vectors and tests their alignment
% using tolerances t1 and t2 (rads) to draw a probability from the distribution
% pf (which takes a fraction as an argument). That probability is then used to decide whether to return 1 (bond) or
% 0 (not bond).
x = dot(orb1,orb2);
m1 = norm(orb1);
m2 = norm(orb2);
% calculate fraction of alignment if they fall within tolerance
if x >= -m1*m2 && x <= -m1*m2*cos(t1+t2)
    fa = (x + m1*m2*cos(t1+t2))/(-m1*m2);   % fraction of alignment, 100% is antiparallel, 0% is on the edge of tolerance
    if pf(fa) == 1
        orb_bond_test = 1;
    else
        orb_bond_test = 0;
    end

else
    orb_bond_test = 0;
end  %if in tolerance

end %fun
