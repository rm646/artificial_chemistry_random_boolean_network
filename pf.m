function prob_sigmoid = pf(fraction)
% function takes argument between 0 and 1, which is used to select a value
% from a graph specified in the function (tanh here) which is compared to a
% random number, if > rand then return 1
shift = 0.5;
stretch = 6;
probability = 0.5 + 0.5*tanh(stretch*(fraction-shift));
compare_to = rand;
if probability > compare_to
    prob_sigmoid = 1;
else
    prob_sigmoid = 0;
end
end %fun