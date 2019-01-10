%A function to calculate the probability for division when the probability
%increases linearly starting from BeginProb and reches 1 at EndProb
function [Prob]= DivisionProb(BeginProb, EndProb, t, last_div)
Prob=max(0, (t-last_div)/(EndProb-BeginProb) + BeginProb/(BeginProb-EndProb));
end