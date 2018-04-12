%% evaluate a generalized plant by plugging in theta value 
% by Pan Zhao, March 31, 2014
% control input: 2;
% output: 4, 2 for error signal, 2 for control input;
% measured output: 2;
function Gaev = AugPltEval(Ga,Par)
%% Ga is augmented plant with symbol variable;
%% Ganv is augmented plant evaluated
if size(Par,1) == 1
    Gaev.A = double(subs(Ga.A,Par));   
    Gaev.B1 = double(subs(Ga.B1,Par));   
    Gaev.B2 = double(subs(Ga.B2,Par));   
    Gaev.C1 = double(subs(Ga.C1,Par));   
    Gaev.D11 = double(subs(Ga.D11,Par));   
    Gaev.D12 = double(subs(Ga.D12,Par));   
    Gaev.C2 = double(subs(Ga.C2,Par));   
    Gaev.D21 = double(subs(Ga.D21,Par));   
    Gaev.D22 = double(subs(Ga.D22,Par));   
elseif size(Par,1) == 2
    theta1 = Par(1,1); theta2 = Par(2,1);
    Gaev.A = double(subs(Ga.A));   
    Gaev.B1 = double(subs(Ga.B1));   
    Gaev.B2 = double(subs(Ga.B2));   
    Gaev.C1 = double(subs(Ga.C1));   
    Gaev.D11 = double(subs(Ga.D11));   
    Gaev.D12 = double(subs(Ga.D12));   
    Gaev.C2 = double(subs(Ga.C2));   
    Gaev.D21 = double(subs(Ga.D21));   
    Gaev.D22 = double(subs(Ga.D22));   
elseif size(Par,1) == 3
    theta1 = Par(1,1); theta2 = Par(2,1); theta3 = Par(3,1);
    Gaev.A = double(subs(Ga.A));   
    Gaev.B1 = double(subs(Ga.B1));   
    Gaev.B2 = double(subs(Ga.B2));   
    Gaev.C1 = double(subs(Ga.C1));   
    Gaev.D11 = double(subs(Ga.D11));   
    Gaev.D12 = double(subs(Ga.D12));   
    Gaev.C2 = double(subs(Ga.C2));   
    Gaev.D21 = double(subs(Ga.D21));   
    Gaev.D22 = double(subs(Ga.D22));   
end
