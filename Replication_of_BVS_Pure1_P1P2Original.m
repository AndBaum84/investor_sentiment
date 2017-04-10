clear all
close all


%% Number of Simulation Periods
Periods = 10;

%% Parameter Variation Loop

VaryParams = 0;                             % Set 1 if parameters shall vary. Set 0 otherwise (if only one set of params to be tested.)
VariationNumber = 20;
VariationPiL1Lower = 0;
VariationPiL1Upper = 0.5;
VariationPiH1Lower = 0.5;
VariationPiH1Upper = 1;
VariationBoundOverreaction = 0.7;           % Dominance of Overreaction Pattern from -1 (dominated) to 1 (dominating)
VariationBoundUnderreaction = 0.7;           % Dominance of Underreaction Pattern from -1 (dominated) to 1 (dominating)
VariationName1 = 'piH1';
VariationName2 = 'piL1';

VariationValues1 = linspace(VariationPiL1Lower,VariationPiL1Upper,VariationNumber);
VariationValues2 = linspace(VariationPiH1Upper,VariationPiH1Lower,VariationNumber);

if VaryParams == 0 VariationNumber = 1; end

for k = 1:VariationNumber          % Variation of first Param
for j = 1:VariationNumber          % Variation of second Param
    

%% Sentiment Investor Parameters

% Probabilities for Model 1 (Mean-Reversion) --> piL1
% Probabilities for Model 2 (Momentum) --> piH1
if VaryParams == 1
    piL1 = VariationValues1(j);
    piH1 = VariationValues2(k);
else
    piL1 = 0.15;     % P(y(i) == y(i-1))
    piH1 = 0.85;     % P(y(i) == y(i-1))
end
    piL2 = 1 - piL1;
    piH2 = 1 - piH1;

% Transition Probabilities
% Higher probability to switch from Model 2 to Model 1 than vice versa (lambda2_1 > lambda1_1)

lambda1_1 = 0.1; % P(Switch from 1 to 2)     
lambda1_2 = 1 - lambda1_1;

lambda2_1 = 0.3; % P(Switch from 2 to 1)     
lambda2_2 = 1 - lambda2_1;

% Initial Probability to be in Model 1
clear q
q(1) = 0.5;

%% Sentiment Pricing Formula params p1 and p2
%PriceBehav = PriceFund + Return.*(p1 - p2 * q)

% Interest Rate (constant)
r = 0.1;

Q = [lambda1_2*piL1 lambda1_2*piL2 lambda2_1*piL1 lambda2_1*piL2;...
     lambda1_2*piL2 lambda1_2*piL1 lambda2_1*piL2 lambda2_1*piL1;...
     lambda1_1*piH1 lambda1_1*piH2 lambda2_2*piH1 lambda2_2*piH2;...
     lambda1_1*piH2 lambda1_1*piH1 lambda2_2*piH2 lambda2_2*piH1];

I = eye(4,4)*(1+r);

gamma0 = [1 -1 1 -1];
gamma1 = [0 0 1 0];
gamma2 = [1 0 -1 0];

p1 = (1/r) * (gamma0 * (1+r) * ((I - Q) \ Q) * gamma1');
p2 = -(1/r) * (gamma0 * (1+r) * ((I - Q) \ Q) * gamma2');

%% Definition of Fundamentals

% Size of Dividend Changes
ChangeSize = 3;

% Starting Dividend and rational Price
Dividend(1) = 30;
PriceFund(1) = Dividend(1) / r;
PriceBehav(1) = PriceFund(1);

% Fake initial change
DivChange(1) = ChangeSize;

%% Simulation Loop
for i = 2:Periods

if rand < 0.5
   y = 1;
else
   y = -1;
end

% Fundamentals
DivChange(i) = y * ChangeSize;
Dividend(i) = Dividend(i-1) + DivChange(i);
PriceFund(i) = Dividend(i) / r;

if sign(DivChange(i-1)) == sign(DivChange(i))
    
    q(i) = ((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL1) / (((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL1) + ((lambda1_1 * q(i-1) + lambda2_2 * (1 - q(i-1))) * piH1));
    
elseif sign(DivChange(i-1)) ~= sign(DivChange(i))
    
    q(i) = ((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL2) / (((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL2) + ((lambda1_1 * q(i-1) + lambda2_2 * (1 - q(i-1))) * piH2));
    
end

% Sentiment investor Pricing
PriceBehav(i) = PriceFund(i) + DivChange(i).*(p1 - p2 * q(i));

end


%% Measuring Pricing Differences

PriceDeviation = (PriceFund - PriceBehav) .* sign(DivChange);
UnderOverRatio(j,k) = mean(sign(PriceDeviation));


%% Variation Output

VariationOutput(j,k) = UnderOverRatio(j,k);
if VariationOutput(j,k) < VariationBoundUnderreaction && VariationOutput(j,k) > -VariationBoundOverreaction && piL1 < 0.5 && piH1 > 0.5
    VariationUnderOverExistence(j,k) = 1;
else
    VariationUnderOverExistence(j,k) = 0;
end

end
end

figure(1)
subplot(3,1,1)
plot([1:Periods],PriceFund,'-o',[1:Periods],PriceBehav,'r -o')
legend('Fundamental Price','Sentiment Price','location','southoutside','Orientation','horizontal')
%ylim([PriceFund(1) - ChangeSize/r*Periods*(1/3) PriceFund(1) + ChangeSize/r*Periods*(1/3)])
subplot(3,1,2)
plot([1:Periods],PriceDeviation,[1:Periods],zeros(1,Periods),'r')
ylabel(['Under(+)/Over(-) ',num2str(round(UnderOverRatio(j,k),2))])
subplot(3,1,3)
plot(q,'-o')

%% Plot for Impact of Variation

figure(2)
spy(VariationUnderOverExistence')
title(['Existence of Under/Over for param Variation ', VariationName1,' and ',VariationName2,' within dominance intervall ',num2str([-VariationBoundOverreaction VariationBoundUnderreaction])])
xlabel(VariationName2)
ylabel(VariationName1)