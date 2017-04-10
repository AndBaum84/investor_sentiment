clear all
close all


%% Number of Simulation Periods
Periods = 30;

%% Parameter Variation Loop

VaryParams = 1;                             % Set 1 if parameters shall vary. Set 0 otherwise (if only one set of params to be tested.)
VariationNumber = 20;
VariationBoundOverreaction = 0.8;           % Dominance of Overreaction Pattern from -1 (dominated) to 1 (dominating)
VariationBoundUnderreaction = 0.8;           % Dominance of Underreaction Pattern from -1 (dominated) to 1 (dominating)


if VaryParams == 0 VariationNumber = 1; elseif VaryParams == 1 Periods = 300; end

for l = 1:4
    for m = 1:4
    
    VaryMain1 = [0.85 0.75 0.65 0.6];
    VaryMain2 = [0.05 0.15 0.2 0.3];
for k = 1:VariationNumber          % Variation of first Param
for j = 1:VariationNumber          % Variation of second Param
    

%% Sentiment Investor Parameters

% Probabilities for Model 1 (Mean-Reversion) --> piL1
% Probabilities for Model 2 (Momentum) --> piH1
if VaryParams == 1
    
    
    piL1 = VaryMain2(m);
    piH1 = VaryMain1(l);
    
    VariationValues1 = linspace(piL1,0.5,VariationNumber);
    VariationValues2 = linspace(piH1,0.5,VariationNumber);
    
    piL2 = VariationValues1(j);
    piL3 = 1 - piL2;
    piL4 = 1 - piL1;
    
    piH2 = VariationValues2(k);     % P(y(i) == y(i-1)) < 0.5
    piH3 = 1 - piH2;                % P(y(i) == y(i-1)) > 0.5
    piH4 = 1 - piH1;
else
    
    piL1 = 0.15;
    piL2 = 0.25;
    piL3 = 1 - piL2;
    piL4 = 1 - piL1;
    
    piH1 = 0.85;
    piH2 = 0.75;
    piH3 = 1 - piH2;
    piH4 = 1 - piH1;
end

% Transition Probabilities
% Higher probability to switch from Model 2 to Model 1 than vice versa (lambda2_1 > lambda1_1)

lambda1_1 = 0.1; % P(Switch from 1 to 2)     
lambda1_2 = 1 - lambda1_1;

lambda2_1 = 0.3; % P(Switch from 2 to 1)     
lambda2_2 = 1 - lambda2_1;

% Initial Probability to be in Model 1
clear q
q(1:2) = 0.5;

%% Sentiment Pricing Formula params p1 and p2
%PriceBehav = PriceFund + Return.*(p1 - p2 * q)

% Interest Rate (constant)
r = 0.1;

Q = [lambda1_2*piL1 lambda1_2*piL2 lambda1_2*piL3 lambda1_2*piL4 lambda2_1*piL1 lambda2_1*piL2 lambda2_1*piL3 lambda2_1*piL4;...
     lambda1_2*piL2 lambda1_2*piL3 lambda1_2*piL4 lambda1_2*piL1 lambda2_1*piL2 lambda2_1*piL3 lambda2_1*piL4 lambda2_1*piL1;...
     lambda1_2*piL3 lambda1_2*piL4 lambda1_2*piL1 lambda1_2*piL2 lambda2_1*piL3 lambda2_1*piL4 lambda2_1*piL1 lambda2_1*piL2;...
     lambda1_2*piL4 lambda1_2*piL1 lambda1_2*piL2 lambda1_2*piL3 lambda2_1*piL4 lambda2_1*piL1 lambda2_1*piL2 lambda2_1*piL3;...     
     lambda1_1*piH1 lambda1_1*piH2 lambda1_1*piH3 lambda1_1*piH4 lambda2_2*piH1 lambda2_2*piH2 lambda2_2*piH3 lambda2_2*piH4;...
     lambda1_1*piH2 lambda1_1*piH3 lambda1_1*piH4 lambda1_1*piH1 lambda2_2*piH2 lambda2_2*piH3 lambda2_2*piH4 lambda2_2*piH1;...
     lambda1_1*piH3 lambda1_1*piH4 lambda1_1*piH1 lambda1_1*piH2 lambda2_2*piH3 lambda2_2*piH4 lambda2_2*piH1 lambda2_2*piH2;...
     lambda1_1*piH4 lambda1_1*piH1 lambda1_1*piH2 lambda1_1*piH3 lambda2_2*piH4 lambda2_2*piH1 lambda2_2*piH2 lambda2_2*piH3];     

I = eye(8,8)*(1+r);

gamma0 = [1 1 -1 -1 1 1 -1 -1];
gamma1 = [0 0 0 0 1 1 0 0];
gamma2 = [1 1 0 0 -1 -1 0 0];

p1 = (1/r) * (gamma0 * (1+r) * ((I - Q) \ Q) * gamma1');
p2 = -(1/r) * (gamma0 * (1+r) * ((I - Q) \ Q) * gamma2');

%% Definition of Fundamentals

% Size of Dividend Changes
ChangeSize = 3;

% Starting Dividend and rational Price
Dividend(1:2) = 30;
PriceFund(1:2) = Dividend(1:2) / r;
PriceBehav(1:2) = PriceFund(1:2);

% Fake initial change
DivChange(1:2) = ChangeSize;

%% Simulation Loop
for i = 3:Periods

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
    if sign(DivChange(i-1)) == sign(DivChange(i-2))
    
        q(i) = ((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL1) / (((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL1) + ((lambda1_1 * q(i-1) + lambda2_2 * (1 - q(i-1))) * piH1));
    
    elseif sign(DivChange(i-1)) ~= sign(DivChange(i-2))
    
        q(i) = ((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL2) / (((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL2) + ((lambda1_1 * q(i-1) + lambda2_2 * (1 - q(i-1))) * piH2));
    
    end
    
elseif sign(DivChange(i-1)) ~= sign(DivChange(i))
    if sign(DivChange(i-1)) ~= sign(DivChange(i-2))
    
        q(i) = ((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL3) / (((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL3) + ((lambda1_1 * q(i-1) + lambda2_2 * (1 - q(i-1))) * piH3));
    
    elseif sign(DivChange(i-1)) == sign(DivChange(i-2))
    
        q(i) = ((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL4) / (((lambda1_2 * q(i-1) + lambda2_1 * (1 - q(i-1))) * piL4) + ((lambda1_1 * q(i-1) + lambda2_2 * (1 - q(i-1))) * piH4));
    
    end

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

%{
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
%}
%% Plot for Impact of Variation

VariationName2 = '\pi_L_2';        % piL2 --> Index k
VariationName1 = '\pi_H_2';        % piH2 --> Index j

figure(l)
subplot(2,2,m)
spy(VariationUnderOverExistence')
%title(['Under/Overreaction for various ',VariationName1,' and ',VariationName2])
xlabel(VariationName2,'FontSize',16)
ylabel(VariationName1,'FontSize',16)
legend([num2str(round(mean(mean(VariationUnderOverExistence)),2)*100),'%'],'location','north')
ax = gca;
ax.XLim = [-(piL1/(0.5-piL1))*VariationNumber-0.5 VariationNumber+0.5];
ax.XTick = [-(piL1/(0.5-piL1))*VariationNumber 0 VariationNumber];
ax.XTickLabel = {'0',['\pi_L_1 = ',num2str(piL1)],'0.5'};
ax.YLim = [-(1-piH1)/(piH1-0.5)*VariationNumber-0.5 VariationNumber+0.5];
ax.YTick = [-(1-piH1)/(piH1-0.5)*VariationNumber 0 VariationNumber];
ax.YTickLabel = {'1',['\pi_H_1 = ',num2str(piH1)],'0.5'};
rectangle('Position',[0 0 21 21])

end
end 