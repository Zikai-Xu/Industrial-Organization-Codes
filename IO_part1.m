%% 1) Reduced-form regressions
% a
% Read data
clc
clear
T = readtable("analysis_data.csv");
%% Run a linear regression
lm = fitlm(T,"investment~region+naics_recode+gravity+HPV_recorded+DAV")

% b
%% Collapse data
tbl = removevars(T,"gravity_str");
varfun(@mean,tbl,'GroupingVariables',{'region','naics_recode','gravity','quarter'});
%% Run a regression using collapsed data
clm = fitlm(tbl,"compliance~inspection+violation+fine")

%% 2) Evidence of dynamic enforcement
% a
ilm = fitlm(T,"inspection~lag_violator_notHPV*lag_HPV_status*DAV")

%% b
blm = fitlm(T,"inspection~lag_violator_notHPV*lag_HPV_status*DAV+region+naics_recode+gravity")

%% c
cclm = fitlm(T,"fine~lag_HPV_status*DAV+region+naics_recode+gravity+inspection+violation")

%% 3)
% a
% Create a table for plants not in compliance
tnic = T(~(T.compliance~=0),:);
% Regression
lm3a = fitlm(tnic,"fine~lag_investment*lag_HPV_status*DAV+region+naics_recode+gravity")

%% b
% Expected fine given investment and regular violator
lag_investment = 1;
naics_recode = 1;
region = 1;
gravity = mean(tnic.gravity);
DAV = mean(tnic.DAV);
lag_HPV_status = 0;
Xnew = table(lag_HPV_status,lag_investment,naics_recode,region,gravity,DAV);
fine_i_v = lm3a.predict(Xnew)

% Expected fine given investment and HPV
Xnew.lag_HPV_status = 1;
fine_i_h = lm3a.predict(Xnew)

% Expected fine given no investment and regular violator
Xnew.lag_HPV_status = 0;
Xnew.lag_investment = 0;
fine_n_h = lm3a.predict(Xnew)

% Expected fine given no investment and HPV
Xnew.lag_HPV_status = 1;
fine_n_h = lm3a.predict(Xnew)

    

%% c
X = tnic;
X.lag_investment(:) = 1;
Xinv = X;
X.lag_investment(:) = 0;
Xnoinv = X;

lag_investment = Xinv.lag_investment;
naics_recode = Xinv.naics_recode;
region = Xinv.region;
gravity = Xinv.gravity;
DAV = Xinv.DAV;
lag_HPV_status = Xinv.lag_HPV_status;
Xinv = table(lag_investment,lag_HPV_status,naics_recode,region,gravity,DAV);

lag_investment = Xnoinv.lag_investment;
naics_recode = Xnoinv.naics_recode;
region = Xnoinv.region;
gravity = Xnoinv.gravity;
DAV = Xnoinv.DAV;
lag_HPV_status = Xnoinv.lag_HPV_status;
Xnoinv = table(lag_investment,lag_HPV_status,naics_recode,region,gravity,DAV);
diffine = lm3a.predict(Xinv)-lm3a.predict(Xnoinv);

%% Probit regression
[probitCoef,dev] = glmfit(diffine,tnic.investment,'binomial',"link",'probit')
x = linspace(min(diffine),max(diffine),104446);
probitFit = glmval(probitCoef,x,'probit');
plot(diffine,tnic.investment,'bs', x,probitFit,'r-');
xlabel('Dif in Fine'); ylabel('Investment');

