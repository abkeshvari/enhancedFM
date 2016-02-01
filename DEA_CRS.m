%%% Please cite the paper properly if you use the code. 
%%% "Keshvari, Abolfazl. 2016. An Enhanced Fourier-Motzkin Method for DEA."
%%% Define the probelm here
load('data.mat'); % Load electricity companies dataset
units=(1:size(data,1))';%The list of units.
m1=1;m2=3; %Specify the number of inputs (m1) and outputs (m2)

uniq=4; %Accuracy setting
algorithms=3; %1: only FM algorithm, 2: only Qhull, 3: both
clc;
if algorithms==1 || algorithms==3
    %%% Run FM algorithm
    D= (data(units,:));
    fprintf('Running enhanced FM... \n')
    tic;
    [Tm]=FM_ver1(D,m1,m2);
    fprintf('Elapsed time for enhanced FM is %f seconds.\n',toc)
    
    T=roundn(normr(Tm(:,1:m1+m2)),-uniq);T(all(T(:,1:m2)==0,2),:)=[];
   %%% Calculate efficiency scores, and the unique facets from FM algorithm
   %%% Uncomment the lines below to get the efficiency scores
%     EffN=Tm(:,1:m2)*data(units,m1+1:end)';
%     EffD=Tm(:,m2+1:m2+m1)*data(units,1:m1)';
%     Eff=EffN./EffD;
%     Eff_FM=max(Eff,[],1)';% %Efficiency scores
end;

if algorithms==2 || algorithms==3
    %%% Initialize Qhull. Make sure the folder "qhull-2012.1" is located next to the script file
    bigM=10^10;
    %Normalize the data: use the commane normr to normalize data, or remove
    %it to use the original data
    D= normr(data(units,:)); 
    Mat=[D;zeros(1,m1+m2);bigM*eye(m1) zeros(m1,m2); zeros(m2,m1) -bigM*eye(m2)];
    dlmwrite('.\qhull-2012.1\data', [m1+m2;length(units)+m1+m2+1] );
    save '.\qhull-2012.1\data' Mat -ASCII  '-append'
    
    %%% Run Qhull
    fprintf('Running Qhull... \n')
    tic;
    [status,results] =dos('.\qhull-2012.1\bin\qconvex.exe s n <".\qhull-2012.1\data" >".\qhull-2012.1\res"');
    %     disp(results)%uncomment this line to see the command line results of Qhull
    fprintf('Elapsed time for Qhull method is %f seconds.\n',toc)
    
    %%% Read the result of Qhull and calculate efficicency scores
    A=dlmread('.\qhull-2012.1\res' );A([1,2],:)=[];
    A(abs(A(:,m1+m2+1))>0.000000000001,:)=[];
    Tq=[A(:,m1+1:m1+m2) -A(:,1:m1);];  Tq=Tq(m1+1:end,:);
    TQ=roundn(normr(Tq(:,1:m1+m2)),-uniq);TQ(all(TQ(:,1:m2)==0,2),:)=[];
    
    %%% Calculate efficiency scores, and the unique facets from Qhull
    %%% Uncomment the lines below to get the efficiency scores
%     EffN=Tq(:,1:m2)*data(units,m1+1:end)';
%     EffD=Tq(:,m2+1:m2+m1)*data(units,1:m1)';
%     Eff=EffN./EffD;
%     Eff_Q=max(Eff,[],1)'; %Efficiency scores
end;
fprintf('\n ****\nNormal vectors: matrix T for enhanced FM, matrix TQ for Qhull.\n')
fprintf('The first m2 columns are coefficients of outputs, and the next m1 columns are coefficients of inputs.\n')
% clearvars A D Eff EffD EffN Mat Tm Tq results status ;
