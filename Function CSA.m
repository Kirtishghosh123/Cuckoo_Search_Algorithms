format short
clear all;
clc;
---------------------------------------------

%% Pharse  : Input Parameter
D = 3;
lb = [-1 -1 -1];
ub = [1 1 1];
N = 10;
n = N;
pa = 0.25
max_iter =100;
----------------------------------------------

%% Pharse 2 : Defining Objective Function
  
function out = fns(X)
    x = X(:,1);
    y = X(:,2);
    z = X(:,3);
    out = (sin(x)./x).*(sin(y)./y).*(sin(z)./z); 
end
=---------------------------------------------------

%% Pharse 3 : Generate Initial Population Range
for i = 1 : N
    for j = 1 : D
        nest(i, j) = lb(:,j) + rand.*(ub(:,j) - lb(:,j))
    end
end
fx = fns(nest);

beta = 3/2;
sigma = (gamma(1 + beta)*sin(pi*beta/2)/(gamma((1 + beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
-------------------------------------------------------------------------------------------------------

for iter = 1 : max_iter
    [fnv, indf] = min(fx);
    best = nest(indf, :);
    
    for j =  1 : N
        s = nest(j,:);
        X = s;
        %%%Levy Flights by Mantegna's algorithms
        u =  randn(size(s))*sigma;
        v = randn(size(s));
        step = u./abs(v).^(1/beta);
        Xnew = X + randn(size(s)).*0.01.*step.*(X - best);
        
        %%Check Bounds
        for kk = 1: size(Xnew, 2)
            if Xnew(kk) > ub(kk)
                Xnew(kk) = ub(kk);
            elseif Xnew(kk) < lb(kk)s
                Xnew(kk) = lb(kk);
            end
        end
        
        %%Perform Greedy Selection
        fnew = fns(Xnew);
        if(fnew<fx(j,:))
            nest(j,:) = Xnew;
            fx(j,:) = fnew;
        end
    end

end
------------------------------------------------------------------------

%% Find the current best
[fmin, K1] = min(fx);
best = nest(K1, :);

-------------------------------------------------------------------------

%% Replace Some nest by constructing new solution
K = rand(size(nest)) < pa;

stepsizeK = rand*(nest(randperm(n), :) - nest(randperm(n), :));
new_nest = nest + stepsizeK.*K;
-----------------------------------------------------------------------

%% Check Bounds
for ii = 1 : size(nest, 1)
    s =new_nest(ii, :);
    for kk = 1 : size(s, 2)
        if s(kk) > ub(kk)
            s(kk) = ub(kk);
        elseif s(kk) < lb(kk)
            s(kk) = lb(kk);
        end
    end
    new_nest(ii,:) = s;
-------------------------------------------------------------------------

%% Perform Greedy Selection
    fnew = fns(s);
    if fnew < fx(ii, :)
        nest(ii, :) = s;
        fx(ii, :) = fnew;
    end
end


%% Finding the optimal Value
[optval, optind] = min(fx);

BestFx(iter) = optval;
BestX(iter, :) = nest(optind, :);


%% Plotting the result
plot(BestFx, 'LineWidth', 2);
xlabel('Iteration Number')
ylabel('Fitness Value')
title('Convergence Vs Iteration')
grid on


