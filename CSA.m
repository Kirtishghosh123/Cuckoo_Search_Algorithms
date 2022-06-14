format short
clear all
clc

%% Pharse 1 : Input Parameter

D = 2; %Dimension fo the problem
lb = [-5, -5]; %Lower Bound of the variable
ub = [5 5]; %Upper bound of the variable
N = 20; %Population Size
n = N; 
pa = 0.25; %Discovery Rate of the alien eggs and solutions
max_iter = 100; %Maximum number of iteration



%% Pharse 3 : Generate Initial Population

for i = 1 : N
    for j = 1 : D
        nest(i, j) = lb(:,j) + rand.*(ub(:,j) - lb(:,j));
    end
end

fx =  fns(nest);
beta = 3/2;
sigma = (gamma(1 + beta)*sin(pi*beta/2)/(gamma((1 + beta)/2)*beta*2^((beta - 1)/2)))^(1/beta);



%% Pharse 4 : Cuckoo Search Main Loop Search

for iter = 3 : max_iter
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


    [fmin,K1] = min(fx);
    best =  nest(K1, :);

    %% Replace some nest by constructing new solution

    K = rand(size(nest)) < 0.25;
    stepsizeK = rand*(nest(randperm(n), :) - nest(randperm(n), :));
    new_nest = nest + stepsizeK.*K;
    


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
        %% Perform Greedy Selection
        fnew = fns(s);
        if fnew < fx(ii, :)
            nest(ii, :) = s;
            fx(ii, :) = fnew;
        end
    end

    

    
end

