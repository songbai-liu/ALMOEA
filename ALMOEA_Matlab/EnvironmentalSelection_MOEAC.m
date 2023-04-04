function [Population,z,znad] = EnvironmentalSelection_MOEAC(Population,N,z,znad,cSize)
%The environmental selection proposed in MOEAC

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    St = find(FrontNo<=MaxFNo);
    
    %% Normalization
    [PopObj,z,znad] = Normalization(Population(St).objs,z,znad);

    %% Clustering-based sorting
    tFrontNo = cNDSort(PopObj,cSize);
    
    %% Selection
    MaxFNo    = find(cumsum(hist(tFrontNo,1:max(tFrontNo)))>=N,1);
    LastFront = find(tFrontNo==MaxFNo);
    LastFront = LastFront(randperm(length(LastFront)));
    tFrontNo(LastFront(1:sum(tFrontNo<=MaxFNo)-N)) = inf;
    Next      = St(tFrontNo<=MaxFNo);
    %% Population for next generation
    Population = Population(Next);
end

function tFrontNo = cNDSort(PopObj,cSize)
% Do clustering-based sorting
    N  = size(PopObj,1);  
    normP = sum(PopObj,2); 
    fit = repmat(normP,1,cSize);
    %% Clustering
    [~,class] = clustering(PopObj,cSize); 
    %% Sort
    tFrontNo = zeros(1,N);
    for i = 1 : cSize
        C = find(class==i);
        [~,rank] = sort(fit(C,i));
        tFrontNo(C(rank)) = 1 : length(C);
    end
end

function [clusters, class] = clustering(PopObj,cSize)
    N  = size(PopObj,1);
    PopObj = PopObj./sum(PopObj,2);
    clusters = cell(1,N);
    for i = 1:N
        clusters{i} = PopObj(i);
        centroids(i) = PopObj(i);
        flag(i) = 0;
        class(i) = i;
    end
    mminD = pdist2(centroids(1), centroids(2));
    index2 = 1;
    index3 = 2;
    for i = 1:N
        rd = randi(N);
        while rd == i
            rd = randi(N);
        end
        minD = pdist2(centroids(i), centroids(rd));
        index0 = i;
        index1 = rd;
        for j = 1:N
            if i ~= j
                dist = pdist2(centroids(i), centroids(j));
                if minD > dist
                    minD = dist;
                    index0 = i;
                    index1 = j;
                end
            end
        end
        minDist(i) = minD;
        minIndex(i) = index1;
        if mminD > minD
            mminD = minD;
            index2 = index0;
            index3 = index1;
        end
    end
    nsize = N;
    while nsize>cSize
        flag(index2) = 1;
        clusters{index3} = [clusters{index3};clusters{index2}];
        centroids(index3) = mean(clusters{index3},1);
        for i = 1:N
            if(class(i) == index2)
                class(i) = index3;
            end
            if  flag(i) == 0 && (minIndex(i) == index3 || minIndex(i) == index2)
                rd = randi(N);
                while rd==i || flag(rd)==1
                    rd = randi(N);
                end
                minD = pdist2(centroids(i), centroids(rd));
                for j = 1:N
                    if i~=j && flag(j) == 0
                        ss = pdist2(centroids(i), centroids(j));
                        if minD > ss
                            minD = ss;
                            rd = j;
                        end
                    end
                end
                minIndex(i) = rd;
                minDist(i) = minD;
            end
        end
        rd = randi(N);
        while flag(rd)==1
            rd = randi(N);
        end
        fDist = minDist(rd);
        index2 = rd;
        index3 = minIndex(rd);
        for i = 1:N
            if flag(i) == 0
                if fDist > minDist(i)
                    fDist = minDist(i);
                    index2 = i;
                    index3 = minIndex(i);
                end
            end
        end
        nsize = nsize-1;
    end

    for i = N:1
        if flag(i) == 1
            clusters(i) = []; 
            for j = N:1
                if class(j) > i
                    class(j) = class(j)-1;
                end
            end
        end 
    end
    class = class';
end