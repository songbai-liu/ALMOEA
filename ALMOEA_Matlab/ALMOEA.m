classdef ALMOEA < ALGORITHM
% <multi/many> <real> <large>
% tec --- 3 --- type of environmental selection. 1 = NSGA-II (Default), 2 = NSGA-III, 3 = TDEA, 4 = MOEAC

%------------------------------- Reference --------------------------------
% S. Liu, J. Li, Q. Lin, Y. Tian and K. C. Tan, Learning to accelerate evolutionary search for large-scale multiobjective optimization, IEEE Transactions
% on Evolutionary Computation, 2023, 27(1): 67-81.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            tec = Algorithm.ParameterSet(3);

            %% Generate the reference points and random population
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population     = Problem.Initialization();
            [z,znad]      = deal(min(Population.objs),max(Population.objs));

            %% Optimization
            while Algorithm.NotTerminated(Population)
                mlp = ModelLearning(Population);
                Offspring = DEgenerator(Problem,Population,mlp);
                if tec==1
                    Population = EnvironmentalSelection_NSGAII([Population,Offspring],Problem.N);
                end

                if tec==2
                    z       = min([z;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                    Population = EnvironmentalSelection_NSGAIII([Population,Offspring],Problem.N,W,z);
                end

                if tec==3
                    [Population,z,znad] = EnvironmentalSelection_TDEA([Population,Offspring],W,Problem.N,z,znad);
                end

                if tec==4
                    [Population,z,znad] = EnvironmentalSelection_MOEAC([Population,Offspring],Problem.N,z,znad,Problem.N);
                end
            end
        end
    end
end