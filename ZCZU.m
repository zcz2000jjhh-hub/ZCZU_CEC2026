classdef ZCZU < ALGORITHM
% <2024> <multi> <real/integer/label/binary/permutation> <constrained>
% Dual-population with dynamic constraint processing and resource allocating

%------------------------------- Reference --------------------------------
% K. Qiao, Z. Chen, B. Qu, K. Yu, C. Yue, K. Chen, and J. Liang. A dual-
% population evolutionary algorithm based on dynamic constraint processing
% and resources allocation for constrained multi-objective optimization
% problems. Expert Systems With Applications, 2024, 238: 121707.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    methods
        function main(Algorithm,Problem)

            %% CEC 2026 CMOP recording variables
            global CEC26_MinIGD CEC26_MCV CEC26_run CEC26_point CEC26_BestIGD CEC26_PF

            %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();

            %% Ideal points for IMTCMO-style offspring generation
            Zmin1 = min(Population1.objs,[],1);
            Zmin2 = min(Population2.objs,[],1);

            %% Initialize reference PF and record initial population
            if ~isempty(CEC26_MinIGD)

                if isempty(CEC26_PF)
                    try
                        CEC26_PF = Problem.GetOptimum(100);
                    catch
                        try
                            CEC26_PF = Problem.optimum;
                        catch
                            error('Cannot obtain reference PF. Please manually set global CEC26_PF before running DPCPRA3.');
                        end
                    end
                end

                CEC26_BestIGD = inf;

                % Point 1: After Initialization
                RecordCEC26CMOP(Population1);

                % If the initialization has already used enough FEs, fill the corresponding sampling point.
                while CEC26_point <= 1001 && Problem.FE >= (CEC26_point-1)*200
                    RecordCEC26CMOP(Population1);
                end
            end

            Fitness1           = CalFitness_pop1(Population1.objs,Population1.cons);
            current_cons       = 0;
            gen                = 0;
            last_gen           = 100;
            change_threshold   = 1e-2;
            change_rate        = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            priority           = [];
            flag               = 0;
            constraint_handing = 0;
            archive            = Population2;
            Fitness2           = CalFitness_pop2(Population2.objs,Population2.cons,priority,current_cons,constraint_handing);
            success_rate1      = 0.5;

            %% Optimization
            while Algorithm.NotTerminated(Population1)

                if flag == 0
                    change_rate = Normalization(Population2,change_rate,ceil(Problem.FE/Problem.N));
                    if Convertion(change_rate,ceil(Problem.FE/Problem.N),gen,last_gen,change_threshold)
                        flag = 1;
                        [priority,evaluatedasible_rate] = Constraint_priority(Population2);
                        Population2 = Problem.Initialization();

                        % Reset auxiliary ideal point after reinitialization
                        Zmin2 = min(Population2.objs,[],1);

                        Fitness2 = CalFitness_pop2(Population2.objs,Population2.cons,priority,current_cons,constraint_handing);
                    end
                else
                    % Judge whether to enter next stage
                    if current_cons == 0
                        CV = Population2.cons;
                        CV = CV(:,priority(1));
                        if length(find(CV>0))/Problem.N > 0
                            current_cons = current_cons + 1;
                            gen = ceil(Problem.FE/Problem.N) + 1;
                        end
                    elseif current_cons <= size(Population2.cons,2)
                        if constraint_handing == 0
                            change_rate = Normalization(Population2,change_rate,ceil(Problem.FE/Problem.N));
                            if Convertion(change_rate,ceil(Problem.FE/Problem.N),gen,last_gen,change_threshold)
                                if current_cons<size(Population2.cons,2) && evaluatedasible_rate(priority(current_cons+1))~=1
                                    current_cons = current_cons+1;
                                elseif current_cons<size(Population2.cons,2) && evaluatedasible_rate(priority(current_cons+1))==1
                                    current_cons = size(Population2.cons,2);
                                elseif current_cons == size(Population2.cons,2)
                                    constraint_handing = 1;
                                end
                                if size(archive,2) == Problem.N
                                    Population2 = archive;
                                    Fitness2    = CalFitness_pop2(Population2.objs,Population2.cons,priority,current_cons,constraint_handing);
                                else
                                    archive     = Archive([archive,Population2],Problem.N,priority,current_cons,size(archive,2));
                                    Population2 = archive;
                                    Fitness2    = CalFitness_pop2(Population2.objs,Population2.cons,priority,current_cons,constraint_handing);
                                end

                                % Update auxiliary ideal point after archive replacement
                                Zmin2 = min(Population2.objs,[],1);

                                gen = ceil(Problem.FE/Problem.N) + 1;
                            end
                        end
                    end
                end

                %% Offspring generation with IMTCMO-style DE operators
                if flag == 0

                    % Keep the original DPCPRA offspring amount:
                    % each population generates Problem.N/2 offspring.
                    OffNum1 = Problem.N/2;
                    OffNum2 = Problem.N/2;

                    Offspring1 = GenerateIMTCMOOffspring(Problem,Population1,OffNum1,Fitness1,Zmin1);
                    Offspring2 = GenerateIMTCMOOffspring(Problem,Population2,OffNum2,Fitness2,Zmin2);

                else

                    % Keep the original resource allocation rule of DPCPRA.
                    if mod(ceil(Problem.N*success_rate1),2)~=0
                        ParentNum1 = ceil(Problem.N*success_rate1) + 1;
                    else
                        ParentNum1 = ceil(Problem.N*success_rate1);
                    end

                    % The original OperatorGAhalf produces half as many offspring as parents.
                    OffNum1 = ParentNum1/2;
                    Offspring1 = GenerateIMTCMOOffspring(Problem,Population1,OffNum1,Fitness1,Zmin1);

                    ParentNum2 = Problem.N - 2*length(Offspring1);
                    OffNum2    = ParentNum2/2;

                    if OffNum2 > 0
                        Offspring2 = GenerateIMTCMOOffspring(Problem,Population2,OffNum2,Fitness2,Zmin2);
                    else
                        Offspring2 = [];
                    end

                end

                %% Update ideal points
                if ~isempty(Offspring1)
                    Zmin1 = min([Zmin1;Offspring1.objs],[],1);
                end
                if ~isempty(Offspring2)
                    Zmin2 = min([Zmin2;Offspring2.objs],[],1);
                end

                %% Update external archive
                if flag==1 && constraint_handing~=1
                    archive = Archive([Offspring2,archive],Problem.N,priority,current_cons);
                end

                %% Environmental selection
                [Population1,Fitness1,success_rate1] = EnvironmentalSelection_pop1( ...
                    [Population1,Offspring1,Offspring2],Problem.N,true,length(Offspring1));

                [Population2,Fitness2,success_rate2] = EnvironmentalSelection_pop2( ...
                    [Population2,Offspring2,Offspring1],Problem.N,priority,current_cons,constraint_handing,length(Offspring2));

                success_rate1 = success_rate1/(success_rate1+success_rate2);

                %% Record CEC 2026 CMOP results every 200 FEs
                if ~isempty(CEC26_MinIGD)
                    while CEC26_point <= 1001 && Problem.FE >= (CEC26_point-1)*200
                        RecordCEC26CMOP(Population1);
                    end
                end

            end
        end
    end
end


function Offspring = GenerateIMTCMOOffspring(Problem,Population,OffNum,Fitness,Zmin)
% Generate offspring using the IMTCMO-style offspring generation mechanism:
%   1) OperatorDE_rand_1
%   2) OperatorDE_pbest_1_main
%
% The total number of generated offspring is OffNum.
% This function preserves DPCPRA's original offspring amount and resource allocation rule.

    OffNum = floor(OffNum);
    PopNum = length(Population);

    if OffNum <= 0 || PopNum <= 3
        Offspring = [];
        return;
    end

    % Split the offspring number into DE/rand/1 and DE/pbest/1 parts
    RandNum  = floor(OffNum/2);
    PbestNum = OffNum - RandNum;

    Offspring = [];

    %% Part 1: DE/rand/1 with neighbor pairing strategy
    if RandNum > 0
        MatingPool = [Population(randsample(PopNum,PopNum))];
        [Mate1,Mate2,Mate3] = Neighbor_Pairing_Strategy(MatingPool,Zmin);

        Offspring_rand = OperatorDE_rand_1(Problem, ...
            Mate1(1:RandNum), Mate2(1:RandNum), Mate3(1:RandNum));

        Offspring = [Offspring,Offspring_rand];
    end

    %% Part 2: DE/pbest/1
    if PbestNum > 0
        Offspring_pbest = OperatorDE_pbest_1_main( ...
            Population, PbestNum, Problem, Fitness, 0.1);

        Offspring = [Offspring,Offspring_pbest];
    end

end


function RecordCEC26CMOP(Population)
% Record Min_IGD and MCV for CEC 2026 CMOP.
%
% point = 1    : after initialization
% point = 2    : 200 FEs
% point = 3    : 400 FEs
% ...
% point = 1001 : 200000 FEs

    global CEC26_MinIGD CEC26_MCV CEC26_run CEC26_point CEC26_BestIGD

    if isempty(CEC26_MinIGD)
        return;
    end

    if CEC26_point > 1001
        return;
    end

    %% MCV
    CEC26_MCV(CEC26_point,CEC26_run) = CalPopMCV(Population);

    %% Current IGD
    currentIGD = CalPopIGD(Population);

    %% Min_IGD
    if isnan(currentIGD)
        CEC26_MinIGD(CEC26_point,CEC26_run) = NaN;
    else
        CEC26_BestIGD = min(CEC26_BestIGD,currentIGD);

        if CEC26_BestIGD < 1.0e-8
            CEC26_MinIGD(CEC26_point,CEC26_run) = 0;
        else
            CEC26_MinIGD(CEC26_point,CEC26_run) = CEC26_BestIGD;
        end
    end

    CEC26_point = CEC26_point + 1;

end


function MCV = CalPopMCV(Population)
% Calculate mean constraint violation.

    cons = Population.cons;

    if isempty(cons)
        MCV = 0;
        return;
    end

    cons(cons < 0) = 0;
    CV  = sum(cons,2);
    MCV = mean(CV);

end


function value = CalPopIGD(Population)
% Calculate IGD using feasible individuals.
% If no feasible solution exists, return NaN.

    global CEC26_PF

    cons = Population.cons;

    if isempty(cons)
        CV = zeros(length(Population),1);
    else
        cons(cons < 0) = 0;
        CV = sum(cons,2);
    end

    feasible = CV <= 0;

    if sum(feasible) == 0
        value = NaN;
        return;
    end

    FeasiblePop = Population(feasible);
    FeasibleObj = FeasiblePop.objs;

    %% Use the first non-dominated front among feasible solutions
    try
        FrontNo = NDSort(FeasibleObj,1);
        FeasibleObj = FeasibleObj(FrontNo == 1,:);
    catch
        % If NDSort is unavailable, use all feasible solutions
    end

    %% Use at most 100 feasible individuals
    if size(FeasibleObj,1) > 100
        try
            CrowdDis = CrowdingDistance(FeasibleObj);
            [~,rank] = sort(CrowdDis,'descend');
            FeasibleObj = FeasibleObj(rank(1:100),:);
        catch
            FeasibleObj = FeasibleObj(1:100,:);
        end
    end

    value = CalculateIGD(FeasibleObj,CEC26_PF);

end


function score = CalculateIGD(PopObj,PF)
% IGD = average distance from each point in PF to the nearest point in PopObj.

    if isempty(PopObj) || isempty(PF)
        score = NaN;
        return;
    end

    Distance = pdist2(PF,PopObj);
    score    = mean(min(Distance,[],2));

end