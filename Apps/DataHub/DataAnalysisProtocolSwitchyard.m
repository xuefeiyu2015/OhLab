function DataAnalysisProtocolSwitchyard(Data, StartTrial, EndTrial,Protocol,EditIndex,ShowFigureFlag,OutputFlag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Basic analysis
if strcmp(Protocol,'Eye movement analysis')
    if EditIndex
        edit EyeMovementAnalysis.m
    else
        EyeMovementAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
    
elseif strcmp(Protocol,'Performance analysis')
    if EditIndex
        edit PerformanceAnalysis.m
    else
        PerformanceAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Performance analysis(Foraging)')
    if EditIndex
        edit PerformanceAnalysis_Foraging.m
    else
        PerformanceAnalysis_Foraging(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Neuron Basic Analysis')
    if EditIndex
        edit Neuron_Basic_Analysis.m
    else
       Neuron_Basic_Analysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis for optostim related tasks%%%%
elseif strcmp(Protocol,'FreeViewEyeMovement')
    if EditIndex
        edit FreeViewEyeMovement.m
    else
       FreeViewEyeMovement(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'FreeViewEyeMovement_txtexport')
    if EditIndex
        edit FreeViewEyeMovement_txtexport.m
    else
       FreeViewEyeMovement_txtexport(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'EyeTrace_rawIllustration')
    if EditIndex
        edit FreeViewEyeSaccadeSelectionIllustration.m
    else
       FreeViewEyeSaccadeSelectionIllustration(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end


elseif strcmp(Protocol,'FreeViewEyeMovement_Vari')
    if EditIndex
        edit FreeViewEyeMovement_Vari.m
    else
       FreeViewEyeMovement_Vari(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
elseif strcmp(Protocol,'FreeViewNeuroResponse')
    if EditIndex
        edit FreeViewNeuroResponse.m
    else
       FreeViewNeuroResponse(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'FreeViewNeuroResponse_ExampleExport')
    if EditIndex
        edit FreeViewNeuroResponse_ExampleExport.m
    else
       FreeViewNeuroResponse_ExampleExport(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'FreeViewEyeMovement_ExampleExport')
    if EditIndex
        edit FreeViewEyeMovement_ExampleExport.m
    else
       FreeViewEyeMovement_ExampleExport(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'OptoWaveformPlot')
    if EditIndex
        edit OptoWaveformPlot.m
    else
       OptoWaveformPlot(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

    

elseif strcmp(Protocol,'FreeViewNeuroResponse_Pop')
    if EditIndex
        edit FreeViewNeuroResponse_Pop.m
    else
       FreeViewNeuroResponse_Pop(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'FreeViewNeuroResponse_Pop_ExampleExport')
    if EditIndex
        edit FreeViewNeuroResponse_Pop_ExampleExport.m
    else
       FreeViewNeuroResponse_Pop_ExampleExport(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'FreeView_Animation')
    if EditIndex
        edit FreeView_Animation.m
    else
       FreeView_Animation(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'VaryAngleStimulation')
    if EditIndex
        edit VaryAngleStimulation.m
    else
       VaryAngleStimulation(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end



    
elseif strcmp(Protocol,'FreeViewNeuroResponse_Vari')
    if EditIndex
        edit FreeViewNeuroResponse_Vari.m
    else
       FreeViewNeuroResponse_Vari(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'FreeViewNeuroResponse_VariPop')
    if EditIndex
        edit FreeViewNeuroResponse_VariPop.m
    else
       FreeViewNeuroResponse_VariPop(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'FreeViewNeuroResponse_RedLight')
    if EditIndex
        edit FreeViewNeuroResponse_RedLight.m
    else
       FreeViewNeuroResponse_RedLight(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

    
elseif strcmp(Protocol,'LaserMeasurement')
    if EditIndex
        edit LaserMeasurement.m
    else
       LaserMeasurement(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'DifferentPeriodStimEffect')
    if EditIndex
        edit DifferentPeriodStimEffect.m
    else
       DifferentPeriodStimEffect(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'SaccadeStimEffect')
    if EditIndex
        edit SaccadeStimEffect.m
    else
       SaccadeStimEffect(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'SaccadeStimEffectExample')
    if EditIndex
        edit SaccadeStimEffectExample.m
    else
       SaccadeStimEffectExample(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'FractalStimEffect')
    if EditIndex
        edit FractalStimEffect.m
    else
       FractalStimEffect(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'SaccadeTuningCompare')
    if EditIndex
        edit SaccadeTuningCompare.m
    else
       SaccadeTuningCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ITI_FirstSaccades')
    if EditIndex
        edit ITI_FirstSaccades.m
    else
       ITI_FirstSaccades(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'VisMotorResponses_Population')
    if EditIndex
        edit VisMotorResponses_Population.m
    else
       VisMotorResponses_Population(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
%%%%%%%%%%%%%%%%%%%% Analysis for Luminance Tuning Task %%%%%%%%%%%%%%%
elseif strcmp(Protocol,'LuminanceTuning_Stim')
    if EditIndex
        edit LuminanceTuning_Stim.m
    else
       LuminanceTuning_Stim(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

%%%%%%%%%%%%%%%%% Analysis for Delay Saccade Stimulation Task
%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%
elseif strcmp(Protocol,'DelaySaccadeVaryPeriodStim')
    if EditIndex
        edit DelaySaccadeVaryPeriodStim.m
    else
       DelaySaccadeVaryPeriodStim(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end


    
%%%%%%%%%%%%%%%%% Analysis for memory saccade %%%%%%%%%%
elseif strcmp(Protocol,'MemorySaccadeAnalysis')
    if EditIndex
        edit MemorySaccadeAnalysis.m
    else
       MemorySaccadeAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'MemorySaccadePopulation')
    if EditIndex
        edit MemorySaccadePopulation.m
    else
       MemorySaccadePopulation(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

%%%%%%%%%%%% Analysis for one DR %%%%
elseif strcmp(Protocol,'OneDRAnalysis')
    if EditIndex
        edit OneDRAnalysis.m
    else
       OneDRAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'OneDR_RT')
    if EditIndex
        edit OneDRAnalysis.m
    else
       OneDRRT(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end    
elseif strcmp(Protocol,'VisualLatencyMeasure')
    if EditIndex
        edit VisualLatencyMeasure.m
    else
       VisualLatencyMeasure(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis for path tasks(squence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tasks)%%%%%
elseif strcmp(Protocol,'Sequence NeuralAnalysis Object')
    if EditIndex
        edit Sequence_NeuralAnalysis_Object.m
    else
         Sequence_NeuralAnalysis_Object(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    

elseif strcmp(Protocol,'LearningOfSequence NeuralResponse')
    if EditIndex
        edit Learning_NeuralResponse.m
    else
         Learning_NeuralResponse(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'LearningOfSequence')
    if EditIndex
        edit LearningOfSequence.m
    else
         LearningOfSequencee(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
elseif strcmp(Protocol,'Sequence Behavior Analysis')
    if EditIndex
        edit Sequence_BehaviorAnalysis.m
    else
         Sequence_BehaviorAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis for object sequence tasks%%%%%
elseif strcmp(Protocol,'ObjectSeq NeuralAnalysis')
    if EditIndex
        edit ObjectSeq_NeuralAnalysis.m
    else
         ObjectSeq_NeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ObjSeq BehaviorAnalysis')
    if EditIndex
        edit ObjSeq_BehaviorAnalysis.m
    else
         ObjSeq_BehaviorAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'VideoReconstructionForPath')
    if EditIndex
        edit VideoReconstructionForPath.m
    else
         VideoReconstructionForPath(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis for the passive tasks    
elseif strcmp(Protocol,'FractalTuning_Passive')
    if EditIndex
        edit FractalTuning_Passive.m
    else
         FractalTuning_Passive(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
elseif strcmp(Protocol,'EachLocFractalTuning_Passive')
    if EditIndex
        edit EachLocFractalTuning_Passive.m
    else
         EachLocFractalTuning_Passive(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
    
elseif strcmp(Protocol,'GoodBadInRF')
    if EditIndex
        edit GoodBadInRF.m
    else
         GoodBadInRF(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
    

    
elseif strcmp(Protocol,'RF_Mapping_Task PSTH')
    if EditIndex
        edit RF_Mapping_Task_PSTH.m
    else
        RF_Mapping_Task_PSTH(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'SceneTuningPassiveView')
    if EditIndex
        edit  SceneTuningPassiveView.m
    else
      SceneTuningPassiveView(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'CatogrySceneResponsePassiveView')
    if EditIndex
        edit  CatogrySceneResponsePassiveView.m
    else
      CatogrySceneResponsePassiveView(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'SceneFractalTuning_PassiveView')
    if EditIndex
        edit  SceneFractalTuning_PassiveView.m
    else
      SceneFractalTuning_PassiveView(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
  
elseif strcmp(Protocol,'ITI_ScenePassiveView')
    if EditIndex
        edit  ITI_ScenePassiveView.m
    else
      ITI_ScenePassiveView(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'SceneTuningSequenceFreeView')
    if EditIndex
        edit  SceneTuningSequenceFreeView.m
    else
      SceneTuningSequenceFreeView(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis for the object directed saccade tasks
elseif strcmp(Protocol,'RTChange_Inactivation')
    if EditIndex
        edit RTChange_Inactivation.m
    else
        RTChange_Inactivation(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'DifferentPeriodCompare')
    if EditIndex
        edit DifferentPeriodCompare.m
    else
        DifferentPeriodCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'ObjSaccadeAnalysis')
    if EditIndex
        edit ObjSaccadeAnalysis.m
    else
        ObjSaccadeAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ObjectValueLocationRegression')
    if EditIndex
        edit ObjectValueLocationRegression.m
    else
        ObjectValueLocationRegression(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'ObjSaccadeAnalysis_TrialHistory')    
    if EditIndex
        edit ObjSaccadeAnalysis_TrialHistory.m
    else
        ObjSaccadeAnalysis_TrialHistory(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ObjectResponsesDifferentLocation')
    if EditIndex
        edit ObjectResponsesDifferentLocation.m
    else
        ObjectResponsesDifferentLocation(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'RTGoodBadAnalysis')
    if EditIndex
        edit RTGoodBadAnalysis.m
    else
        RTGoodBadAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
 elseif strcmp(Protocol,'ObjectOnsetLatency')
    if EditIndex
        edit ObjectOnsetLatency.m
    else
        ObjectOnsetLatency(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis for the flexible value task%%%%%
elseif strcmp(Protocol,'FlexibleValueResponseAnalysis')
    if EditIndex
        edit FlexibleValueResponseAnalysis.m
    else
        FlexibleValueResponseAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis for the object memory saccade tasks   
 elseif strcmp(Protocol,'ObjectMemorySaccadeBehavior')
    if EditIndex
        edit ObjectMemorySaccadeBehavior.m
    else
        ObjectMemorySaccadeBehavior(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis for the sequence memory task%%%%%%%%%%%%%%%%%%%%
 elseif strcmp(Protocol,'SequenceMemoryBehavior')
    if EditIndex
        edit SequenceMemoryBehavior.m
    else
        SequenceMemoryBehavior(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end   
 
 
 
    
    
%%%%%%%%%%%%%%%%%%%%%%Analysis for the delay match to sample object tasks
elseif strcmp(Protocol,'DelayObjectResponsesAnalysis')
    if EditIndex
        edit DelayObjectResponsesAnalysis.m
    else
        DelayObjectResponsesAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'DelayObjectTestsResponsesAnalysis')
    if EditIndex
        edit DelayObjectTestResponsesAnalysis.m
    else
        DelayObjectTestResponsesAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
  
elseif strcmp(Protocol,'MemoryBehaviorAnalysis')
    if EditIndex
        edit MemoryBehaviorAnalysis.m
    else
        MemoryBehaviorAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'ShortLong_MemoryBehaviorAnalysis')
    if EditIndex
        edit ShortLong_MemoryBehaviorAnalysis.m
    else
        ShortLong_MemoryBehaviorAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis for the free viewing task
elseif strcmp(Protocol,'SceneTuningFreeView')
    if EditIndex
        edit SceneTuningFreeView.m
    else
       SceneTuningFreeView(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'CatogrySceneResponseFreeView')
     if EditIndex
        edit CatogrySceneResponseFreeView.m
    else
       CatogrySceneResponseFreeView(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis for the object memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% foraging task

elseif strcmp(Protocol,'SceneObjecrtMemoryAnalysis')
    if EditIndex
        edit SceneObjecrtMemoryAnalysis.m
    else
        SceneObjecrtMemoryAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis for the foraging tasks    
elseif strcmp(Protocol,'Foraging Task PSTH')
    if EditIndex
        edit ForagingTask_PSTH.m
    else
        ForagingTask_PSTH(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'Foraging Task Location Selectivity')
    if EditIndex
        edit ForagingTask_LocationSelectivity.m
    else
       ForagingTask_LocationSelectivity(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'SceneTuning')
    if EditIndex
        edit SceneTuning.m
    else
       SceneTuning(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end   

elseif strcmp(Protocol,'FractalTuning')
    if EditIndex
        edit FractalTuning.m
    else
       FractalTuning(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol, 'ForagingRasterCheck')
    if EditIndex
        edit ForagingRasterCheck.m
    else
        ForagingRasterCheck(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);     
    end
elseif strcmp(Protocol,'BonusTrialAnalysis')
    if EditIndex
        edit BonusTrialAnalysis.m
    else
       BonusTrialAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'SceneFractalAssociationCheck')
    if EditIndex
        edit SceneFractalAssociationCheck.m
    else
       SceneFractalAssociationCheck(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
    
    
elseif strcmp(Protocol,'MotivationAndReactionTime_Analysis')
    if EditIndex
        edit MotivationAndReactionTime_Analysis.m
    else
      MotivationAndReactionTime_Analysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'RandomSeqFixSeqCompare')
    if EditIndex
        edit RandomSeqFixSeqCompare.m
    else
      RandomSeqFixSeqCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ITI_RandomSeqFixSeqCompare')
    if EditIndex
        edit ITI_RandomSeqFixSeqCompare.m
    else
        ITI_RandomSeqFixSeqCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'EventSegregationAnalysis')
    if EditIndex
        edit EventSegregationAnalysis.m
    else
        EventSegregationAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'RandomFixOrderSceneSep')
    if EditIndex
        edit RandomFixOrderSceneSep.m
    else
        RandomFixOrderSceneSep(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'TargetSaccadeAnalysis')
    if EditIndex
        edit TargetSaccadeAnalysis.m
    else
        TargetSaccadeAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'SaccadeAnalysis')
    if EditIndex
        edit SaccadeAnalysis.m
    else
       SaccadeAnalysis(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'RewardResponseAnalysis')
    if EditIndex
       edit RewardResponseAnalysis.m
    else       
       RewardResponseAnalysis(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'FractalOrderNeuralAnalysis')
    if EditIndex
       edit FractalOrderNeuralAnalysis.m
    else       
       FractalOrderNeuralAnalysis(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
 %%%%%%%%%%%%%%%%%%%%%%%%RichPoorSceneSeq
elseif strcmp(Protocol,'SequenceCompareEachBlock')
    if EditIndex
       edit SequenceCompareEachBlock.m
    else       
       SequenceCompareEachBlock(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Foraging Centeral Periperal Task
elseif strcmp(Protocol,'CenterPeriperalTargetResponseAnalysis')
    if EditIndex
        edit CenterPeriperalTargetResponseAnalysis.m
    else
        CenterPeriperalTargetResponseAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Below is for foraging explore task
elseif strcmp(Protocol,'ExploreBehaviorAnalysis')
    if EditIndex
        edit ExploreBehaviorAnalysis.m
    else
        ExploreBehaviorAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ExploreLearningCurve')
    if EditIndex
        edit ExploreLearningCurve.m
    else
        ExploreLearningCurve(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ExploreBehaviorPuipleAnalysis')
    if EditIndex
        edit ExploreBehaviorPuipleAnalysis.m
    else
        ExploreBehaviorPuipleAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
    
elseif strcmp(Protocol,'GazeLocationSelectivityAnalysis')
    if EditIndex
        edit GazeLocationSelectivityAnalysis.m
    else
       GazeLocationSelectivityAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
elseif strcmp(Protocol,'CrossModeNeuralAnalysis')
    if EditIndex
        edit CrossModeNeuralAnalysis.m
    else
        CrossModeNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'FractalEffortNeuralAnalysis')
    if EditIndex
        edit FractalEffortNeuralAnalysis.m
    else
        FractalEffortNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'EventSegAcrossModeWithRT')
    if EditIndex
        edit EventSegAcrossModeWithRT.m
    else
        EventSegAcrossModeWithRT(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
elseif strcmp(Protocol,'SceneOnResponseAcrossModes')
    if EditIndex
        edit SceneOnResponseAcrossModes.m
    else
        SceneOnResponseAcrossModes(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'FixOffPredictabilityNeuralAnalysis')
    if EditIndex
        edit FixOffPredictabilityNeuralAnalysis.m
    else
        FixOffPredictabilityNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
    
elseif strcmp(Protocol,'ExploreNeuralAnalysis')
    if EditIndex
        edit ExploreNeuralAnalysis.m
    else
        ExploreNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ExploreRTNeuralCorrelatesAnalysis')
    if EditIndex
        edit ExploreRTNeuralCorrelatesAnalysis.m
    else
        ExploreRTNeuralCorrelatesAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ExploreSacGazNeuralAnalysis')
    if EditIndex
        edit ExploreSacGazNeuralAnalysis.m
    else
        ExploreSacGazNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ExploreFilterSacNeuralAnalysis')
    if EditIndex
        edit ExploreFilterSacNeuralAnalysis.m
    else
        ExploreFilterSacNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
elseif strcmp(Protocol,'SceneOnExplore_NeuralAnalysis')
    if EditIndex
        edit SceneOnExploreNeuralAnalysis.m
    else
       SceneOnExploreNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Explore_NeuralRaster_Check')
    if EditIndex
        edit ExploreRasterCheck.m
    else
        ExploreRasterCheck(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    
elseif strcmp(Protocol,'VideoReConstructionForExplore')
    if EditIndex
        edit VideoReConstructionForExplore.m
    else
        VideoReConstructionForExplore(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
   
elseif strcmp(Protocol,'ExploreSaccadeAnalysis')
    if EditIndex
        edit ExploreSaccadeAnalysis.m
    else
        ExploreSaccadeAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'RejectReactionAnalysis')
    if EditIndex
        edit RejectReactionAnalysis.m
    else
        RejectReactionAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'BeforeAfterSceneResponsesWithSeq')
    if EditIndex
        edit BeforeAfterSceneResponsesWithSeq.m
    else
        BeforeAfterSceneResponsesWithSeq(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ExploreBehaviorAnalysis_VarEye')
    if EditIndex
        edit ExploreBehaviorAnalysis_VarEye.m
    else
        ExploreBehaviorAnalysis_VarEye(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ExploreBehaviorScrambleCompare')
    if EditIndex
        edit ExploreBehaviorScrambleCompare.m
    else
        ExploreBehaviorScrambleCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'FourBlockTypeSceneResponseCompare')
    if EditIndex
        edit FourBlockTypeSceneResponseCompare.m
    else
        FourBlockTypeSceneResponseCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'ITI_CorrAnalysis')
    if EditIndex
        edit ITI_CorrAnalysis.m
    else
        ITI_CorrAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'RepetitionChangeAnalysis')
    if EditIndex
        edit RepetitionChangeAnalysis.m
    else
        RepetitionChangeAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'SceneObjResponseSummary')
    if EditIndex
        edit SceneObjResponseSummary.m
    else
        SceneObjResponseSummary(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

    
    
%New learning analysis
elseif strcmp(Protocol,'NewLearningAnalysis')
    if EditIndex
        edit NewLearningAnalysis.m
    else
        NewLearningAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
%%%%%%%%%%%%Analysis for the Free Foraging Task %%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Protocol,'FreeForagingBehaviorAnalysis')
    if EditIndex
        edit FreeForagingBehaviorAnalysis.m
    else
        FreeForagingBehaviorAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'GazeFieldMeasurementAnalysis')
    if EditIndex
        edit GazeFieldMeasurementAnalysis.m
    else
        GazeFieldMeasurementAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
%%%%%%%%%%%%%%%%%%%%%Analysis for the reward order task %%%%%%%%%%%%%%%%%
elseif strcmp(Protocol,'RewardOrderChoiceAnalysis')
    if EditIndex
        edit RewardOrderChoiceAnalysis.m
    else
        RewardOrderChoiceAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'RewardOrderNeuralAnalysis')
    if EditIndex
        edit RewardOrderNeuralAnalysis.m
    else
        RewardOrderNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end



    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(Protocol,'Batch_SaccadeAnalysis') 
    if EditIndex
       edit Batch_SaccadeAnalysis.m
    else       
       Batch_SaccadeAnalysis(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'Batch_SceneTuning_Foraging') 
    if EditIndex
       edit SceneTuning_Batch_Foraging.m
    else       
       SceneTuning_Batch_Foraging(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
    

elseif strcmp(Protocol,'NeuronClassification_Foraging') 
    if EditIndex
        edit NeuronClassification_Foraging.m
    else
        NeuronClassification_Foraging(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Batch_FractalTuning_Foraging')
  if EditIndex
        edit Batch_FractalTuning_Foraging.m
  else
        Batch_FractalTuning_Foraging(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
  end
  
elseif strcmp(Protocol,'Batch_SceneFractalCorrelation_Foraging')
    if EditIndex
        edit Batch_SceneFractalCorrelation_Foraging.m
    else
        Batch_SceneFractalCorrelation_Foraging(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'Batch_Saccade_Foraging')
    
    if EditIndex
        edit Batch_Saccade_Foraging.m
    else
        Batch_Saccade_Foraging(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'Batch_EventRT_Exploration')
    if EditIndex
        edit Batch_EventRT_Exploration.m
    else
        Batch_EventRT_Exploration(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'Batch_FourBlockCmp_Exploration')
    if EditIndex
        edit Batch_FourBlockCmp_Exploration.m
    else
        Batch_FourBlockCmp_Exploration(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'Batch_VideoFreeView_Neural')
    if EditIndex
        edit Batch_VideoFreeView_Neural.m
    else
        Batch_VideoFreeView_Neural(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'Batch_MemorySaccade')
    if EditIndex
        edit Batch_MemorySaccade.m
    else
        Batch_MemorySaccade(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'Batch_MemorySaccadeLatency')
    if EditIndex
        edit Batch_MemorySaccadeLatency.m
    else
        Batch_MemorySaccadeLatency(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'Batch_ObjectLatency')
    if EditIndex
        edit Batch_ObjectLatency.m
    else
        Batch_ObjectLatency(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end

elseif strcmp(Protocol,'Batch_SequenceLearning')
    if EditIndex
        edit Batch_SequenceLearning.m
    else
        Batch_SequenceLearning(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end


elseif strcmp(Protocol,'Batch_FEF_SC_VideoFreeView')
    if EditIndex
        edit Batch_FEF_SC_VideoFreeView.m
    else
        Batch_FEF_SC_VideoFreeView(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'Batch_FEF_SC_VideoFreeView_ExampleExport')
    if EditIndex
        edit Batch_FEF_SC_VideoFreeView_ExampleExport.m
    else
        Batch_FEF_SC_VideoFreeView_ExampleExport(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
   
  elseif strcmp(Protocol,'Batch_OptogeneticsFreeView')
    if EditIndex
        edit Batch_OptogeneticsFreeView.m
    else
        Batch_OptogeneticsFreeView(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag); 
    end

  elseif strcmp(Protocol,'Batch_OptogeneticsFreeView_SC')
    if EditIndex
        edit Batch_OptogeneticsFreeView_SC.m
    else
        Batch_OptogeneticsFreeView_SC(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag); 
    end


elseif strcmp(Protocol,'Batch_Optogenetics')
    if EditIndex
        edit Batch_Optogenetics.m
    else
        Batch_Optogenetics(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end

elseif strcmp(Protocol,'Batch_Optogenetics_FEF')
    if EditIndex
        edit Batch_Optogenetics_FEF.m
    else
        Batch_Optogenetics_FEF(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'Batch_Optogenetics_FEF_Depth')
    if EditIndex
        edit Batch_Optogenetics_FEF_Depth.m
    else
        Batch_Optogenetics_FEF_Depth(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end

elseif strcmp(Protocol,'Batch_Optogenetics_SNr')
    if EditIndex
        edit Batch_Optogenetics_SNr.m
    else
        Batch_Optogenetics_SNr(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end


elseif strcmp(Protocol,'Batch_Opto_ValueTasks')
    if EditIndex
        edit Batch_Opto_ValueTasks.m
    else
        Batch_Opto_ValueTasks(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end

elseif strcmp(Protocol,'Batch_FreeViewBehavior')
    if EditIndex
        edit Batch_FreeViewBehavior.m
    else
        Batch_FreeViewBehavior(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end
elseif strcmp(Protocol,'Batch_FreeViewBehavior_Grid')
    if EditIndex
        edit Batch_FreeViewBehavior_Grid.m
    else
        Batch_FreeViewBehavior_Grid(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end


elseif strcmp(Protocol,'Batch_FreeViewBehavior_VaryDurAmp')
    if EditIndex
        edit Batch_FreeViewBehavior_VaryDurAmp.m
    else
        Batch_FreeViewBehavior_VaryDurAmp(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
        
    end

elseif strcmp(Protocol,'Batch_OneDR')
    if EditIndex
        edit Batch_OneDR.m
    else
        Batch_OneDR(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Batch_OptoOneDR')
    if EditIndex
        edit Batch_OptoOneDR.m
    else
        Batch_OptoOneDR(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Batch_StableValue_Opto')
    if EditIndex
        edit Batch_StableValue.m
    else
        Batch_StableValue(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Batch_FlexibleValueOpto')
    if EditIndex
        edit Batch_FlexibleValueOpto.m
    else
        Batch_FlexibleValueOpto(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Batch_StableValue_Behavior')
    if EditIndex
        edit Batch_StableValue_Behavior.m
    else
        Batch_StableValue_Behavior(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Batch_StableOneDR_Opto')
    if EditIndex
        edit Batch_StableOneDR_Opto.m
    else
        Batch_StableOneDR_Opto(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'Batch_FixStimGapStimCompare')
    if EditIndex
        edit Batch_FixStimGapStimCompare.m
    else
        Batch_FixStimGapStimCompare(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Batch_DelaySaccadeStim')
    if EditIndex
        edit Batch_DelaySaccadeStim.m
    else
        Batch_DelaySaccadeStim(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Batch_Opto_WaveFormClustering')
    if EditIndex
        edit Batch_Opto_WaveFormClustering.m
    else
        Batch_Opto_WaveFormClustering(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Batch_OneDR_Behavior')
    if EditIndex
        edit Batch_OneDR_Behavior.m
    else
        Batch_OneDR_Behavior(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'Batch_Optogenetics_Control')
    if EditIndex
        edit Batch_Optogenetics_Control.m
    else
        Batch_Optogenetics_Control(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end
elseif strcmp(Protocol,'Batch_StableValue_All')
    if EditIndex
        edit Batch_StableValue_All.m
    else
        Batch_StableValue_All(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end

elseif strcmp(Protocol,'Batch_FixStimFixation')
    if EditIndex
        edit Batch_FixStimGapStimCompare.m
    else
        Batch_FixStimFixation(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
    end




   
end





end