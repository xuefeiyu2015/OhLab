## Project goal

Test whether object value and location value were presented by the same population or different among frontal eye field(FEF) and superior colliculus(SC)

## Task

- Object Value Task/ Stable Value Task

The task starts with an fixation point, after around 800ms, one fractal will appear either on the left or the right visual field. After a 400-500ms delay, the fixation point will disappear, indicating the goal signal. Then the monkey need make an saccadic eyemovement to the target with 400ms, after gazing onto the fractal for 400ms, the monkey will be given reward, the target will continue to appear for ~200ms before disappear. 

There are two types of objects: high value objects--good objects, which will give monkey high amount of reward; low value objects--bad objects, which will give monkey low amount of reward. Each session typically has 4 high value objects and 4 low value objects

Behavior quantification for object value task: measure the reaction time, peak velocity for saccades to fractals, and gaze duration after saccade. Those saccade parameters should be different between high and low value objects.

The batch file for all object value task is @/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/FEFSC_BothRA_StableValue.txt

- Location Value Task/one direction reward task

The task procedural is exactly the same as object value task, but the value association is with the spatial location of the object rather than the object feature. For one block, when a target appeared on the left, the target was associated with high reward (good location), whereas when a target appeared on the right, the target was associated with low reward (bad location). After 15-20 trials, the value-location association will be reversed.

Behavior quantification for location value task is the same as object value task.

The batch file for all location value is @/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/FEFSC_BothRA_OneDR.txt

## Neural analysis

All the analysis are focus on the responses inside the receptive field(RF). Only include neurons with has both object value task and location value task.

The batch file for object value for FEF neurons is @/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/FEF_BothRA_StableValue.txt

The batch file for object value for SC neurons is @/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/SC_BothRA_StableValue.txt

-- single unit level:

Use ROC analysis to quantify value responses difference inside RF. Compare the ROC between high value and low value responses inside RF for both object value and location value.

And compare whether they are the consistent or different, whether the same group of neurons respond to both or different value type was represented by different popuation of neurons.

The batch file for location value for FEF neurons is @/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/FEF_BothRA_StableValue.txt

The batch file for location value for SC neurons is @/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/SC_BothRA_StableValue.txt

-- population level:

Construct decoders separately for object value and location value

1. Construct within value type value axis within object value condition and location value condition, for example, randomly chose particular number of trials with high value condition and trials with low value condition and construct value coding axis for each value type (object value v.s location value).

2. Project the population responses for object value onto object value axis and location value axis separately. And check the value reponses separation on each value axis for two value types.

3. Calculate the within value type axis angle and between value type axis angle, within value axis anlge should be small, close to zero, and between value type axis angle should be large, approaching 90 deg. Which indicate the two value axis are almost orthogonal

-- Circuit level optogenetics

We have used retrograde optogenetics to identify sc projecting FEF neurons and FEF receipent SC neurons. We may repeat the above analysis for those projection neurons.



