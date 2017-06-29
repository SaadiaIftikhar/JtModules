%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%

[CurrentModule, CurrentModuleNum, ModuleName] = CPwhichmodule(handles);

%textVAR01 = What do you want to call the objects identified by this module?
%defaultVAR01 = Nuclei
%infotypeVAR01 = objectgroup indep
ObjectName = char(handles.Settings.VariableValues{CurrentModuleNum,1});

%textVAR02 = What did you call the intensity image that should be used for object identification?
%infotypeVAR02 = imagegroup
ImageName = char(handles.Settings.VariableValues{CurrentModuleNum,2});
%inputtypeVAR02 = popupmenu

%textVAR03 = Intensity thresholding: Threshold correction factor
%defaultVAR03 = 0.98
ThresholdCorrection = str2num(char(handles.Settings.VariableValues{CurrentModuleNum,3}));

%textVAR04 = Intensity thresholding: Lower and upper bounds on threshold, in the range [0,1]
%defaultVAR04 = 0.0019,0.02
ThresholdRange = char(handles.Settings.VariableValues{CurrentModuleNum,4});

%textVAR05 = Cutting passes (0 = no cutting)
%defaultVAR05 = 2
CuttingPasses = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,5}));

%textVAR06 = Debug mode: Show individual steps while iterating over clumped objects
%choiceVAR06 = Off
%choiceVAR06 = On
DebugMode = char(handles.Settings.VariableValues{CurrentModuleNum,6});
%inputtypeVAR06 = popupmenu

%textVAR07 = Object selection: Maximal SOLIDITY of objects, which should be cut (1 = solidity independent)
%defaultVAR07 = 0.92
SolidityThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,7}));

%textVAR08 = Object selection: Minimal FORM FACTOR of objects, which should be cut (0 = form factor independent)
%defaultVAR08 = 0.40
FormFactorThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,8}));

%textVAR09 = Object selection: Maximal AREA of objects, which should be cut (0 = area independent)
%defaultVAR09 = 500000
UpperSizeThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,9}));

%textVAR10 = Object selection: Minimal AREA of objects, which should be cut (0 = area independent)
%defaultVAR10 = 5000
LowerSizeThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,10}));

%textVAR11 = Minimal area that cut objects should have.
%defaultVAR11 = 2000
LowerSizeCutThres = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,11}));

%textVAR12 = Test mode for object selection: solidity, area, form factor
%choiceVAR12 = No
%choiceVAR12 = Yes
TestMode2 = char(handles.Settings.VariableValues{CurrentModuleNum,12});
%inputtypeVAR12 = popupmenu

%textVAR13 = Perimeter analysis: SLIDING WINDOW size for curvature calculation
%defaultVAR13 = 9
WindowSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,13}));

%textVAR14 = Perimeter analysis: FILTER SIZE for smoothing objects
%defaultVAR14 = 15
smoothingDiskSize = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,14}));

%textVAR15 = Perimeter analysis: Maximum concave region equivalent RADIUS
%defaultVAR15 = 20
PerimSegEqRadius = str2double(char(handles.Settings.VariableValues{CurrentModuleNum,15}));

%textVAR16 = Perimeter analysis: Minimum concave region equivalent CIRCULAR SEGMENT (degree)
%defaultVAR16 = 6
PerimSegEqSegment = degtorad(str2double(char(handles.Settings.VariableValues{CurrentModuleNum,16})));

%textVAR17 = Test mode for perimeter analysis: overlay curvature etc. on objects
%choiceVAR17 = No
%choiceVAR17 = Yes
TestMode = char(handles.Settings.VariableValues{CurrentModuleNum,17});
%inputtypeVAR17 = popupmenu

%%%VariableRevisionNumber = 15