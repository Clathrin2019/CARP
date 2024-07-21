Installation:
1. Add the CMEAnalysis package to the Matlab path.
2. Set this folder as the working directory.

Usage:
>> FusionFrequency('usefulChNum', 3)
% Suitable for three channels: EGFP-Sensor, AP2-TagRFP, and CLTA-nano670.
>> FusionFrequency('usefulChNum', 2, 'AP2ch', 2)
% Suitable for two channels. If AP2 is in the 561 channel, set the parameter to 2; if in the 647 channel, set it to 3.

Examples:
>> FusionFrequency('usefulChNum', 2, 'AP2ch', 2) % AP2-mScarlet-I
>> FusionFrequency('usefulChNum', 2, 'AP2ch', 3) % AP2-nano670

