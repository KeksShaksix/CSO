function Hd = Filter
%FILTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.14 and Signal Processing Toolbox 9.2.
% Generated on: 23-Oct-2024 17:40:06

Numerator   = [0.028 0.053 0.071 0.053 0.028];  % Numerator coefficient
                                                % vector
Denominator = [1 -2.026 2.148 -1.159 0.279];    % Denominator coefficient
                                                % vector

Hd = dfilt.df2t(Numerator, Denominator);


% [EOF]
