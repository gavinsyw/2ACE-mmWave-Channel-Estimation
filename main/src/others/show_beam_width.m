%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2019 Yi Zhang and The University of Texas at Austin 
%  
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code or any (modified) part of it in any publication,
% please cite:
%
% Yi Zhang, Kartik Patel, Sanjay Shakkottai, and Robert W. Heath Jr.. 2019. 
% Side-information-aided Non-coherent Beam Alignment Design for Millimeter 
% Wave Systems. In MobiHoc '19: The Twentieth ACM International Symposium 
% on Mobile Ad Hoc Networking and Computing, July 02-05, 2019, Catania, 
% Italy. ACM, New York, NY, USA, 10 pages.
%
% Author: Yi Zhang
% Contact email: yi.zhang.cn@utexas.edu 
% Last modified: Apr. 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function description:
% This function plots the directional beam patterns generated by a ULA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% ULA: a structure array that groups the related parameters on ULA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function show_beam_width(ULA)
    lambda = ULA.lambda;
    d = ULA.d;
    N = ULA.Nt;
    G = ULA.NQt;
    Phase_Bit = 2;
    Grid = [0:180/G:180*(G-1)/G];
    figure
    for i = round(80/(180/G))
        Beampattern = transpose(exp(1i*2*pi*d/lambda*cosd(Grid(i)).*(0:1:N-1)));
        BeamformerQuant = Quantize_PS(Beampattern,Phase_Bit);
        rad_range = pi/3:0.0005:pi*2/3;
        gain = zeros(length(rad_range),1);
        ind = 1;
        for rad = rad_range
            Channel = transpose(exp(1i*2*pi*d/lambda*cos(rad).*(0:1:N-1)));
            gain(ind) = abs(BeamformerQuant'*Channel);
            gain(find(gain<=0.001)) = 0.001;
            ind = ind + 1;
        end
        Gain_dB = 10*log10(gain);
        M_Gain = max(Gain_dB);
        M_Gain_minus3dB = M_Gain-3;
        plot(rad2deg(rad_range),M_Gain*ones(1,length(rad_range)));hold on;
        plot(rad2deg(rad_range),M_Gain_minus3dB*ones(1,length(rad_range)));hold on;
        plot(rad2deg(rad_range),Gain_dB,'r'); 
        hold on;
    end
end