clear vars
close all
clc

dim = 3;

seq = fastaread('seq_1.fasta');
seq1 = fastaread('seq_2.fasta');

a = cgrDft(seq.Sequence);
b = cgrDft(seq1.Sequence);


function [ PSZ ] = cgrDft( seq )
len = length(seq);
x(1) = 0.0; y(1) = 0.0;
z(1) = x(1) + 1i * y(1);
for j = 1:len
    
    a = seq(j);
    
    switch a
    
    case 'A'
            v = [0,0];
    case 'C'
            v = [0,1];
    case 'G'
            v = [1,0];
    case 'T'
            v = [1,1];
    end
    
    x(j+1) = 0.5 * (x(j) + v(1));
    y(j+1) = 0.5 * (y(j) + v(2));
    z(j+1) = x(j+1) + 1i * y(j+1);  
end
%[x;y]'
 figure;
 %plot(x,y);
 scatter(x,y,'.');
 xlim([0,1]);
 ylim([0,1]);
%Discrete Fourier Transforms
Z = fft(z);
%Power spectrums
PSZ = abs(Z).^2;
PSZ(1)=[]; %Exclude the first term
%figure;
%plot(PSZ);
end

function [Tm] = evenScaleVector(Tn,M)
% Function to even scale number serie Tn (length N) to Tm (length M, N<M)
% Inputs: numberical serie: Tn, length to be scaled: M
% Output: even scaled numberical series of length M
%
% Changchuan Yin
% University of Illinois at Chicago
% Email: cyinbox@gmail.com
% Last update 06/18/2015
%
% Citation:
% Hoang,T., Yin, C., & Yau, S. S. T. (2016). Numerical encoding of DNA sequences by Chaos Game Representation
% with application in similarity comparison. Genomics, Vol 107, 2016, Elsevier Inc. Journal of Theoretical Biology.
 N = length(Tn);
 Tm = zeros(1,M);
 Tm(1) = Tn(1);
 
 for k = 2:M
   Q = k*N/M;
   R = floor(k*N/M);
  
   if R == 0
      R = 1;
   end
    
   if Q == R  %Q is an integer  
      Tm(k) = Tn(Q);
   else       %Q is not an integer  
      Tm(k) = Tn(R) + (Q-R)*(Tn(R+1)-Tn(R)); 
   end
end
end 
