
ffasearch

Implementation of the Fast Folding Algorithm. For more information
about description and implementation, see Appendix of my paper
Kondratiev et al. 2009, ApJ, 702, 692 
(http://adsabs.harvard.edu/abs/2009ApJ...702..692K)

INSTALL: just run 'make' command

(c) Vlad Kondratiev, email: vlad dot kondratiev at gmail dot com


Fast Folding Algorithm (FFA):
------------------------------

1. Originally developped by D.H. Staelin (Proceedings of the IEEE, 
April 1969, vol. 57, issue 4, 724-725)
Number of samples N and base period P0 can be any.
The only thing required is log2(N/P0) should be integer, or N/P0 is power of 2,
but there are words that algorithm can be modified for other values of N.
He gives example of FFA, for N=12 and P0=3.
Search is doing between P0 and P0+1 (in samples)

2. Implementation of FFA by Lovelace (Lovelace et al., 1969)
N is power of 2, N=2^k
P0 is power of 2, P0=2^l
and for sure N/P0 is also power of 2
Search is doing between P and P+1 , where P is between P0/2 and P0


Ideas:
-------

In general, if we have N samples in raw-data (tim file) and we want
to run FFA for the periods between P0 (in samples) and P0+1 than
the ratio M = N/P0 should be power of 2. If it is not the case then
someone needs to get less samples N from the file to make this true.
To search for a pulsar in a broader range of periods from P1 to P2
you have to run FFA for every integer period P0 from this range.

Searching for a periods between P0 and P0+1, 
period resolution is deltaT/M, where deltaT is the time resolution of the
input time series, and M is the N/P0 and it is power of 2. For different
periods between P1 and P2, M also could be different, smaller for larger
periods.

Coding the FFA we can implement different kinds of rebinning.
1) decimate the raw-data stream (tim-file) to reduce original time resolution.
This is good when the pulsar is strong enough and pulse-width is large.
It is especially worth while to do when raw-file is large and we need to save 
some memory. As I understood correctly, this is what is implemented in 
Michael' ffa program with -n option.
This is implemented also in my ffasearch program (-n or --nbins option)

After this step we will have time resolution of _nbins_ times worse than the original
time resolution. Lets mark it DeltaT1 = DeltaT * nbins;

With this DeltaT1 we are starting FFA... and we will have period resolution in FFA of
equal to DeltaT1 / M, depeneding on M.

2) if we do not actually need the period resolution of DeltaT1 / M and if we want to
have it worse or fix it, we can act in two ways:

a) as we already talked before about this...
   we can read not all the file but just a block of data from it. Then we will have less
   number of samples N, and hence smaller M, and hence worse period resolution.
   After, we can read next block from the file, search again, calculating SNRs of the
   candidates as square root of the sum of their squares from every block.

   This is partly implemented. Now it is possible to fix the period resolution
   with --periodres or --mres options. But only first block from the file will be read.

b) alternatively, we can rebin the data in the FFA itself having hot just P0 samples but
   less. But this rebinning will be different from the decimating of the raw-data, because
   here we will rebin every period P0 keeping the phase.
   In this case the speed of the algorithm will increase, because we won't need to make
   calculations for every sample in PO. If we want to add together _mbins_ then we will
   have just P0/_mbins_ samples in the period.
   Now, we will have worse time resolution DeltaT2 = DeltaT1 * _mbins_
   and hence the period resolution will be also worse  DeltaT2/M

   This variant is good if we do not need very good period resolution and pulsar is weak
   to be detected in the time series with time resolution DeltaT1.

   I have implemented this possibility already in ffasearch. The corresponding option
   is -m or --mbins

3) extra (or final) rebinning.
   we can use FFA with time resolution of DeltaT1 (without using 2a) or 2b))
   but finally when we will get folded profiles for periods between P0 and P0+1,
   we can rebin them to increase SNR.

   This case allows to keep good period resolution and to increase SNR as in the case 2) also.
   Of course, one can use all of 3 types of rebinning.

   The third case is should actually the case for RRAT 0848-43. I have tried to find it 
   using 2b) and it is seen in the folded profiles, but the period resolution was very bad (order
   of several milliseconds) and the value of period of the candidate is wrong.

   So, for this RRAT I guess we can try to use 2b) with reasonable _mbins_ number to have 
   period resolution we need, and after that to rebin the folded profiles to look for
   the candidates with high SNRs.

   I have implemented this possibility already in ffasearch. The corresponding option
   is -p or --pbins
