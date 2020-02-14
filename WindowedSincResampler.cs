using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace QDSS
{
    public class WindowedSincReSampler
    {

        //    Notes:
        //fmax should be less than half of fsr, and less than half of new_fsr(the reciprocal of the x step size).
        //Filter quality increases with a larger window width.The wider the window, the closer fmax can approach half of fsr or new_fsr.
        //Several operations inside the FOR loop can be pre-calculated.
        //There are more optimal windows than the von Hann window.
        //If the x step size is rational the same Window and Sinc values will be recalculated repeatedly.Therefore these values can either be cached, 
        // or pre-calculated and stored in a table (polyphase interpolation); or interpolated from a smaller pre-calculated table; or computed from a 
        //set of low-order polynomials fitted to each section or lobe between zero-crossings of the windowed Sinc(Farrow). (Performance optimization
        //is left as an exercise for the student). 


        // - QDSS Windowed-Sinc ReSampling subroutine in Basic
        //
        // - This function can also be used for interpolation of FFT results
        //
        // function parameters
        // : x      = new sample point location(relative to old indexes)
        //(e.g.every other integer for 0.5x decimation)
        // : indat  = original data array
        // : alim   = size of data array
        // : fmax   = low pass filter cutoff frequency
        // : fsr    = sample rate
        // : wnwdth = width of windowed Sinc used as the low pass filter
        //
        // resamp() returns a filtered new sample point
        public static double Resample(double x, List<double> indat, double fmax, double fsr, int wnwdth)
        {
            int alim = indat.Count;
            double fmaxDivFsr = fmax / fsr;
            double pi2 = 2 * Math.PI;
            double r_g = 2 * fmaxDivFsr; // Calc gain correction factor
            double r_y = 0;
            int ic = 0;
            int wnwdthHalf = wnwdth / 2;
            for (int i = -wnwdthHalf; i < wnwdthHalf; i++) // For 1 window width
            {
                // Calc input sample index
                int j = (int)(x + i);
                if ((j >= 0) && (j < alim))
                {
                    // calculate von Hann Window.Scale and calculate Sinc
                    double r_w = 0.5 - 0.5 * Math.Cos(pi2 * (0.5 + (j - x) / wnwdth));
                    double r_a = pi2 * (j - x) * fmaxDivFsr;
                    double r_snc = (r_a != 0) ? (Math.Sin(r_a) / r_a) : 1;
                    r_y += (r_g * r_w * r_snc * indat[j]);
                    ic++;
                }
            }

            return (ic > 0 ? r_y : double.NaN); // Return new filtered sample
        }


        public static List<double> ResamplePARALLEL(List<double> indat, double inX0, double inDX, double outX0, double outDX, int nOut)
        {
            double oldFSR = 1 / inDX;
            double newFSR = 1 / outDX;
            double ratioFreq = newFSR / oldFSR;
            double fmax = (ratioFreq >= 0.5) ? Math.Min(oldFSR, newFSR) * 0.5 : newFSR; // less than half for freq rates

            int wnwdth = 20;
            int wnwdthHalf = wnwdth / 2;
            int alim = indat.Count;
            double fmaxDivFsr = fmax / oldFSR;
            double pi2 = 2 * Math.PI;
            double r_g = 2 * fmaxDivFsr; // Calc gain correction factor


            List<double> outDat = new List<double>();
            for (int k = 0; k < nOut; k++) outDat.Add(0);

            Parallel.For(0, nOut, k =>
            {
                double x = outX0 + k * outDX;
                double relativeSample = (x - inX0) / inDX;

                {
                    int ic = 0;
                    double r_y = 0;

                    for (int i = -wnwdthHalf; i < wnwdthHalf; i++) // For 1 window width
                    {
                        // Calc input sample index
                        int j = (int)(relativeSample + i);
                        if ((j >= 0) && (j < alim))
                        {
                            // calculate von Hann Window.Scale and calculate Sinc
                            double r_w = 0.5 - 0.5 * Math.Cos(pi2 * (0.5 + (j - relativeSample) / wnwdth));
                            double r_a = pi2 * (j - relativeSample) * fmaxDivFsr;
                            double r_snc = (r_a != 0) ? (Math.Sin(r_a) / r_a) : 1;
                            r_y += (r_g * r_w * r_snc * indat[j]);
                            ic++;
                        }
                    }

                    outDat[k] = (ic > 0 ? r_y : double.NaN);
                }

            }); // Parallel.For



            return (outDat);
        }





        public static List<double> ResampleSEQUENTIAL(List<double> indat, double inX0, double inDX, double outX0, double outDX, int nOut)
        {
            double oldFSR = 1 / inDX;
            double newFSR = 1 / outDX;
            double ratioFreq = newFSR / oldFSR;
            double fmax = (ratioFreq >= 0.5) ? Math.Min(oldFSR, newFSR) * 0.5 : newFSR; // less than half for freq rates

            int wnwdth = 20;
            int wnwdthHalf = wnwdth / 2;
            int alim = indat.Count;
            double fmaxDivFsr = fmax / oldFSR;
            double pi2 = 2 * Math.PI;
            double r_g = 2 * fmaxDivFsr; // Calc gain correction factor


            List<double> outDat = new List<double>();
            for (int k = 0; k < nOut; k++)
            {
                double x = outX0 + k * outDX;
                double relativeSample = (x - inX0) / inDX;

                {
                    int ic = 0;
                    double r_y = 0;

                    for (int i = -wnwdthHalf; i < wnwdthHalf; i++) // For 1 window width
                    {
                        // Calc input sample index
                        int j = (int)(relativeSample + i);
                        if ((j >= 0) && (j < alim))
                        {
                            // calculate von Hann Window.Scale and calculate Sinc
                            double r_w = 0.5 - 0.5 * Math.Cos(pi2 * (0.5 + (j - relativeSample) / wnwdth));
                            double r_a = pi2 * (j - relativeSample) * fmaxDivFsr;
                            double r_snc = (r_a != 0) ? (Math.Sin(r_a) / r_a) : 1;
                            r_y += (r_g * r_w * r_snc * indat[j]);
                            ic++;
                        }
                    }

                    outDat.Add(ic > 0 ? r_y : double.NaN);
                }

            }



            return (outDat);
        }








        // - Ron Nicholson's QDSS ReSampler cookbook recipe
        //           QDSS = Quick, Dirty, Simple and Short
        //           Version 0.1b - 2007 - Aug - 01
        //           Copyright 2007 Ronald H.Nicholson Jr.
        // No warranties implied.  Error checking, optimization, and
        //           quality assessment of the "results" is left as an exercise
        //              for the student.
        //(consider this code Open Source under a BSD style license)
        //           IMHO.YMMV.http://www.nicholson.com/rhn/dsp.html




    }
}
