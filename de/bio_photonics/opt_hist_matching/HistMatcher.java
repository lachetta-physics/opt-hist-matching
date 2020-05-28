/*
 * Copyright (C) 2019 m.lachetta
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.bio_photonics.opt_hist_matching;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

/**
 *
 * @author m.lachetta
 */
public class HistMatcher implements PlugIn {
    static int NBINS = 128;
    @Override
    public void run(String arg) {
        ImagePlus ip = ij.WindowManager.getCurrentImage();
        ImageStack is = ip.getImageStack();
        ImageStack newIs = new ImageStack(is.getWidth(), is.getHeight(), is.getSize());
        for (int slice = 1; slice <= is.size(); slice++) {
            ImageProcessor processor = is.getProcessor(slice);
            double min = processor.getMin();
            double max = processor.getMax();
            float binWidth = (float) (max / NBINS);
            int[] histogram = processor.getHistogram(NBINS);
            
            //find histogram peak
            int peakIndex = -1;
            int peakHeight = Integer.MIN_VALUE;
            for (int i = 0; i < NBINS; i++) {
                if (histogram[i] > peakHeight) {
                    peakIndex = i;
                    peakHeight = histogram[i];
                }
            }
            
            //find lower border
            int lowerBorder = 0;
            int upperBin = peakIndex;
            int lowerBin = upperBin - 1;
            int lowerCounter = 0;
            if (lowerBin >= 0) {
                int upperHeight = histogram[upperBin];
                int lowerHeight = histogram[lowerBin];
                while (lowerHeight <= upperHeight && lowerBin > 0) {
                    upperBin--;
                    lowerBin--;
                    upperHeight = histogram[upperBin];
                    lowerHeight = histogram[lowerBin];
                    lowerCounter++;
                }
                lowerBorder = upperBin;
            }
            
            //find upper border
            int upperBorder = 0;
            upperBin = peakIndex + 1;
            lowerBin = peakIndex;
            int upperCounter = 0;
            if (upperBin <= NBINS) {
                int upperHeight = histogram[upperBin];
                int lowerHeight = histogram[lowerBin];
                while (lowerHeight >= upperHeight && upperBin < NBINS-1) {
                    upperBin++;
                    lowerBin++;
                    upperHeight = histogram[upperBin];
                    lowerHeight = histogram[lowerBin];
                    upperCounter++;
                }
                upperBorder = upperBin;
            }
            
            //gauss fitting
            IJ.log("Fitting slice " + slice);
            GaussianCurveFitter fitter = GaussianCurveFitter.create();

            WeightedObservedPoints obs = new WeightedObservedPoints();

            for (int index = lowerBorder; index < upperBin; index++) {
                obs.add(index, histogram[index]);
            }

            double[] bestFit = fitter.fit(obs.toList());
            double gaussHeight = bestFit[0];
            double gaussMean = bestFit[1];
            double gaussSigma = bestFit[2];
            
            
            // normalize with gauss parameters
            float[][] floatArray = processor.getFloatArray();
            for (int i = 0; i < floatArray.length; i++) for (int k = 0; k < floatArray[i].length; k++) {
                floatArray[i][k] -= gaussMean * binWidth;
                floatArray[i][k] = (float) (floatArray[i][k] / gaussSigma);
            }
            FloatProcessor fp = new FloatProcessor(floatArray);
            newIs.setProcessor(fp, slice);
        }
        ImagePlus newIp = new ImagePlus(ip.getTitle() + "_hm", newIs);
        newIp.show();
    }
        
    /** main method for testing*/
    public static void main( String [] args ) {
	HistMatcher a  = new HistMatcher();
        ImageJ imageJ = new ij.ImageJ( ij.ImageJ.EMBEDDED);
        ImagePlus ip = new ImagePlus("D:\\OPT-images\\wt-hirn-red\\binnedx4-Stack Reconstruction.tif");
        ip.show();
	a.run("");
    }

}
