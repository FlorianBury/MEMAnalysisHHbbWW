import os
import sys
import copy
import json
import argparse
import itertools
import numpy as np
import pandas as pd
import shutil
import ROOT
from IPython import embed
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import ndimage

from data_helper import Observable, Data
from numpy_hist import NumpyHist

EPSILON = 1e-8

class JetWeight:
    def __init__(self,jetPath,outputDir):
        # Make output directory #
        self.outputDir = outputDir
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        # Open TFile #
        jetFile = ROOT.TFile(jetPath)
        # Produce ratio histograms #
        self.nph_ak4_bjet = self._makeRatio(F        = jetFile,
                                  numNames = ['N_ak4_truth_bjets'],
                                  denNames = ['N_ak4_truth_bjets',
                                              'N_ak4_truth_cjets',
                                              'N_ak4_truth_lightjets'])
        self.nph_ak4_cjet = self._makeRatio(F        = jetFile,
                                  numNames = ['N_ak4_truth_cjets'],
                                  denNames = ['N_ak4_truth_bjets',
                                              'N_ak4_truth_cjets',
                                              'N_ak4_truth_lightjets'])
        self.nph_ak4_ljet = self._makeRatio(F        = jetFile,
                                  numNames = ['N_ak4_truth_lightjets'],
                                  denNames = ['N_ak4_truth_bjets',
                                              'N_ak4_truth_cjets',
                                              'N_ak4_truth_lightjets'])
        self.nph_ak8_bjet = self._makeRatio(F        = jetFile,
                                  numNames = ['N_ak8_truth_subJet1_bjets',
                                              'N_ak8_truth_subJet2_bjets'],
                                  denNames = ['N_ak8_truth_subJet1_bjets',
                                              'N_ak8_truth_subJet2_bjets',
                                              'N_ak8_truth_subJet1_cjets',
                                              'N_ak8_truth_subJet2_cjets',
                                              'N_ak8_truth_subJet1_lightjets',
                                              'N_ak8_truth_subJet2_lightjets'])
        self.nph_ak8_cjet = self._makeRatio(F        = jetFile,
                                  numNames = ['N_ak8_truth_subJet1_cjets',
                                              'N_ak8_truth_subJet2_cjets'],
                                  denNames = ['N_ak8_truth_subJet1_bjets',
                                              'N_ak8_truth_subJet2_bjets',
                                              'N_ak8_truth_subJet1_cjets',
                                              'N_ak8_truth_subJet2_cjets',
                                              'N_ak8_truth_subJet1_lightjets',
                                              'N_ak8_truth_subJet2_lightjets'])
        self.nph_ak8_ljet = self._makeRatio(F        = jetFile,
                                  numNames = ['N_ak8_truth_subJet1_lightjets',
                                              'N_ak8_truth_subJet2_lightjets'],
                                  denNames = ['N_ak8_truth_subJet1_bjets',
                                              'N_ak8_truth_subJet2_bjets',
                                              'N_ak8_truth_subJet1_cjets',
                                              'N_ak8_truth_subJet2_cjets',
                                              'N_ak8_truth_subJet1_lightjets',
                                              'N_ak8_truth_subJet2_lightjets'])
        jetFile.Close()
        # Rebin #
        rebin_ak4 = [
            np.array([0,1.2,2.0,2.5]),
            np.array([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1040, 1060, 1080, 1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280, 1300, 1320, 1340, 1360, 1380, 1400, 1420, 1440, 1460, 1480, 1500, 1520, 1540, 1560, 1580, 1600, 1620, 1640, 1660, 1680, 1700, 1720, 1740, 1760, 1780, 1800, 1820, 1840, 1860, 1880, 1900, 1920, 1940, 1960, 1980, 2000, 2020, 2040, 2060, 2080, 2100, 2120, 2140, 2160, 2180, 2200, 2220, 2240, 2260, 2280, 2300, 2320, 2340, 2360, 2380, 2400, 2420, 2440, 2460, 2480, 2500, 2520, 2540, 2560, 2580, 2600, 2620, 2640, 2660, 2680, 2700, 2720, 2740, 2760, 2780, 2800, 2820, 2840, 2860, 2880, 2900, 2920, 2940, 2960, 2980, 3000]),
            np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0])
        ]
        rebin_ak8 = [
            np.array([0,1.2,2.0,2.5]),
            np.array([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500, 2550, 2600, 2650, 2700, 2750, 2800, 2850, 2900, 2950, 3000]),
            np.array([0.0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0])
        ]

        self.nph_ak4_bjet = self.nph_ak4_bjet.rebin(rebin_ak4,average_bins=True)
        self.nph_ak4_cjet = self.nph_ak4_cjet.rebin(rebin_ak4,average_bins=True)
        self.nph_ak4_ljet = self.nph_ak4_ljet.rebin(rebin_ak4,average_bins=True)
        self.nph_ak8_bjet = self.nph_ak8_bjet.rebin(rebin_ak8,average_bins=True)
        self.nph_ak8_cjet = self.nph_ak8_cjet.rebin(rebin_ak8,average_bins=True)
        self.nph_ak8_ljet = self.nph_ak8_ljet.rebin(rebin_ak8,average_bins=True)
        # Correct histograms #
        self._correctHist(self.nph_ak4_bjet)
        self._correctHist(self.nph_ak4_cjet)
        self._correctHist(self.nph_ak4_ljet)
        self._correctHist(self.nph_ak8_bjet)
        self._correctHist(self.nph_ak8_cjet)
        self._correctHist(self.nph_ak8_ljet)
        # Make observables #
        self.data_ak4_bjet = self._makeObservable(self.nph_ak4_bjet)
        self.data_ak4_cjet = self._makeObservable(self.nph_ak4_cjet)
        self.data_ak4_ljet = self._makeObservable(self.nph_ak4_ljet)
        self.data_ak8_bjet = self._makeObservable(self.nph_ak8_bjet)
        self.data_ak8_cjet = self._makeObservable(self.nph_ak8_cjet)
        self.data_ak8_ljet = self._makeObservable(self.nph_ak8_ljet)

    def _makeRatio(self,F,numNames,denNames):
        # Make ratio #
        h_num = self._getHist(F,numNames)
        h_den = self._getHist(F,denNames) + EPSILON
        h_ratio = h_num / h_den
        # Correct weird values #
        h_ratio._w[h_ratio._w>1] = 1.
        h_ratio._w[h_ratio._w<0.] = 0.
        # Rounding issue #
        h_ratio._e = [h_ratio._e[0].round(3),h_ratio._e[1].round(3),h_ratio._e[2].round(3)]
        return h_ratio

    @staticmethod
    def _getHist(F,names):
        hsum = None
        for name in names:
            h = F.Get(name)
            nph = NumpyHist.getFromRoot(h)
            if hsum is None:
                hsum = nph
            else:
                hsum.add(nph)
        return hsum

    @staticmethod
    def _correctHist(nph):
        for i in range(nph.w.shape[0]):
            nph.w[i,:,:] = ndimage.uniform_filter(nph.w[i,:,:])

    @staticmethod
    def _makeObservable(nph):
        centers = [((e[1:]+e[:-1])/2).round(3) for e in nph.e]
        return Observable(nph.w,{lname:lval.round(3) for lname,lval in zip(['eta','pt','score'],centers)})

    def makeROOT(self):
        subdir = os.path.join(self.outputDir,'root')
        if not os.path.exists(subdir):
            os.makedirs(subdir)
        h_ak4_bjet = self.nph_ak4_bjet.fillHistogram('ak4_bjet_prob')
        h_ak4_cjet = self.nph_ak4_cjet.fillHistogram('ak4_cjet_prob')
        h_ak4_ljet = self.nph_ak4_ljet.fillHistogram('ak4_ljet_prob')
        h_ak8_bjet = self.nph_ak8_bjet.fillHistogram('ak8_bjet_prob')
        h_ak8_cjet = self.nph_ak8_cjet.fillHistogram('ak8_cjet_prob')
        h_ak8_ljet = self.nph_ak8_ljet.fillHistogram('ak8_ljet_prob')
    
        filename = os.path.join(subdir,'jet_prob.root')
        if os.path.exists(filename):
            shutil.rmtree(filename)

        F = ROOT.TFile(filename,"RECREATE")

        h_ak4_bjet.Write() 
        h_ak4_cjet.Write()
        h_ak4_ljet.Write()
        h_ak8_bjet.Write()
        h_ak8_cjet.Write()
        h_ak8_ljet.Write()

        F.Write()
        F.Close()

        print (f'Saved {filename}')

    def makePlots(self):
        self._makePlots(self.data_ak4_bjet,suffix='AK4_bjet',label='AK4 b-jet probability')
        self._makePlots(self.data_ak4_cjet,suffix='AK4_cjet',label='AK4 c-jet probability')
        self._makePlots(self.data_ak4_ljet,suffix='AK4_ljet',label='AK4 light-jet probability')
        self._makePlots(self.data_ak8_bjet,suffix='AK8_bjet',label='AK8 b-subjet probability')
        self._makePlots(self.data_ak8_cjet,suffix='AK8_cjet',label='AK8 c-subjet probability')
        self._makePlots(self.data_ak8_ljet,suffix='AK8_ljet',label='AK8 light-subjet probability')

    def _makePlots(self,observable,suffix,label):
        from scipy import ndimage, misc
        axes_names = {'eta':'Jet $\eta$','pt':'Jet $p_{T}$','score':'Jet btag score'}
        pdfDir = os.path.join(self.outputDir,'pdf')
        #observable.array = ndimage.uniform_filter(observable.array,size=3,mode='nearest')
        if not os.path.exists(pdfDir):
            os.makedirs(pdfDir)
        for primary_name in ['eta']:
            print (f'Looking at {primary_name}')
            for primary_value in observable.GetLabels()[primary_name]:
                title = f'{axes_names[primary_name]} = {primary_value:.3f}'
                obs2d = observable.GetSlice(**{primary_name:primary_value})
                # Plot 1D #
                for secondary_label in obs2d.GetLabels().keys():
                    tertiary_label = list(set(obs2d.GetLabels().keys())-set([secondary_label]))[0]
                    fig, ax = plt.subplots(figsize=(8,7))
                    plt.subplots_adjust(left=0.15, right=0.90, top=0.85, bottom=0.12)
                    secondary_values = obs2d.GetLabels()[secondary_label]
                    colors = cm.jet(np.linspace(0,1,secondary_values.shape[0]))
                    for secondary_val,color in zip(secondary_values,colors):
                        obs1d = obs2d.GetSlice(**{secondary_label:secondary_val})
                        obs1d.Pyplot1D(ax,color=color)
                    sm = cm.ScalarMappable(cmap=cm.rainbow, norm=plt.Normalize(vmin=secondary_values.min(), vmax=secondary_values.max()))
                    cbar = plt.colorbar(sm)
                    cbar.set_label(axes_names[secondary_label],fontsize=18,labelpad=20)
                    cbar.ax.tick_params(labelsize=14)
                    plt.xlabel(axes_names[tertiary_label],fontsize=18,labelpad=20)
                    plt.ylabel(label,fontsize=18,labelpad=10)
                    plt.ylim(0,1)
                    plt.title(title,fontsize=20,pad=25)
                    plt.tick_params(axis='both', which='major', labelsize=14)
                    pngName = os.path.join(pdfDir,f'{suffix}_1D_{primary_name}_{primary_value:.3f}_axis_{tertiary_label}').replace(' ','_').replace('$','').replace('.','p')+'.png'
                    fig.savefig(pngName)
                    print (f'\tProduced {pngName}')
                    plt.close()

                # Plot 2D #
                fig,ax = plt.subplots(figsize=(8,7))
                other_axes = [key for key in observable.GetLabels().keys() if key != primary_name]
                plt.subplots_adjust(left=0.15, right=0.95, top=0.85, bottom=0.1)
                #obs2d.Pyplot2D(x=other_axes[0],y=other_axes[1], ax=ax, shading='auto', linewidth=0,rasterized=True, vmin=0,vmax=1)
                obs2d.Pyplot2D(x=other_axes[0],y=other_axes[1], ax=ax, shading='auto', linewidth=0,rasterized=True)
                cbar = fig.colorbar(ax.collections[0])
                cbar.set_label(label,fontsize=18,labelpad=20)
                cbar.ax.tick_params(labelsize=14)
                ax.set_xlabel(axes_names[other_axes[0]],fontsize=18,labelpad=10)
                ax.set_ylabel(axes_names[other_axes[1]],fontsize=18,labelpad=20)
                ax.set_title(title,fontsize=20,pad=25)
                ax.tick_params(axis='both', which='major', labelsize=14)
                pngName = os.path.join(pdfDir,f'{suffix}_2D_{primary_name}_{primary_value:.3f}').replace(' ','_').replace('$','').replace('.','p')+'.png'
                fig.savefig(pngName)
                print (f'\tProduced {pngName}')
                plt.close()
    
        
        
parser = argparse.ArgumentParser(description='Produce b prob files and plots')
parser.add_argument('--data', action='store', required=True, type=str, 
                    help='Path to root file')
parser.add_argument('--out', action='store', required=True, type=str, 
                    help='Output dir')
parser.add_argument('--plots', action='store_true', required=False, default=False,
                    help='Produce plots')
parser.add_argument('--root', action='store_true', required=False, default=False,
                    help='Produce root file')
#parser.add_argument('--json', action='store_true', required=False, default=False,
#                    help='Produce json')
#parser.add_argument('--root', action='store_true', required=False, default=False,
#                    help='Produce root')    
args = parser.parse_args()

instance = JetWeight(args.data,args.out)
if args.plots:
    instance.makePlots()
if args.root:
    instance.makeROOT()

