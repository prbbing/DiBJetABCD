#!/usr/bin/python

import os
from ROOT import *
import re
import math
from array import array

def Round2DHistograms(histogram, precision):
  newHistogram = histogram.Clone()
  for j in range(1, histogram.GetYaxis().GetNbins() + 1):
    for i in range(1, histogram.GetXaxis().GetNbins() + 1): 
      newHistogram.SetBinContent(i, j, round(histogram.GetBinContent(i,j), precision))
      newHistogram.SetBinError(i, j, round(histogram.GetBinError(i,j), precision))
  return newHistogram


def rebin2D(hist,option):
  etaBins = array("f",[0,0.4,0.8,1.2,1.6,2.0])
  if option == "normal":
    ptBins = array("f",[0,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,1100,1200,1400,1600,5000])
  else:
    ptBins = array("f",[0,50,100,150,200,250,300,350,400,450,500,600,700,800,900,1000,1200,5000])
  newHist = TH2F(hist.GetTitle(),hist.GetTitle(),len(ptBins) -1 ,ptBins,len(etaBins) -1,etaBins)
  newHist.Sumw2() 
  for i in range(1, hist.GetXaxis().GetNbins()):
    for j in range(1,hist.GetYaxis().GetNbins()):
      newHist.Fill(hist.GetXaxis().GetBinCenter(i),abs(hist.GetYaxis().GetBinCenter(j)),hist.GetBinContent(i,j))
  for i in range(1, newHist.GetXaxis().GetNbins()):
    for j in range(1,newHist.GetYaxis().GetNbins()):
      newHist.SetBinError(i,j,math.sqrt(newHist.GetBinContent(i,j)))
  return newHist


def getRatioError(num,numerror,den,denerror):
  error = math.pow(math.pow(numerror,2)/math.pow(den,2) + math.pow(denerror,2)*math.pow(num,2)/math.pow(den,4),0.5)
  return error
#inputFile = TFile("/users/bingxuan.liu/DiBjet/2016Frozen/analysis_sys_dijet/out/wilkEff/data16_25ns_HLT_j380.root")
#inputFile = TFile("/users/bingxuan.liu/DiBjet/2016Frozen/analysis_sys_dijet_lowmass/out/wilkEff/data16_25ns_HLT_j380_PeriodI_5p89ifb_test.root")
inputFile = TFile("/users/bingxuan.liu/DiBjet/2016Frozen/analysis_sys_dijet/out/QCD_eff/sys0/merged.root")
outputFile = TFile("efficiency_mc_highmass.root","RECREATE")
outputFile.cd()
leadingJetHistNum = inputFile.Get("LeadingJetEtavsPt_b").Clone()
leadingJetHistNum.SetDirectory(0)
leadingJetHistNum.Sumw2()
secondJetHistNum = inputFile.Get("SecondJetEtavsPt_bj").Clone()
secondJetHistNum.SetDirectory(0)
secondJetHistNum.Sumw2()
leadingJetHistDen = inputFile.Get("LeadingJetEtavsPt").Clone()
leadingJetHistDen.SetDirectory(0)
leadingJetHistDen.Sumw2()
secondJetHistDen = inputFile.Get("SecondJetEtavsPt").Clone()
secondJetHistDen.SetDirectory(0)
secondJetHistDen.Sumw2()
secondJetHistDen_con = inputFile.Get("SecondJetEtavsPt_b").Clone()
secondJetHistDen_con.SetDirectory(0)
secondJetHistDen_con.Sumw2()
secondJetHistNum_con = inputFile.Get("SecondJetEtavsPt_bb").Clone()
secondJetHistNum_con.SetDirectory(0)
secondJetHistNum_con.Sumw2()

JetHistNum = inputFile.Get("LeadingJetEtavsPt_b").Clone()
#JetHistNum.Add(secondJetHistNum)
JetHistDen = inputFile.Get("LeadingJetEtavsPt").Clone()
#JetHistDen.Add(secondJetHistDen)
JetHistNum = rebin2D(JetHistNum,"normal")
JetHistDen = rebin2D(JetHistDen,"normal")
JetHistNum.Divide(JetHistNum ,JetHistDen,1,1,"B")
JetHistNum.SetName("JetEff")
JetHistNum.SetTitle("JetEff")

leadingJetHistNum = rebin2D(leadingJetHistNum,"normal")
leadingJetHistDen = rebin2D(leadingJetHistDen,"normal")
secondJetHistNum_con = rebin2D(secondJetHistNum_con,"abnormal")
secondJetHistDen_con = rebin2D(secondJetHistDen_con,"abnormal")
secondJetHistNum_con.Divide(secondJetHistNum_con,secondJetHistDen_con,1,1,"B")
secondJetHistNum_con.SetName("secondJetEff_con")
secondJetHistNum_con.SetTitle("secondJetEff_con")
leadingJetHistNum.Divide(leadingJetHistNum ,leadingJetHistDen,1,1,"B")
leadingJetHistNum.SetName("leadingJetEff")
leadingJetHistNum.SetTitle("leadingJetEff")
secondJetHistNum = rebin2D(secondJetHistNum,"normal")
secondJetHistDen = rebin2D(secondJetHistDen,"normal")
secondJetHistNum.Divide(secondJetHistNum,secondJetHistDen,1,1,"B")
secondJetHistNum.SetName("secondJetEff")
secondJetHistNum.SetTitle("secondJetEff")

leadingJetHistNum.Write()
secondJetHistNum.Write()
secondJetHistNum_con.Write()
JetHistNum.Write()
leadingJetCanvas = TCanvas("leadingJetEff","leadingJetEff",10,10,600,600)
leadingJetCanvas.Divide(1,1,0.008,0.007)
leadingJetCanvas.cd(1)
leadingJetHistNum = Round2DHistograms(leadingJetHistNum,3)
leadingJetHistNum.Draw("TEXTE COLZ")
secondJetCanvas = TCanvas("secondJetEff","secondJetEff",10,10,600,600)
secondJetCanvas.Divide(1,1,0.008,0.007)
secondJetCanvas.cd(1)
secondJetHistNum=Round2DHistograms(secondJetHistNum,3)
secondJetHistNum.Draw("TEXTE COLZ")
JetCanvas = TCanvas("JetEff","JetEff",10,10,600,600)
JetCanvas.Divide(1,1,0.008,0.007)
JetCanvas.cd(1)
JetHistNum=Round2DHistograms(JetHistNum,3)
JetHistNum.Draw("TEXTE COLZ")
secondJet_con_Canvas = TCanvas("secondJetEff_con","secondJetEff_con",10,10,600,600)
secondJet_con_Canvas.Divide(1,1,0.008,0.007)
secondJet_con_Canvas.cd(1)
secondJetHistNum_con=Round2DHistograms(secondJetHistNum_con,3)
secondJetHistNum_con.Draw("TEXTE COLZ")
#leadingJetCanvas.Write()
#secondJetCanvas.Write()
#secondJet_con_Canvas.Write()
#JetCanvas.Write()

outputFile.Close()
