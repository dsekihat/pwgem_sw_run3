import numpy as np
import math
import ROOT
from ROOT import TH1D, TH2D, TH3D, TMath, TH1F, TProfile

#______________________________________________________________________
def convert_dn2n(h1org):
    #simply convert differential yield to the number of counts. i.e. dn/dx -> n as a function of x
    h1 = h1org.Clone("h1");
    h1.Reset();
    for i in range(0,h1.GetNbinsX()):
        dndx     = h1org.GetBinContent(i+1);
        dndx_err = h1org.GetBinError(i+1);
        dx       = h1org.GetBinWidth(i+1);
        h1.SetBinContent(i+1,dndx * dx);
        h1.SetBinError(i+1,dndx_err * dx);
    return h1;
#______________________________________________________________________
def rebin_histogram_2d(h2org, arrX, arrY):
    #please make sure that bin edges coincide between old and new histograms

    nx = len(arrX);
    ny = len(arrY);
    h2 = TH2D("{0}_rebin".format(h2org.GetName()), h2org.GetTitle(), nx-1, arrX, ny-1, arrY);

    h2.SetXTitle(h2org.GetXaxis().GetTitle());
    h2.SetYTitle(h2org.GetYaxis().GetTitle());

    for ix in range(0, nx-1):
        x1 = arrX[ix];
        x2 = arrX[ix+1];
        bin_x1org = h2org.GetXaxis().FindBin(x1 + 1e-6);
        bin_x2org = h2org.GetXaxis().FindBin(x2 - 1e-6);
        for iy in range(0, ny-1):
            y1 = arrY[iy];
            y2 = arrY[iy+1];

            bin_y1org = h2org.GetYaxis().FindBin(y1 + 1e-6);
            bin_y2org = h2org.GetYaxis().FindBin(y2 - 1e-6);

            n = 0;
            n_err = 0;
            for ix_org in range(bin_x1org, bin_x2org+1):
                for iy_org in range(bin_y1org, bin_y2org+1):
                    n += h2org.GetBinContent(ix_org, iy_org);
                    n_err += math.pow(h2org.GetBinError(ix_org, iy_org), 2); #only stat. err is supported.

            h2.SetBinContent(ix+1, iy+1, n);
            h2.SetBinError(ix+1, iy+1, math.sqrt(n_err));

    return h2;
#______________________________________________________________________
def rebin_histogram(h1, arrX, isdiff, is_pt_correlated = False):
    h1tmp = h1.Clone("h1tmp");
    h1rebin = h1tmp.Rebin(len(arrX)-1,"h1rebin",arrX);
    h1rebin.SetName("{0}_rebin".format(h1.GetName()));
    #h1rebin.Sumw2();

    if is_pt_correlated:
        for i in range(0,h1rebin.GetNbinsX()):
            x_min = h1rebin.GetBinLowEdge(i+1);
            x_max = h1rebin.GetBinLowEdge(i+2);
            bin0 = h1.FindBin(x_min + 1e-6);
            bin1 = h1.FindBin(x_max - 1e-6);
            y = 0;
            err = 0;
            for j in range(bin0, bin1+1):
                y += h1.GetBinContent(j);
                err += h1.GetBinError(j) / h1.GetBinContent(j) * h1.GetBinContent(j);
            h1rebin.SetBinContent(i+1, y);
            h1rebin.SetBinError(i+1, err);
            #print("check 0 rel syst. = ", h1rebin.GetBinError(i+1) / h1rebin.GetBinContent(i+1));

    if isdiff and (h1rebin.Class() == TH1D.Class() or h1rebin.Class() == TH1F.Class() ) : #do you want differential histogram? e.g. dN/dpT, dN/dm.
        h1rebin.Scale(1.,"width");
    return h1rebin;
#______________________________________________________________________
def rebin_profile(h1, arr):
    h1tmp = h1.Clone("h1tmp");
    h1rebin = h1tmp.Rebin(len(arr)-1, "h1rebin", arr);
    return h1rebin;
#______________________________________________________________________
def slice_histogram(h2,x0,x1,axis,isdiff):
    h1 = 0;
    delta = 1e-6;
    if "x" in axis.lower():
        bin0 = h2.GetYaxis().FindBin(x0 + delta);
        bin1 = h2.GetYaxis().FindBin(x1 - delta);
        h1 = h2.ProjectionX("h1prjx_{0}".format(h2.GetName()),bin0,bin1,"");
    elif "y" in axis.lower():
        bin0 = h2.GetXaxis().FindBin(x0 + delta);
        bin1 = h2.GetXaxis().FindBin(x1 - delta);
        h1 = h2.ProjectionY("h1prjy_{0}".format(h2.GetName()),bin0,bin1,"");

    if isdiff and h1.Class() == TH1D.Class(): #do you want differential histogram? e.g. dN/dpT, dN/dm.
        h1.Scale(1.,"width");
    return h1;
#______________________________________________________________________
def slice_profile(h2,x0,x1,axis,isdiff=False):
    h1 = 0;
    delta = 1e-6;
    if "x" in axis.lower():
        bin0 = h2.GetYaxis().FindBin(x0 + delta);
        bin1 = h2.GetYaxis().FindBin(x1 - delta);
        h1 = h2.ProfileX("h1prfx_{0}".format(h2.GetName()),bin0,bin1,"");
    elif "y" in axis.lower():
        bin0 = h2.GetXaxis().FindBin(x0 + delta);
        bin1 = h2.GetXaxis().FindBin(x1 - delta);
        h1 = h2.ProfileY("h1prfy_{0}".format(h2.GetName()),bin0,bin1,"");

    if isdiff and h1.Class() == TProfile.Class(): #do you want differential profile?
        h1.Scale(1.,"width");
    return h1;
#______________________________________________________________________
def get_bkg_subtracted(h1all, h1bkg):
    h1sig = h1all.Clone("h1sig");
    #h1sig.Sumw2();
    h1sig.Add(h1bkg,-1);
    return h1sig;
#______________________________________________________________________
def get_ratio(h1sig, h1bkg, option="B"):
    h1r = h1sig.Clone("h1r");
    h1r.Reset();
    h1r.Divide(h1sig,h1bkg,1.,1., option);
    return h1r;
#______________________________________________________________________
def get_cumulative_histogram(h1org, isdiff, is_x_correlated = False):
    h1 = h1org.Clone("h1");
    h1.Reset();
    h1.Sumw2();
    n_integral = 0;
    n_integral_err = 0;
    for i in range(1, h1org.GetNbinsX() + 1):
        width = h1org.GetBinWidth(i);
        n = h1org.GetBinContent(i);
        n_err = h1org.GetBinError(i);
        if isdiff:
            n_integral += n * width;
            n_integral_err += math.pow(n_err * width, 2);
            h1.SetBinContent(i, n_integral);
            h1.SetBinError(i, math.sqrt(n_integral_err));
        else:
            n_integral += n;
            n_integral_err += math.pow(n_err, 2);
            h1.SetBinContent(i, n_integral);
            h1.SetBinError(i, n_integral_err);
    return h1;
#______________________________________________________________________
