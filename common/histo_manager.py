import numpy as np
import math
import ROOT
from ROOT import TH1D, TH2D, TH3D, TMath, TH1F

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
    h1rebin.Sumw2();

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
def get_significance(h1sig, h1bkg):
    h1r = h1sig.Clone("h1r");
    h1r.Reset();

    n = h1r.GetNbinsX();
    for i in range(0,n):
        s = h1sig.GetBinContent(i+1);
        b = h1bkg.GetBinContent(i+1);
        s_err = h1sig.GetBinError(i+1);
        b_err = h1bkg.GetBinError(i+1);

        sig = 0.0;
        sig_err = 0.0;
        if s < 0.5 and b < 0.5:
            sig = 0.0;
            sig_err = 0.0;
        else:
            sig = s/math.sqrt(s + 2.*b);
            #sig_err = math.sqrt( pow(pow(s+2*b,-1/2.) - s/2.*pow(s+2.*b, -3/2) ,2) *pow(s_err,2) + pow( 2*s * pow(s+2*b,-3/2.) ,2) * pow(b_err,2)  );
            sig_err = 0.0;
        h1r.SetBinContent(i+1, sig);
        h1r.SetBinError(i+1, sig_err);
    return h1r;
#______________________________________________________________________
def get_corrected_bkg_simple(R, R_err, h1m_LSpp_same, h1m_LSnn_same):
    h1bkg = h1m_LSpp_same.Clone("h1bkg");
    h1bkg.Reset();
    #h1bkg.Sumw2();

    Nm = h1bkg.GetNbinsX();
    for im in range(0,Nm):
        lspp     = h1m_LSpp_same.GetBinContent(im+1);
        lspp_err = h1m_LSpp_same.GetBinError(im+1);
        lsnn     = h1m_LSnn_same.GetBinContent(im+1);
        lsnn_err = h1m_LSnn_same.GetBinError(im+1);

        bkg = 0;
        bkg_err = 0;
        if lspp * lsnn > 1e-6 :
            bkg = 2.0 * R * TMath.Sqrt(lspp * lsnn);
            bkg_err = TMath.Sqrt( R*R * ( pow(lspp*lsnn_err,2) + pow(lsnn*lspp_err,2) ) / (lspp*lsnn) );
        elif lspp + lsnn > 1e-6 : #arithmetic mean
            bkg = 2.0 * R * 0.5 * (lspp + lsnn);
            bkg_err = TMath.Sqrt( R*R * (lspp_err*lspp_err + lsnn_err*lsnn_err) );
        h1bkg.SetBinContent(im+1,bkg);
        h1bkg.SetBinError(im+1,bkg_err);

    return h1bkg;
#______________________________________________________________________
def get_R_factor(h1m_ULSnp_mix, h1m_ULSpn_mix, h1m_LSpp_mix, h1m_LSnn_mix):
    h1R   = h1m_ULSnp_mix.Clone("h1R");
    h1R.Reset();
    #h1R.Sumw2();

    h1m_ULS_mix = h1m_ULSnp_mix.Clone("h1m_ULS_mix");#sum of np+pn
    if h1m_ULSpn_mix is not None:
        h1m_ULS_mix.Add(h1m_ULSpn_mix,1.);

    Nm = h1R.GetNbinsX();
    for im in range(0,Nm):
        uls      = h1m_ULS_mix .GetBinContent(im+1);
        uls_err  = h1m_ULS_mix .GetBinError(im+1);
        lspp     = h1m_LSpp_mix.GetBinContent(im+1);
        lspp_err = h1m_LSpp_mix.GetBinError(im+1);
        lsnn     = h1m_LSnn_mix.GetBinContent(im+1);
        lsnn_err = h1m_LSnn_mix.GetBinError(im+1);
        R = 1.;
        R_err = 0.;

        if uls > 1e-6:
            if lspp * lsnn > 1e-6 :
                R = uls / ( 2 * TMath.Sqrt( lspp*lsnn ) );
                R_err = TMath.Sqrt( ( pow(lspp*lsnn_err*uls,2) + pow(lspp_err*lsnn*uls,2) + 4*pow(lspp*lsnn*uls_err,2) ) / ( 16 * pow(lspp*lsnn,3) ) );
            elif lspp + lsnn > 1e-6 :
                R = uls / ( 2 * 0.5 * (lspp + lsnn) );
                R_err = TMath.Sqrt( (pow(uls*lspp_err,2) + pow(uls*lsnn_err,2) + pow(lspp+lsnn,2)*pow(uls_err,2) ) / pow(lspp+lsnn,4) );
        else:
            R = 1.;
            R_err = 0.;
        #print("mix im+1 = {0} , x = {1} , R = {2} , uls = {3} , lspp = {4} , lsnn = {5}".format(im+1,h1R.GetBinCenter(im+1),R,uls,lspp,lsnn));
        h1R.SetBinContent(im+1, R);
        h1R.SetBinError(im+1, R_err);
    return h1R;
#______________________________________________________________________
def get_corrected_bkg(h1R, h1m_LSpp_same, h1m_LSnn_same):
    h1bkg = h1m_LSpp_same.Clone("h1bkg");
    h1bkg.Reset();
    #h1bkg.Sumw2();

    Nm = h1bkg.GetNbinsX();
    for im in range(0,Nm):
        lspp     = h1m_LSpp_same.GetBinContent(im+1);
        lspp_err = h1m_LSpp_same.GetBinError(im+1);
        lsnn     = h1m_LSnn_same.GetBinContent(im+1);
        lsnn_err = h1m_LSnn_same.GetBinError(im+1);
        R        = h1R.GetBinContent(im+1);
        R_err    = h1R.GetBinError(im+1);

        #print("im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsnn = {4}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsnn));

        bkg = 0;
        bkg_err = 0;
        if lspp * lsnn > 1e-6 : #geometric mean
            bkg = 2.0 * R * TMath.Sqrt(lspp * lsnn);
            bkg_err = TMath.Sqrt( pow(R_err/R,2) + 1./4*pow(lspp_err/lspp,2) + 1./4*pow(lsnn_err/lsnn,2) ) * bkg;
            #print("im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsnn = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsnn,bkg_err));
            #bkg_err = TMath.Sqrt( R*R * ( pow(lspp*lsnn_err,2) + pow(lsnn*lspp_err,2) ) / (lspp*lsnn) );
            #print("old im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsnn = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsnn,bkg_err));
        elif lspp + lsnn > 1e-6 : #arithmetic mean
            bkg = 2.0 * R * 0.5 * (lspp + lsnn);
            bkg_err = TMath.Sqrt( pow(R_err/R,2) + (pow(lspp_err,2)+pow(lsnn_err,2))/pow(lspp+lsnn,2) ) * bkg;
            #print("im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsnn = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsnn,bkg_err));
            #bkg_err = TMath.Sqrt( R*R * (lspp_err*lspp_err + lsnn_err*lsnn_err) );
            #print("old im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsnn = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsnn,bkg_err));
        h1bkg.SetBinContent(im+1,bkg);
        h1bkg.SetBinError(im+1,bkg_err);

    return h1bkg;

#______________________________________________________________________
#______________________________________________________________________
