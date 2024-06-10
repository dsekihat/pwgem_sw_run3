import numpy as np
import math
import ROOT
from ROOT import TH1D, TH2D, TH3D, TMath, TH1F

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
def get_R_factor(h1m_ULSnp_mix, h1m_ULSpn_mix, h1m_lspp_mix, h1m_lsmm_mix):
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
        lspp     = h1m_lspp_mix.GetBinContent(im+1);
        lspp_err = h1m_lspp_mix.GetBinError(im+1);
        lsmm     = h1m_lsmm_mix.GetBinContent(im+1);
        lsmm_err = h1m_lsmm_mix.GetBinError(im+1);
        R = 1.;
        R_err = 0.;

        if uls > 1e-6:
            if lspp * lsmm > 1e-6 :
                R = uls / ( 2 * TMath.Sqrt( lspp*lsmm ) );
                R_err = TMath.Sqrt( ( pow(lspp*lsmm_err*uls,2) + pow(lspp_err*lsmm*uls,2) + 4*pow(lspp*lsmm*uls_err,2) ) / ( 16 * pow(lspp*lsmm,3) ) );
            elif lspp + lsmm > 1e-6 :
                R = uls / ( 2 * 0.5 * (lspp + lsmm) );
                R_err = TMath.Sqrt( (pow(uls*lspp_err,2) + pow(uls*lsmm_err,2) + pow(lspp+lsmm,2)*pow(uls_err,2) ) / pow(lspp+lsmm,4) );
        else:
            R = 1.;
            R_err = 0.;
        #print("mix im+1 = {0} , x = {1} , R = {2} , uls = {3} , lspp = {4} , lsmm = {5}".format(im+1,h1R.GetBinCenter(im+1),R,uls,lspp,lsmm));
        h1R.SetBinContent(im+1, R);
        h1R.SetBinError(im+1, R_err);
    return h1R;
#______________________________________________________________________
def get_corrected_bkg(h1R, h1m_lspp_same, h1m_lsmm_same):
    h1bkg = h1m_lspp_same.Clone("h1bkg");
    h1bkg.Reset();
    #h1bkg.Sumw2();

    Nm = h1bkg.GetNbinsX();
    for im in range(0,Nm):
        lspp     = h1m_lspp_same.GetBinContent(im+1);
        lspp_err = h1m_lspp_same.GetBinError(im+1);
        lsmm     = h1m_lsmm_same.GetBinContent(im+1);
        lsmm_err = h1m_lsmm_same.GetBinError(im+1);
        R        = h1R.GetBinContent(im+1);
        R_err    = h1R.GetBinError(im+1);

        #print("im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsmm = {4}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsmm));

        bkg = 0;
        bkg_err = 0;
        if lspp * lsmm > 1e-6 : #geometric mean
            bkg = 2.0 * R * TMath.Sqrt(lspp * lsmm);
            bkg_err = TMath.Sqrt( pow(R_err/R,2) + 1./4*pow(lspp_err/lspp,2) + 1./4*pow(lsmm_err/lsmm,2) ) * bkg;
            #print("im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsmm = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsmm,bkg_err));
            #bkg_err = TMath.Sqrt( R*R * ( pow(lspp*lsmm_err,2) + pow(lsmm*lspp_err,2) ) / (lspp*lsmm) );
            #print("old im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsmm = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsmm,bkg_err));
        elif lspp + lsmm > 1e-6 : #arithmetic mean
            bkg = 2.0 * R * 0.5 * (lspp + lsmm);
            bkg_err = TMath.Sqrt( pow(R_err/R,2) + (pow(lspp_err,2)+pow(lsmm_err,2))/pow(lspp+lsmm,2) ) * bkg;
            #print("im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsmm = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsmm,bkg_err));
            #bkg_err = TMath.Sqrt( R*R * (lspp_err*lspp_err + lsmm_err*lsmm_err) );
            #print("old im+1 = {0} , x = {1} , R = {2} , lspp = {3} , lsmm = {4}, bkg_err = {5}".format(im+1,h1R.GetBinCenter(im+1),R,lspp,lsmm,bkg_err));
        h1bkg.SetBinContent(im+1,bkg);
        h1bkg.SetBinError(im+1,bkg_err);

    return h1bkg;
#______________________________________________________________________
def get_corrected_bkg_simple(R, R_err, h1m_lspp_same, h1m_lsmm_same):
    h1bkg = h1m_lspp_same.Clone("h1bkg");
    h1bkg.Reset();
    #h1bkg.Sumw2();

    Nm = h1bkg.GetNbinsX();
    for im in range(0,Nm):
        lspp     = h1m_lspp_same.GetBinContent(im+1);
        lspp_err = h1m_lspp_same.GetBinError(im+1);
        lsmm     = h1m_lsmm_same.GetBinContent(im+1);
        lsmm_err = h1m_lsmm_same.GetBinError(im+1);

        bkg = 0;
        bkg_err = 0;
        if lspp * lsmm > 1e-6 :
            bkg = 2.0 * R * TMath.Sqrt(lspp * lsmm);
            bkg_err = TMath.Sqrt( R*R * ( pow(lspp*lsmm_err,2) + pow(lsmm*lspp_err,2) ) / (lspp*lsmm) );
        elif lspp + lsmm > 1e-6 : #arithmetic mean
            bkg = 2.0 * R * 0.5 * (lspp + lsmm);
            bkg_err = TMath.Sqrt( R*R * (lspp_err*lspp_err + lsmm_err*lsmm_err) );
        h1bkg.SetBinContent(im+1,bkg);
        h1bkg.SetBinError(im+1,bkg_err);

    return h1bkg;
#______________________________________________________________________
