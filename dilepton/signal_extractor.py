import numpy as np
import math
import ROOT
from ROOT import TH1D, TH2D, TH3D, TMath, TH1F, THnSparseF, THnSparseD

#__________________________________________________________________
def get_Rn(h1prf_sp12, h1prf_sp13, h1prf_sp23):
    h1Rn = h1prf_sp12.Clone("h1prf_Rn").ProjectionX("h1Rn");
    sp12_detname = h1prf_sp12.GetYaxis().GetTitle();
    sp13_detname = h1prf_sp13.GetYaxis().GetTitle();
    sp23_detname = h1prf_sp23.GetYaxis().GetTitle();
    R2_title = "R_{{2}} =  #sqrt{{#frac{{<{0}> <{1}>}}{{<{2}>}}}}".format(sp12_detname, sp13_detname, sp23_detname);
    h1Rn.GetYaxis().SetTitle(R2_title);

    nbin = h1Rn.GetNbinsX();
    for i in range(0, nbin):
        sp12 = h1prf_sp12.GetBinContent(i+1);
        sp13 = h1prf_sp13.GetBinContent(i+1);
        sp23 = h1prf_sp23.GetBinContent(i+1);
        sp12_err = h1prf_sp12.GetBinError(i+1);
        sp13_err = h1prf_sp13.GetBinError(i+1);
        sp23_err = h1prf_sp23.GetBinError(i+1);

        if sp23 > 1e-6:
            R2 = math.sqrt(sp12 * sp23 / sp23);
            R2_err = R2 * math.sqrt(math.pow(sp12_err/sp12, 2) + math.pow(sp13_err/sp13, 2) + math.pow(sp23_err/sp23, 2));
            h1Rn.SetBinContent(i+1, R2);
            h1Rn.SetBinError(i+1, R2_err);
        else:
            h1Rn.SetBinContent(i+1, 0);
            h1Rn.SetBinError(i+1, 0);
    return h1Rn;
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
        R = 0.0;
        R_err = 0.0;

        if uls > 0.5:
            if lspp * lsmm > 0.5:
                R = uls / ( 2 * TMath.Sqrt( lspp*lsmm ) );
                #R_err = TMath.Sqrt( ( pow(lspp*lsmm_err*uls,2) + pow(lspp_err*lsmm*uls,2) + 4*pow(lspp*lsmm*uls_err,2) ) / ( 16 * pow(lspp*lsmm,3) ) );
                R_err = abs(R) * math.sqrt(math.pow(lspp_err/lspp/2, 2) + math.pow(lsmm_err/lsmm/2, 2) + math.pow(uls_err/uls, 2));
            elif lspp + lsmm > 0.5 :
                R = uls / ( 2 * 0.5 * (lspp + lsmm) );
                R_err = TMath.Sqrt( (pow(uls*lspp_err,2) + pow(uls*lsmm_err,2) + pow(lspp+lsmm,2)*pow(uls_err,2) ) / pow(lspp+lsmm,4) );
        else:
            R = 0.0;
            R_err = 0.0;
        if uls < 1.5 and lspp < 1.5 and lsmm < 1.5:
            R = 0.0;
            R_err = 0.0;
        #print("mix im+1 = {0} , x = {1} , R = {2} , uls = {3} , lspp = {4} , lsmm = {5}".format(im+1,h1R.GetBinCenter(im+1),R,uls,lspp,lsmm));
        h1R.SetBinContent(im+1, R);
        h1R.SetBinError(im+1, R_err);
    return h1R;
#______________________________________________________________________
def get_R_factor_1number(n_uls_mix, n_lspp_mix, n_lsmm_mix, n_uls_mix_err, n_lspp_mix_err, n_lsmm_mix_err):
    R = 0.0;
    R_err = 0.0;
    if n_uls_mix > 0.5:
        if n_lspp_mix * n_lsmm_mix > 0.5: #geometrical mean
            R = n_uls_mix / ( 2 * math.sqrt( n_lspp_mix*n_lsmm_mix ) );
            R_err = abs(R) * math.sqrt(math.pow(n_lspp_mix_err/n_lspp_mix/2, 2) + math.pow(n_lsmm_mix_err/n_lsmm_mix/2, 2) + math.pow(n_uls_mix_err/n_uls_mix, 2));
        elif n_lspp_mix + n_lsmm_mix > 0.5: #arithmetic mean
            R = n_uls_mix / ( 2 * 0.5 * (n_lspp_mix + n_lsmm_mix) );
            R_err = math.sqrt( (pow(n_uls_mix*n_lspp_mix_err,2) + pow(n_uls_mix*n_lsmm_mix_err,2) + pow(n_lspp_mix+n_lsmm_mix,2)*pow(n_uls_mix_err,2) ) / pow(n_lspp_mix+n_lsmm_mix,4) );
    else:
        R = 0.0;
        R_err = 0.0;
    if n_uls_mix < 1.5 and n_lspp_mix < 1.5 and n_lsmm_mix < 1.5: #only for protection
        R = 0.0;
        R_err = 0.0;
    return [R, R_err];
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
        if R > 1e-6:
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
        h1bkg.SetBinContent(im+1, bkg);
        h1bkg.SetBinError(im+1, bkg_err);
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
        if lspp * lsmm > 1e-6 : #geometrical mean
            bkg = 2.0 * R * TMath.Sqrt(lspp * lsmm);
            bkg_err = TMath.Sqrt( pow(2 * math.sqrt(lspp*lsmm) * R_err, 2) + (R*R)/(lspp*lsmm) * (pow(lspp*lsmm_err,2) + pow(lsmm*lspp_err,2) ) );
        elif lspp + lsmm > 1e-6 : #arithmetic mean
            bkg = 2.0 * R * 0.5 * (lspp + lsmm);
            bkg_err = TMath.Sqrt(  pow(lspp + lsmm, 2) * pow(R_err, 2) + R*R * (lspp_err*lspp_err + lsmm_err*lsmm_err) );
        h1bkg.SetBinContent(im+1,bkg);
        h1bkg.SetBinError(im+1,bkg_err);

    return h1bkg;
#______________________________________________________________________
def rebin_flow_thn(hs_yield, hs_vn, xmin, xmax, ymin, ymax, zmin, zmax):
    #the number of bins has to be consistent between hs_yield and hs_vn.
    list_n = [];
    list_n_err = [];
    list_vn = [];
    list_vn_err = [];

    for ix in range (0, hs_yield.GetAxis(0).GetNbins()):
        x_center = hs_yield.GetAxis(0).GetBinCenter(ix+1);
        if x_center < xmin or xmax < x_center:
            continue;

        for iy in range (0, hs_yield.GetAxis(1).GetNbins()):
            y_center = hs_yield.GetAxis(1).GetBinCenter(iy+1);
            if y_center < ymin or ymax < y_center:
                continue;

            for iz in range (0, hs_yield.GetAxis(2).GetNbins()):
                z_center = hs_yield.GetAxis(2).GetBinCenter(iz+1);
                if z_center < zmin or zmax < z_center:
                    continue;
                n     = hs_yield.GetBinContent(np.array([ix+1, iy+1, iz+1], dtype=np.int32));
                n_err = hs_yield.GetBinError(np.array([ix+1, iy+1, iz+1], dtype=np.int32));
                vn     = hs_vn.GetBinContent(np.array([ix+1, iy+1, iz+1], dtype=np.int32));
                vn_err = hs_vn.GetBinError(np.array([ix+1, iy+1, iz+1], dtype=np.int32));
                list_n.append(n);
                list_n_err.append(n_err);
                list_vn.append(vn);
                list_vn_err.append(vn_err);

    #print("list_n = ", list_n);
    #print("list_n_err = ", list_n_err);
    #print("list_vn = ", list_vn);
    #print("list_vn_err = ", list_vn_err);

    arr_n_vn = np.multiply(np.array(list_n, dtype=float), np.array(list_vn, dtype=float));
    arr_n2 = np.multiply(np.array(list_n, dtype=float), np.array(list_n, dtype=float));
    arr_vn_err2 = np.multiply(np.array(list_vn_err, dtype=float), np.array(list_vn_err, dtype=float));
    arr_n2_vn_err2 = np.multiply(arr_n2, arr_vn_err2);

    n_sum = np.sum(np.array(list_n, dtype=float));
    #print("arr_n_vn = ", arr_n_vn);
    vn = np.sum(arr_n_vn) / n_sum;
    #print("vn, vn_err = ", vn, vn_err);
    vn_err = math.sqrt(np.sum(arr_n2_vn_err2)) / n_sum;
    return [vn, vn_err];
#______________________________________________________________________
#______________________________________________________________________
#______________________________________________________________________

