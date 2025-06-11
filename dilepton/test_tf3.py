import ROOT
from ROOT import TF3, TH3D, TF2

f3 = TF3("f3", "1 + [0] * exp(-[1]*[1]*x*x/0.197/0.197 -[2]*[2]*y*y/0.197/0.197 -[3]*[3]*z*z/0.197/0.197)", -0.1, +0.1, -0.1, +0.1, -0.1, +0.1);
#f3 = TF2("f3", "1 + [0] * exp(-[1]*[1]*x*x/0.197/0.197 -[2]*[2]*y*y/0.197/0.197)", -0.1, +0.1, -0.1, +0.1);
f3.SetParameters(1, 5, 5, 5)
f3.Draw();
ROOT.SetOwnership(f3, False);

f3.SetNpx(200);
f3.SetNpy(200);
f3.SetNpz(200);

#qmin = -0.03;
#qmax = +0.03;
#
#h3_fit = f3.CreateHistogram();
#print(h3_fit.GetName());
#h3_fit.SetName("h3_fit");
#
#bin0 = h3_fit.GetXaxis().FindBin(qmin + 1e-6);
#bin1 = h3_fit.GetXaxis().FindBin(qmax - 1e-6);
#nbin = bin1 - bin0 + 1; #note : nbin is the same for x,y,z
#print("nbin_fit = ",nbin);
#
#print(h3_fit.GetNbinsX(), h3_fit.GetXaxis().GetXmin(), h3_fit.GetXaxis().GetXmax());
#print(h3_fit.GetNbinsY(), h3_fit.GetYaxis().GetXmin(), h3_fit.GetYaxis().GetXmax());
#print(h3_fit.GetNbinsZ(), h3_fit.GetZaxis().GetXmin(), h3_fit.GetZaxis().GetXmax());
#
#for ix in range(0, h3_fit.GetNbinsX()):
#    for iy in range(0, h3_fit.GetNbinsY()):
#        for iz in range(0, h3_fit.GetNbinsZ()):
#            x = h3_fit.GetXaxis().GetBinCenter(ix+1);
#            y = h3_fit.GetYaxis().GetBinCenter(iy+1);
#            z = h3_fit.GetZaxis().GetBinCenter(iz+1);
#            cf = f3.Eval(x,y,z);
#            h3_fit.SetBinContent(ix+1, iy+1, iz+1, cf);
#            h3_fit.SetBinError(ix+1, iy+1, iz+1, 0);
#h1out_fit  = h3_fit.ProjectionY("h1out_fit" , bin0, bin1, bin0, bin1, "");
#h1side_fit = h3_fit.ProjectionZ("h1side_fit", bin0, bin1, bin0, bin1, "");
#h1long_fit = h3_fit.ProjectionX("h1long_fit", bin0, bin1, bin0, bin1, "");
#
#h1out_fit.SetDirectory(0);
#h1out_fit.Scale(1/nbin/nbin);
#ROOT.SetOwnership(h1out_fit, False);
#
#h1side_fit.SetDirectory(0);
#h1side_fit.Scale(1/nbin/nbin);
#ROOT.SetOwnership(h1side_fit, False);
#
#h1long_fit.SetDirectory(0);
#h1long_fit.Scale(1/nbin/nbin);
#ROOT.SetOwnership(h1long_fit, False);
#
#h1out_fit.Draw("");
##h1side_fit.Draw("");
##h1long_fit.Draw("");
