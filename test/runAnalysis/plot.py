import ROOT




h_masseff=ROOT.TH1D("masseff","masseff",200,0,20000)
h_CSeff=ROOT.TH1D("masseff","masseff",100,-1,1)
f_deno=ROOT.TFile("deno.root")
f_no=ROOT.TFile("no.root")
mass_deno=f_deno.Get("Our2018MuonsPlusMuonsMinusHistos/DimuonMassVertexConstrained")
mass_no=f_no.Get("Our2018MuonsPlusMuonsMinusHistos/DimuonMassVertexConstrained")
CS_deno=f_deno.Get("Our2018MuonsPlusMuonsMinusHistos/CosThetaStarDilepton")
CS_no=f_no.Get("Our2018MuonsPlusMuonsMinusHistos/CosThetaStarDilepton")
mass_noBin=mass_no.Rebin(100)
mass_denoBin=mass_deno.Rebin(100)
h_masseff.Add(mass_noBin)
h_masseff.Divide(mass_denoBin)
h_CSeff.Add(CS_no)
h_CSeff.Divide(CS_deno)
c1=ROOT.TCanvas("c1","c1",800,800)
c1.SetLogx()
h_masseff.GetXaxis().SetTitle("M [GeV]")
h_masseff.GetXaxis().SetRangeUser(4500,6000)
h_masseff.GetYaxis().SetTitle("eff")
h_masseff.GetXaxis().SetMoreLogLabels()
h_masseff.Draw("hist")
c1.Print("masseff.pdf")
c2=ROOT.TCanvas("c2","c2",800,800)
h_CSeff.GetXaxis().SetTitle("Cos #theta")
h_CSeff.GetYaxis().SetTitle("eff")
h_CSeff.Draw("hist")
c2.Print("CSeff.pdf")

