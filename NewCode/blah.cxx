Double_t timePerTurn[6] = {3.454825, 10.463932, 32.299520, 85.955160, 193.332472, 379.186997};
Double_t lines[6] = {137, 321, 694, 1369, 2501, 4285};
Double_t order[6] = {2,3,4,5,6,7};

TGraph* gr1;
TGraph* gr2;
TGraph* gr3;

void doit()
{
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->Divide(3,1);
  c1->cd(1);
  gr1 = new TGraph(6,order,lines);
  gr1->SetMarkerStyle(20);
  gr1->Draw("ap");
  c1->cd(2);
  gr2 = new TGraph(6,order,timePerTurn);
  gr2->SetMarkerStyle(20);
  gr2->Draw("ap");
  c1->cd(3);
  gr3 = new TGraph(6,lines,timePerTurn);
  gr3->SetMarkerStyle(20);
  gr3->Draw("ap");

  TNtupleD* nt = new TNtupleD("nt","nt","o:l:t");
  nt->SetMarkerStyle(20);
  for (UInt_t i = 0; i<6; i++)
  {
    nt->Fill(order[i],lines[i],timePerTurn[i]);
  }
}
