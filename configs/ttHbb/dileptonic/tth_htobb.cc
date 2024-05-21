#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TObject.h"
#include "TMath.h"
#include "TString.h"

#include "config.h"

int main(int argc, char* argv[]) {

    if (argc<3) {
      std::cout << "Usage of the program : tth_htobb [config_file] [file_list]" << std::endl;
      std::cout << "config_file : file with all necessary configurations to be used by the program" << std::endl;
      std::cout << "file_list : file with list of ROOT files to be processed, name will be used for output files" << std::endl;
      exit(1);
    }

    config cfg(argv[1]);

    std::ifstream file_list(argv[2]);

    std::string file_out_name(argv[2]);

    gErrorIgnoreLevel = kError;

    /////////////////////// *** Configuration *** /////////////////////////////
    const int number_of_events_to_process                             = cfg.get<int>("number_of_events_to_process");

    ////////////////////////// *** Events *** /////////////////////////////////
    TChain Events("Events", "Events");
    std::string iFile;
    while (file_list >> iFile) Events.Add(TString(iFile));

    // create the reader
    TTreeReader Reader(&Events);
    Reader.Restart();
    if ( number_of_events_to_process>0 ) Reader.SetEntriesRange(0,number_of_events_to_process);

    std::cout << std::endl;
    std::cout << "   ==> " << Events.GetNtrees() << " files to be processed in '" << file_out_name << "' file list ... " << std::endl;
    std::cout << "   ==> " << int((number_of_events_to_process>0) ? TMath::Min(number_of_events_to_process,int(Events.GetEntries())) : Events.GetEntries()) << " events to be processed out of " << Events.GetEntries() << " in '" << Events.GetName() << "' chain ... " << std::endl;

    ////////////////////////// *** Variable Declaration *** /////////////////////////////////
    TTreeReaderValue<UInt_t> run = {Reader, "run"};
    TTreeReaderValue<UInt_t> luminosityBlock = {Reader, "luminosityBlock"};
    TTreeReaderValue<ULong64_t> event = {Reader, "event"};

    TTreeReaderValue<Float_t> genWeight = {Reader, "genWeight"};

    TTreeReaderValue<UInt_t> nGenPart = {Reader, "nGenPart"};
    TTreeReaderArray<Float_t> GenPart_pt = {Reader, "GenPart_pt"};
    TTreeReaderArray<Float_t> GenPart_phi = {Reader, "GenPart_phi"};
    TTreeReaderArray<Float_t> GenPart_eta = {Reader, "GenPart_eta"};
    TTreeReaderArray<Float_t> GenPart_mass = {Reader, "GenPart_mass"};
    TTreeReaderArray<Int_t> GenPart_pdgId = {Reader, "GenPart_pdgId"};
    TTreeReaderArray<Int_t> GenPart_statusFlags = {Reader, "GenPart_statusFlags"};
    TTreeReaderArray<Int_t> GenPart_genPartIdxMother = {Reader, "GenPart_genPartIdxMother"};
    
    /*
    TTreeReaderValue<UInt_t> nGenPart = {Reader, "nLHEPart"};
    TTreeReaderArray<Float_t> GenPart_pt = {Reader, "LHEPart_pt"};
    TTreeReaderArray<Float_t> GenPart_phi = {Reader, "LHEPart_phi"};
    TTreeReaderArray<Float_t> GenPart_eta = {Reader, "LHEPart_eta"};
    TTreeReaderArray<Float_t> GenPart_mass = {Reader, "LHEPart_mass"};
    TTreeReaderArray<Int_t> GenPart_pdgId = {Reader, "LHEPart_pdgId"};
    TTreeReaderArray<Int_t> GenPart_statusFlags = {Reader, "LHEPart_status"};
    */
    
    

    std::vector<Int_t> particle_type = {5,11,12,13,14,-5,-11,-12,-13,-14};
    const std::map<Int_t, Double_t> particle_mass = {{5,4.18},{11,0.0005110},{12,0.0},{13,0.105658},{14,0.0},};
    std::vector<TString> particles = {"top_b","top_b_bar","h_b","h_b_bar","l","l_bar","nu","nu_bar"};
    const Int_t n_particles = particles.size();
    
    // create output file
    TFile* file_out = TFile::Open(TString(file_out_name)+".root","RECREATE");
    file_out->cd();

    // add the histograms
    TH1::SetDefaultSumw2();
    
    // Top quarks
    TH1D* Histo_Top_Mass = new TH1D("Histo_Top_Mass", "Histo_Top_Mass", 100, 0., 300.);
    TH1D* Histo_AntiTop_Mass = new TH1D("Histo_AntiTop_Mass", "Histo_AntiTop_Mass", 100, 0., 300.);
    TH1D* Histo_Top_Pt = new TH1D("Histo_Top_Pt", "Histo_Top_Pt", 20, 0., 500.);
    TH1D* Histo_AntiTop_Pt = new TH1D("Histo_AntiTop_Pt", "Histo_AntiTop_Pt", 20, 0., 500.);
    TH1D* Histo_Top_Eta = new TH1D("Histo_Top_Eta", "Histo_Top_Eta", 20, -5., 5.);
    TH1D* Histo_AntiTop_Eta = new TH1D("Histo_AntiTop_Eta", "Histo_AntiTop_Eta", 20, -5., 5.);
    
    // W bosons
    TH1D* Histo_Wp_Mass = new TH1D("Histo_Wp_Mass", "Histo_Wp_Mass", 100, 0., 200.);
    TH1D* Histo_Wm_Mass = new TH1D("Histo_Wm_Mass", "Histo_Wm_Mass", 100, 0., 200.);
    TH1D* Histo_Wp_Pt = new TH1D("Histo_Wp_Pt", "Histo_Wp_Pt", 20, 0., 500.);
    TH1D* Histo_Wm_Pt = new TH1D("Histo_Wm_Pt", "Histo_Wm_Pt", 20, 0., 500.);
    TH1D* Histo_Wp_Eta = new TH1D("Histo_Wp_Eta", "Histo_Wp_Eta", 20, -5., 5.);
    TH1D* Histo_Wm_Eta = new TH1D("Histo_Wm_Eta", "Histo_Wm_Eta", 20, -5., 5.);
    
    // B quarks
    TH1D* Histo_B_Top_Pt = new TH1D("Histo_B_Top_Pt", "Histo_B_Top_Pt", 20, 0., 500.);
    TH1D* Histo_AntiB_Top_Pt = new TH1D("Histo_AntiB_Top_Pt", "Histo_AntiB_Top_Pt", 20, 0., 500.);
    TH1D* Histo_B_Top_Eta = new TH1D("Histo_B_Top_Eta", "Histo_B_Top_Eta", 20, -5., 5.);
    TH1D* Histo_AntiB_Top_Eta = new TH1D("Histo_AntiB_Top_Eta", "Histo_AntiB_Top_Eta", 20, -5., 5.);
    TH1D* Histo_B_H_Pt = new TH1D("Histo_B_H_Pt", "Histo_B_H_Pt", 20, 0., 500.);
    TH1D* Histo_AntiB_H_Pt = new TH1D("Histo_AntiB_H_Pt", "Histo_AntiB_H_Pt", 20, 0., 500.);
    TH1D* Histo_B_H_Eta = new TH1D("Histo_B_H_Eta", "Histo_B_H_Eta", 20, -5., 5.);
    TH1D* Histo_AntiB_H_Eta = new TH1D("Histo_AntiB_H_Eta", "Histo_AntiB_H_Eta", 20, -5., 5.);
    
    // Higgs/Gluon
    TH1D* Histo_H_Mass = new TH1D("Histo_H_Mass", "Histo_H_Mass", 100, 0., 200.);
    TH1D* Histo_H_Pt = new TH1D("Histo_H_Pt", "Histo_H_Pt", 20, 0., 500.);
    TH1D* Histo_H_Eta = new TH1D("Histo_H_Eta", "Histo_H_Eta", 20, -5., 5.);
    
    // Reconstruction
    TH1D* Histo_BAntiB_H_DeltaR = new TH1D("Histo_BAntiB_H_DeltaR", "Histo_BAntiB_H_DeltaR", 100, 0., 10.);
    TH1D* Histo_CHEL = new TH1D("Histo_CHEL", "Histo_CHEL", 100, -1., 1.);
    TH1D* Histo_B_H_CosPhi = new TH1D("Histo_B_H_CosPhi", "Histo_B_H_CosPhi", 100, -1., 1.);
    
    

    //// *** Running the Event Loop *** ////
    while (Reader.Next()) { // Loop Over Events

        if ( (Reader.GetCurrentEntry()+1)%100000==0 ) std::cout << "       --> processing event number " << Reader.GetCurrentEntry()+1 <<  std::endl;

        double weight = (*genWeight.Get());

        std::map<TString, Int_t> particle_map = {};
        for (int iPart=0; iPart<*nGenPart.Get(); iPart++) {    // Loop over all particles in given event 

            if ( particle_map.size() >= n_particles ) break;

            if ( std::find(particle_type.begin(), particle_type.end(), GenPart_pdgId[iPart]) != particle_type.end() ) {

                if ( (GenPart_statusFlags[iPart] & 256) && (GenPart_statusFlags[iPart] & 4096) ) {
                    
                    // b quark from top
                    if (particle_map.count("top_b") == 0) {
                        if ( GenPart_pdgId[iPart] == 5) {
                            Int_t mother_idx = GenPart_genPartIdxMother[iPart];
                            while (GenPart_genPartIdxMother[mother_idx]>0 && abs(GenPart_pdgId[mother_idx]) != 6) mother_idx = GenPart_genPartIdxMother[mother_idx];
                            if ( abs(GenPart_pdgId[mother_idx]) == 6 ) { particle_map["top_b"] = iPart; continue; };
                        }
                    }
                    // anti-b quark from anti-top
                    if (particle_map.count("top_b_bar") == 0) {
                        if ( GenPart_pdgId[iPart] == -5) {
                            Int_t mother_idx = GenPart_genPartIdxMother[iPart];
                            while (GenPart_genPartIdxMother[mother_idx]>0 && abs(GenPart_pdgId[mother_idx]) != 6) mother_idx = GenPart_genPartIdxMother[mother_idx];
                            if ( abs(GenPart_pdgId[mother_idx]) == 6 ) { particle_map["top_b_bar"] = iPart; continue; };
                        }
                    }
                    // b quark from boson
                    if (particle_map.count("h_b") == 0) {
                        if ( GenPart_pdgId[iPart] == 5) {
                            Int_t mother_idx = GenPart_genPartIdxMother[iPart];
                            while (GenPart_genPartIdxMother[mother_idx]>0) mother_idx = GenPart_genPartIdxMother[mother_idx];
                            if ( abs(GenPart_pdgId[mother_idx]) != 6 ) { particle_map["h_b"] = iPart; continue; };
                        }
                    }
                    // anti-b quark from boson
                    if (particle_map.count("h_b_bar") == 0) {
                        if ( GenPart_pdgId[iPart] == -5) {
                            Int_t mother_idx = GenPart_genPartIdxMother[iPart];
                            while (GenPart_genPartIdxMother[mother_idx]>0) mother_idx = GenPart_genPartIdxMother[mother_idx];
                            if ( abs(GenPart_pdgId[mother_idx]) != 6 ) { particle_map["h_b_bar"] = iPart; continue; };
                        }
                    }
                    // lepton
                    if (particle_map.count("l") == 0) {
                        if ( GenPart_pdgId[iPart] == 11 || GenPart_pdgId[iPart] == 13 ) {
                            Int_t mother_idx = GenPart_genPartIdxMother[iPart];
                            while (GenPart_genPartIdxMother[mother_idx]>0 && abs(GenPart_pdgId[mother_idx]) != 6) mother_idx = GenPart_genPartIdxMother[mother_idx];
                            if ( abs(GenPart_pdgId[mother_idx]) == 6 ) { particle_map["l"] = iPart; continue; };
                        }
                    }
                    // anti lepton
                    if (particle_map.count("l_bar") == 0) {
                        if ( GenPart_pdgId[iPart] == -11 || GenPart_pdgId[iPart] == -13 ) {
                            Int_t mother_idx = GenPart_genPartIdxMother[iPart];
                            while (GenPart_genPartIdxMother[mother_idx]>0 && abs(GenPart_pdgId[mother_idx]) != 6) mother_idx = GenPart_genPartIdxMother[mother_idx];
                            if ( abs(GenPart_pdgId[mother_idx]) == 6 ) { particle_map["l_bar"] = iPart; continue; };
                        }
                    }
                    // neutrino
                    if (particle_map.count("nu") == 0) {
                        if ( GenPart_pdgId[iPart] == 12 || GenPart_pdgId[iPart] == 14 ) {
                            Int_t mother_idx = GenPart_genPartIdxMother[iPart];
                            while (GenPart_genPartIdxMother[mother_idx]>0 && abs(GenPart_pdgId[mother_idx]) != 6) mother_idx = GenPart_genPartIdxMother[mother_idx];
                            if ( abs(GenPart_pdgId[mother_idx]) == 6 ) { particle_map["nu"] = iPart; continue; };
                        }
                    }
                    // anti neutrino
                    if (particle_map.count("nu_bar") == 0) {
                        if ( GenPart_pdgId[iPart] == -12 || GenPart_pdgId[iPart] == -14 ) {
                            Int_t mother_idx = GenPart_genPartIdxMother[iPart];
                            while (GenPart_genPartIdxMother[mother_idx]>0 && abs(GenPart_pdgId[mother_idx]) != 6) mother_idx = GenPart_genPartIdxMother[mother_idx];
                            if ( abs(GenPart_pdgId[mother_idx]) == 6 ) { particle_map["nu_bar"] = iPart; continue; };
                        }
                    }
                    
                }

            }

        }

        if ( particle_map.size() != n_particles ) continue;

        TLorentzVector B_Top;
        B_Top.SetPtEtaPhiM(GenPart_pt[particle_map["top_b"]],GenPart_eta[particle_map["top_b"]],GenPart_phi[particle_map["top_b"]],particle_mass.at(abs(GenPart_pdgId[particle_map["top_b"]])));
        TLorentzVector AntiB_Top;
        AntiB_Top.SetPtEtaPhiM(GenPart_pt[particle_map["top_b_bar"]],GenPart_eta[particle_map["top_b_bar"]],GenPart_phi[particle_map["top_b_bar"]],particle_mass.at(abs(GenPart_pdgId[particle_map["top_b_bar"]])));
        TLorentzVector B_H;
        B_H.SetPtEtaPhiM(GenPart_pt[particle_map["h_b"]],GenPart_eta[particle_map["h_b"]],GenPart_phi[particle_map["h_b"]],particle_mass.at(abs(GenPart_pdgId[particle_map["h_b"]])));
        TLorentzVector AntiB_H;
        AntiB_H.SetPtEtaPhiM(GenPart_pt[particle_map["h_b_bar"]],GenPart_eta[particle_map["h_b_bar"]],GenPart_phi[particle_map["h_b_bar"]],particle_mass.at(abs(GenPart_pdgId[particle_map["h_b_bar"]])));
        TLorentzVector Lepton;
        Lepton.SetPtEtaPhiM(GenPart_pt[particle_map["l"]],GenPart_eta[particle_map["l"]],GenPart_phi[particle_map["l"]],particle_mass.at(abs(GenPart_pdgId[particle_map["l"]])));
        TLorentzVector AntiLepton;
        AntiLepton.SetPtEtaPhiM(GenPart_pt[particle_map["l_bar"]],GenPart_eta[particle_map["l_bar"]],GenPart_phi[particle_map["l_bar"]],particle_mass.at(abs(GenPart_pdgId[particle_map["l_bar"]])));
        TLorentzVector Neutrino;
        Neutrino.SetPtEtaPhiM(GenPart_pt[particle_map["nu"]],GenPart_eta[particle_map["nu"]],GenPart_phi[particle_map["nu"]],particle_mass.at(abs(GenPart_pdgId[particle_map["nu"]])));
        TLorentzVector AntiNeutrino;
        AntiNeutrino.SetPtEtaPhiM(GenPart_pt[particle_map["nu_bar"]],GenPart_eta[particle_map["nu_bar"]],GenPart_phi[particle_map["nu_bar"]],particle_mass.at(abs(GenPart_pdgId[particle_map["nu_bar"]])));
        
        // Reconstructing intermidiate particles
        TLorentzVector Wp = AntiLepton + Neutrino;
        TLorentzVector Wm = Lepton + AntiNeutrino;
        TLorentzVector Top = B_Top + Wp;
        TLorentzVector AntiTop = AntiB_Top + Wm;
        TLorentzVector H = B_H + AntiB_H;

        // Spin polarization
        auto ttbar_boost_vector = -1 * (Top+AntiTop).BoostVector();

        auto top_ZMF = Top;
        auto atop_ZMF = AntiTop;
        auto lepton_atop_ZMF = Lepton;
        auto alepton_top_ZMF = AntiLepton;

        top_ZMF.Boost(ttbar_boost_vector);
        atop_ZMF.Boost(ttbar_boost_vector);
        lepton_atop_ZMF.Boost(ttbar_boost_vector);
        alepton_top_ZMF.Boost(ttbar_boost_vector);

        lepton_atop_ZMF.Boost(-1 * atop_ZMF.BoostVector());
        alepton_top_ZMF.Boost(-1 * top_ZMF.BoostVector());
        
        Double_t chel = lepton_atop_ZMF.Vect().Unit().Dot(alepton_top_ZMF.Vect().Unit());
        
        auto h_boostv3 = -H.BoostVector();

        auto b_h_ZMF = B_H;

        b_h_ZMF.Boost(h_boostv3);

        double cosphi_B_H = H.Vect().Unit().Dot(b_h_ZMF.Vect().Unit());

        // filling histograms
        
        Histo_Top_Mass->Fill(Top.M(), weight);
        Histo_AntiTop_Mass->Fill(AntiTop.M(), weight);
        Histo_Top_Pt->Fill(Top.Pt(), weight);
        Histo_AntiTop_Pt->Fill(AntiTop.Pt(), weight);
        Histo_Top_Eta->Fill(Top.Eta(), weight);
        Histo_AntiTop_Eta->Fill(AntiTop.Eta(), weight);
        
        Histo_Wp_Mass->Fill(Wp.M(), weight);
        Histo_Wm_Mass->Fill(Wm.M(), weight);
        Histo_Wp_Pt->Fill(Wp.Pt(), weight);
        Histo_Wm_Pt->Fill(Wm.Pt(), weight);
        Histo_Wp_Eta->Fill(Wp.Eta(), weight);
        Histo_Wm_Eta->Fill(Wm.Eta(), weight);
        
        Histo_B_Top_Pt->Fill(B_Top.Pt(), weight);
        Histo_AntiB_Top_Pt->Fill(AntiB_Top.Pt(), weight);
        Histo_B_Top_Eta->Fill(B_Top.Eta(), weight);
        Histo_AntiB_Top_Eta->Fill(AntiB_Top.Eta(), weight);
        Histo_B_H_Pt->Fill(B_H.Pt(), weight);
        Histo_AntiB_H_Pt->Fill(AntiB_H.Pt(), weight);
        Histo_B_H_Eta->Fill(B_H.Eta(), weight);
        Histo_AntiB_H_Eta->Fill(AntiB_H.Eta(), weight);
        
        Histo_H_Mass->Fill(H.M(), weight);
        Histo_H_Pt->Fill(H.Pt(), weight);
        Histo_H_Eta->Fill(H.Eta(), weight);
        
        Histo_BAntiB_H_DeltaR->Fill(B_H.DeltaR(AntiB_H), weight);
        Histo_CHEL->Fill(chel, weight);
        Histo_B_H_CosPhi->Fill(cosphi_B_H, weight);
       
    }

    file_out->Write();
    file_out->Close();

    std::cout <<"   ==> done." << std::endl;
    std::cout << std::endl;

    return 0;
}

