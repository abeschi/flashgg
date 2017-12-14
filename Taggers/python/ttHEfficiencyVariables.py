dipho_variables=[
        "dipho_sumpt                     := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.sumPt : -100",
        "dipho_mass                      := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.mass : -100",
        "dipho_leadPt                    := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.pt : -100",
        "dipho_leadEt                    := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.et : -100",
        "dipho_leadEta                   := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.eta : -100",
        "dipho_leadPhi                   := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.phi : -100",
        "dipho_leadEnergy                := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.energy : -100",
        "dipho_lead_sieie                := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.sigmaIetaIeta : -100",
        "dipho_lead_hoe                  := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.hadronicOverEm : -100",
        "dipho_lead_sigmaEoE             := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.sigEOverE : -100",
        "dipho_lead_ptoM                 := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.pt/GetTTH[0].diPhoton.mass : -100",
        "dipho_leadR9                    := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.full5x5_r9 : -100",
        "dipho_leadEgChargedHadronIso    := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.egChargedHadronIso : -100",
        "dipho_subleadPt                 := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.pt : -100",
        "dipho_subleadEt                 := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.et : -100",
        "dipho_subleadEta                := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.eta : -100",
        "dipho_subleadPhi                := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.phi : -100",
        "dipho_subleadEnergy             := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.energy : -100",
        "dipho_sublead_sieie             := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.sigmaIetaIeta : -100",
        "dipho_sublead_hoe               := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.hadronicOverEm : -100",
        "dipho_sublead_sigmaEoE          := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.sigEOverE : -100",
        "dipho_sublead_ptoM              := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.pt/GetTTH[0].diPhoton.mass : -100",
        "dipho_subleadR9                 := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.full5x5_r9 : -100",
        "dipho_subleadEgChargedHadronIso := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.egChargedHadronIso : -100",
        "dipho_leadIDMVA                 := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingView.phoIdMvaWrtChosenVtx : -100",
        "dipho_subleadIDMVA              := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingView.phoIdMvaWrtChosenVtx : -100",
        "dipho_lead_elveto               := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.leadingPhoton.passElectronVeto : -100",
        "dipho_sublead_elveto            := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.subLeadingPhoton.passElectronVeto : -100",
        "dipho_mva                       := ? GetTTH.size()>0 ? GetTTH[0].diPhotonMVA.result : -100",
        "passHLT                         := ? GetTTH.size()>0 ? GetTTH[0].passHLT : -100",
        "passPreselection                := ? GetTTH.size()>0 ? GetTTH[0].passPreselection : -100",
        "dipho_vertexZ                   := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.vtx.z : -100",
        "dipho_vertexSigmaZ              := ? GetTTH.size()>0 ? GetTTH[0].diPhoton.vtx.zError : -100",
    ]

hadronic_variables=[
        "jet_pt1            :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>0 ? GetTTH[0].jetVector[0].pt : -100 ",
	"jet_eta1           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>0 ? GetTTH[0].jetVector[0].eta : -100 ",
	"jet_phi1           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>0 ? GetTTH[0].jetVector[0].phi : -100 ",
	"jet_energy1        :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>0 ? GetTTH[0].jetVector[0].energy : -100 ",
	"jet_bdiscriminant_CSV1 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>0 ? GetTTH[0].jetVector[0].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_bdiscriminant_MVACSV1 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>0 ? GetTTH[0].jetVector[0].bDiscriminator('pfCombinedMVAV2BJetTags') : -100",
	"jet_bdiscriminant_DeepCSV1 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>0 ? GetTTH[0].jetVector[0].bDiscriminator('pfDeepCSVJetTags') : -100",
	"jet_hadronFlavour1 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>0 ? GetTTH[0].jetVector[0].hadronFlavour : -100",
	"jet_pt2            :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>1 ? GetTTH[0].jetVector[1].pt : -100 ",
	"jet_eta2           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>1 ? GetTTH[0].jetVector[1].eta : -100 ",
	"jet_phi2           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>1 ? GetTTH[0].jetVector[1].phi : -100 ",
	"jet_energy2        :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>1 ? GetTTH[0].jetVector[1].energy : -100 ",
	"jet_bdiscriminant2 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>1 ? GetTTH[0].jetVector[1].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_hadronFlavour2 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>1 ? GetTTH[0].jetVector[1].hadronFlavour : -100",
	"jet_pt3            :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>2 ? GetTTH[0].jetVector[2].pt : -100 ",
	"jet_eta3           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>2 ? GetTTH[0].jetVector[2].eta : -100 ",
	"jet_phi3           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>2 ? GetTTH[0].jetVector[2].phi : -100 ",
	"jet_energy3        :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>2 ? GetTTH[0].jetVector[2].energy : -100 ",
	"jet_bdiscriminant3 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>2 ? GetTTH[0].jetVector[2].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour3 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>2 ? GetTTH[0].jetVector[2].hadronFlavour : -100",
	"jet_pt4            :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>3 ? GetTTH[0].jetVector[3].pt : -100 ",
	"jet_eta4           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>3 ? GetTTH[0].jetVector[3].eta : -100 ",
	"jet_phi4           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>3 ? GetTTH[0].jetVector[3].phi : -100 ",
	"jet_energy4        :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>3 ? GetTTH[0].jetVector[3].energy : -100 ",
	"jet_bdiscriminant4 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>3 ? GetTTH[0].jetVector[3].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour4 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>3 ? GetTTH[0].jetVector[3].hadronFlavour : -100",
	"jet_pt5            :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>4 ? GetTTH[0].jetVector[4].pt : -100 ",
	"jet_eta5           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>4 ? GetTTH[0].jetVector[4].eta : -100 ",
	"jet_phi5           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>4 ? GetTTH[0].jetVector[4].phi : -100 ",
	"jet_energy5        :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>4 ? GetTTH[0].jetVector[4].energy : -100 ",
	"jet_bdiscriminant5 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>4 ? GetTTH[0].jetVector[4].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour5 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>4 ? GetTTH[0].jetVector[4].hadronFlavour : -100",
	"jet_pt6            :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>5 ? GetTTH[0].jetVector[5].pt : -100 ",
	"jet_eta6           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>5 ? GetTTH[0].jetVector[5].eta : -100 ",
	"jet_phi6           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>5 ? GetTTH[0].jetVector[5].phi : -100",
	"jet_energy6        :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>5 ? GetTTH[0].jetVector[5].energy : -100 ",
	"jet_bdiscriminant6 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>5 ? GetTTH[0].jetVector[5].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour6 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>5 ? GetTTH[0].jetVector[5].hadronFlavour : -100",
	"jet_pt7            :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>6 ? GetTTH[0].jetVector[6].pt : -100 ",
	"jet_eta7           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>6 ? GetTTH[0].jetVector[6].eta : -100 ",
	"jet_phi7           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>6 ? GetTTH[0].jetVector[6].phi : -100",
	"jet_energy7        :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>6 ? GetTTH[0].jetVector[6].energy : -100 ",
	"jet_bdiscriminant7 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>6 ? GetTTH[0].jetVector[6].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour7 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>6 ? GetTTH[0].jetVector[6].hadronFlavour : -100",
	"jet_pt8            :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>7 ? GetTTH[0].jetVector[7].pt : -100 ",
	"jet_eta8           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>7 ? GetTTH[0].jetVector[7].eta : -100 ",
	"jet_phi8           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>7 ? GetTTH[0].jetVector[7].phi : -100",
	"jet_energy8        :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>7 ? GetTTH[0].jetVector[7].energy : -100 ",
	"jet_bdiscriminant8 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>7 ? GetTTH[0].jetVector[7].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour8 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>7 ? GetTTH[0].jetVector[7].hadronFlavour : -100",
	"jet_pt9            :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>8 ? GetTTH[0].jetVector[8].pt : -100 ",
	"jet_eta9           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>8 ? GetTTH[0].jetVector[8].eta : -100 ",
	"jet_phi9           :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>8 ? GetTTH[0].jetVector[8].phi : -100",
	"jet_energy9        :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>8 ? GetTTH[0].jetVector[8].energy : -100 ",
	"jet_bdiscriminant9 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>8 ? GetTTH[0].jetVector[8].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour9 :=  ? GetTTH.size()>0 && GetTTH[0].jetVector.size()>8 ? GetTTH[0].jetVector[8].hadronFlavour : -100",
	"njet               :=  ? GetTTH.size()>0 ? GetTTH[0].nJet : -100",
	"nbjetLoose         :=  ? GetTTH.size()>0 ? GetTTH[0].nBLoose : -100",
	"nbjetMedium        :=  ? GetTTH.size()>0 ? GetTTH[0].nBMedium : -100",
	"nbjetTight         :=  ? GetTTH.size()>0 ? GetTTH[0].nBTight : -100",
	"ttHMVA             :=  ? GetTTH.size()>0 ? GetTTH[0].tthMvaRes : -100",
	"MetPt              :=  ? GetTTH.size()>0 ? GetTTH[0].MetPt : -100",
	"MetPhi             :=  ? GetTTH.size()>0 ? GetTTH[0].MetPhi : -100"
]

truth_photon_variables=[
	"GenLeadPhton_pt         :=  LeadingPhotonPt ",
	"GenLeadPhton_eta        :=  LeadingPhotonEta ",
	"GenLeadPhton_phi        :=  LeadingPhotonPhi ",
	"GenLeadPhton_energy     :=  LeadingPhotonEnergy ",
	"GenSubleadPhton_pt      :=  SubleadingPhotonPt ",
	"GenSubleadPhton_eta     :=  SubleadingPhotonEta ",
	"GenSubleadPhton_phi     :=  SubleadingPhotonPhi ",
	"GenSubleadPhton_energy  :=  SubleadingPhotonEnergy ",
	"HiggsVtxZ               :=  higgsVertex.z"
    ]



leptonic_variables=[
	"mu_pt1             :=  ? GetTTH.size()>0 && GetTTH[0].muons.size()>0 ? GetTTH[0].muons[0].pt() : -100 ",
	"mu_eta1            :=  ? GetTTH.size()>0 && GetTTH[0].muons.size()>0 ? GetTTH[0].muons[0].eta() : -100 ",
	"mu_phi1            :=  ? GetTTH.size()>0 && GetTTH[0].muons.size()>0 ? GetTTH[0].muons[0].phi() : -100 ",
	"mu_pt2             :=  ? GetTTH.size()>0 && GetTTH[0].muons.size()>1 ? GetTTH[0].muons[1].pt() : -100 ",
	"mu_eta2            :=  ? GetTTH.size()>0 && GetTTH[0].muons.size()>1 ? GetTTH[0].muons[1].eta() : -100 ",
	"mu_phi2            :=  ? GetTTH.size()>0 && GetTTH[0].muons.size()>1 ? GetTTH[0].muons[1].phi() : -100 ",
	"ele_pt1            :=  ? GetTTH.size()>0 && GetTTH[0].electrons.size()>0 ? GetTTH[0].electrons[0].pt() : -100 ",
	"ele_eta1           :=  ? GetTTH.size()>0 && GetTTH[0].electrons.size()>0 ? GetTTH[0].electrons[0].eta() : -100 ",
	"ele_phi1           :=  ? GetTTH.size()>0 && GetTTH[0].electrons.size()>0 ? GetTTH[0].electrons[0].phi() : -100 ",
	"ele_pt2            :=  ? GetTTH.size()>0 && GetTTH[0].electrons.size()>1 ? GetTTH[0].electrons[1].pt() : -100 ",
	"ele_eta2           :=  ? GetTTH.size()>0 && GetTTH[0].electrons.size()>1 ? GetTTH[0].electrons[1].eta() : -100 ",
	"ele_phi2           :=  ? GetTTH.size()>0 && GetTTH[0].electrons.size()>1 ? GetTTH[0].electrons[1].phi() : -100 ",
	"jet_pt1            :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>0 ? GetTTH[0].jets[0].pt() : -100 ",
	"jet_eta1           :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>0 ? GetTTH[0].jets[0].eta() : -100 ",
	"jet_phi1           :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>0 ? GetTTH[0].jets[0].phi() : -100 ",
	"jet_energy1        :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>0 ? GetTTH[0].jets[0].energy() : -100 ",
	"jet_bdiscriminant1 :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>0 ? GetTTH[0].jets[0].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_pt2            :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>1 ? GetTTH[0].jets[1].pt() : -100 ",
	"jet_eta2           :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>1 ? GetTTH[0].jets[1].eta() : -100 ",
	"jet_phi2           :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>1 ? GetTTH[0].jets[1].phi() : -100 ",
	"jet_energy2        :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>1 ? GetTTH[0].jets[1].energy() : -100 ",
	"jet_bdiscriminant2 :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>1 ? GetTTH[0].jets[1].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_pt3            :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>2 ? GetTTH[0].jets[2].pt() : -100 ",
	"jet_eta3           :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>2 ? GetTTH[0].jets[2].eta() : -100 ",
	"jet_phi3           :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>2 ? GetTTH[0].jets[2].phi() : -100 ",
	"jet_energy3        :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>2 ? GetTTH[0].jets[2].energy() : -100 ",
	"jet_bdiscriminant3 :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>2 ? GetTTH[0].jets[2].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt4            :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>3 ? GetTTH[0].jets[3].pt() : -100 ",
	"jet_eta4           :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>3 ? GetTTH[0].jets[3].eta() : -100 ",
	"jet_phi4           :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>3 ? GetTTH[0].jets[3].phi() : -100 ",
	"jet_energy4        :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>3 ? GetTTH[0].jets[3].energy() : -100 ",
	"jet_bdiscriminant4 :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>3 ? GetTTH[0].jets[3].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt5            :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>4 ? GetTTH[0].jets[4].pt() : -100 ",
	"jet_eta5           :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>4 ? GetTTH[0].jets[4].eta() : -100 ",
	"jet_phi5           :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>4 ? GetTTH[0].jets[4].phi() : -100 ",
	"jet_energy5        :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>4 ? GetTTH[0].jets[4].energy() : -100 ",
	"jet_bdiscriminant5 :=  ? GetTTH.size()>0 && GetTTH[0].jets.size()>4 ? GetTTH[0].jets[4].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"njet               :=  ? GetTTH.size()>0 ? GetTTH[0].nJet : -100",
	"nbjetLoose         :=  ? GetTTH.size()>0 ? GetTTH[0].nBLoose : -100",
	"nbjetMedium        :=  ? GetTTH.size()>0 ? GetTTH[0].nBMedium : -100",
	"nbjetTight         :=  ? GetTTH.size()>0 ? GetTTH[0].nBTight : -100",
	"MetPt              :=  ? GetTTH.size()>0 ? GetTTH[0].MetPt : -100",
	"MetPhi             :=  ? GetTTH.size()>0 ? GetTTH[0].MetPhi : -100"
]



efficiency_hadronic_variables=[
        "GenJet_pt1            :=  ? GenJets.size() > 0 ? GenJets[0].pt : -100 ",
	"GenJet_eta1           :=  ? GenJets.size() > 0 ? GenJets[0].eta : -100 ",
	"GenJet_phi1           :=  ? GenJets.size() > 0 ? GenJets[0].phi : -100 ",
	"GenJet_energy1        :=  ? GenJets.size() > 0 ? GenJets[0].energy : -100 ",
        "GenJet_pt2            :=  ? GenJets.size() > 1 ? GenJets[1].pt : -100 ",
	"GenJet_eta2           :=  ? GenJets.size() > 1 ? GenJets[1].eta : -100 ",
	"GenJet_phi2           :=  ? GenJets.size() > 1 ? GenJets[1].phi : -100 ",
	"GenJet_energy2        :=  ? GenJets.size() > 1 ? GenJets[1].energy : -100 ",
        "GenJet_pt3            :=  ? GenJets.size() > 2 ? GenJets[2].pt : -100 ",
	"GenJet_eta3           :=  ? GenJets.size() > 2 ? GenJets[2].eta : -100 ",
	"GenJet_phi3           :=  ? GenJets.size() > 2 ? GenJets[2].phi : -100 ",
	"GenJet_energy3        :=  ? GenJets.size() > 2 ? GenJets[2].energy : -100 ",
        "GenJet_pt4            :=  ? GenJets.size() > 3 ? GenJets[3].pt : -100 ",
	"GenJet_eta4           :=  ? GenJets.size() > 3 ? GenJets[3].eta : -100 ",
	"GenJet_phi4           :=  ? GenJets.size() > 3 ? GenJets[3].phi : -100 ",
	"GenJet_energy4        :=  ? GenJets.size() > 3 ? GenJets[3].energy : -100 ",
	"NGenJets              :=  ? GenJets.size() > 0 ? NGenJets : 0",
	"NGenJetsEta2p4        :=  ? GenJets.size() > 0 ? NGenJetsEta2p4 : 0 ",
	"Vtx0Z                 :=  vertex0.z",
	"Vtx0ZSigma            :=  ? vertex0.z !=-137 ?  vertex0.zError : -100 ",
	"diphotonSize          :=  NDiphotons",
        "hasRecoLepton         := ? GetTTH.size()>0 ? GetTTH[0].HasRecoLeptons : -100",
   ]



efficiency_leptonic_variables=[
        "GenJet_pt1            :=  ? GenJets.size() > 0 ? GenJets[0].pt : -100 ",
	"GenJet_eta1           :=  ? GenJets.size() > 0 ? GenJets[0].eta : -100 ",
	"GenJet_phi1           :=  ? GenJets.size() > 0 ? GenJets[0].phi : -100 ",
	"GenJet_energy1        :=  ? GenJets.size() > 0 ? GenJets[0].energy : -100 ",
        "GenJet_pt2            :=  ? GenJets.size() > 1 ? GenJets[1].pt : -100 ",
	"GenJet_eta2           :=  ? GenJets.size() > 1 ? GenJets[1].eta : -100 ",
	"GenJet_phi2           :=  ? GenJets.size() > 1 ? GenJets[1].phi : -100 ",
	"GenJet_energy2        :=  ? GenJets.size() > 1 ? GenJets[1].energy : -100 ",
        "GenJet_pt3            :=  ? GenJets.size() > 2 ? GenJets[2].pt : -100 ",
	"GenJet_eta3           :=  ? GenJets.size() > 2 ? GenJets[2].eta : -100 ",
	"GenJet_phi3           :=  ? GenJets.size() > 2 ? GenJets[2].phi : -100 ",
	"GenJet_energy3        :=  ? GenJets.size() > 2 ? GenJets[2].energy : -100 ",
        "GenJet_pt4            :=  ? GenJets.size() > 3 ? GenJets[3].pt : -100 ",
	"GenJet_eta4           :=  ? GenJets.size() > 3 ? GenJets[3].eta : -100 ",
	"GenJet_phi4           :=  ? GenJets.size() > 3 ? GenJets[3].phi : -100 ",
	"GenJet_energy4        :=  ? GenJets.size() > 3 ? GenJets[3].energy : -100 ",
	"NGenJets              :=  ? GenJets.size() > 0 ? NGenJets : 0",
	"NGenJetsEta2p4        :=  ? GenJets.size() > 0 ? NGenJetsEta2p4 : 0 ",
	"Vtx0Z                 :=  vertex0.z",
	"Vtx0ZSigma            :=  ? vertex0.z !=-137 ?  vertex0.zError : -100 ",
	"diphotonSize          :=  NDiphotons",
	"GenLepton_pt1         :=  ? GenLeptons.size() > 0 ? GenLeptons[0].pt : -100 ",
	"GenLepton_eta1        :=  ? GenLeptons.size() > 0 ? GenLeptons[0].eta : -100 ",
	"GenLepton_phi1        :=  ? GenLeptons.size() > 0 ? GenLeptons[0].phi : -100 ",
	"GenLepton_energy1     :=  ? GenLeptons.size() > 0 ? GenLeptons[0].energy : -100 ",
	"GenLepton_pdgId1      :=  ? GenLeptons.size() > 0 ? GenLeptons[0].pdgId : -100 ",
	"GenLepton_pt2         :=  ? GenLeptons.size() > 1 ? GenLeptons[1].pt : -100 ",
	"GenLepton_eta2        :=  ? GenLeptons.size() > 1 ? GenLeptons[1].eta : -100 ",
	"GenLepton_phi2        :=  ? GenLeptons.size() > 1 ? GenLeptons[1].phi : -100 ",
	"GenLepton_energy2     :=  ? GenLeptons.size() > 1 ? GenLeptons[1].energy : -100 ",
	"GenLepton_pdgId2      :=  ? GenLeptons.size() > 1 ? GenLeptons[1].pdgId : -100 "
   ]


















































