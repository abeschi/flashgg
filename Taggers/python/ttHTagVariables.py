






dipho_variables=[
    "dipho_sumpt                     := ? diphoIndex!=-1 ? diPhoton.sumPt : -100",
    "dipho_mass                      := ? diphoIndex!=-1 ? diPhoton.mass : -100",
    "dipho_leadPt                    := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.pt : -100",
    "dipho_leadEt                    := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.et : -100",
    "dipho_leadEta                   := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.eta : -100",
    "dipho_leadPhi                   := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.phi : -100",
    "dipho_leadEnergy                := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.energy : -100",
    "dipho_lead_sieie                := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.sigmaIetaIeta : -100",
    "dipho_lead_hoe                  := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.hadronicOverEm : -100",
    "dipho_lead_sigmaEoE             := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.sigEOverE : -100",
    "dipho_lead_ptoM                 := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.pt/diPhoton.mass : -100",
    "dipho_leadR9                    := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.full5x5_r9 : -100",
    "dipho_leadEgChargedHadronIso    := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.egChargedHadronIso : -100",
    "dipho_subleadPt                 := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.pt : -100",
    "dipho_subleadEt                 := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.et : -100",
    "dipho_subleadEta                := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.eta : -100",
    "dipho_subleadPhi                := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.phi : -100",
    "dipho_subleadEnergy             := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.energy : -100",
    "dipho_sublead_sieie             := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.sigmaIetaIeta : -100",
    "dipho_sublead_hoe               := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.hadronicOverEm : -100",
    "dipho_sublead_sigmaEoE          := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.sigEOverE : -100",
    "dipho_sublead_ptoM              := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.pt/diPhoton.mass : -100",
    "dipho_subleadR9                 := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.full5x5_r9 : -100",
    "dipho_subleadEgChargedHadronIso := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.egChargedHadronIso : -100",
    "dipho_leadIDMVA                 := ? diphoIndex!=-1 ? diPhoton.leadingView.phoIdMvaWrtChosenVtx : -100",
    "dipho_subleadIDMVA              := ? diphoIndex!=-1 ? diPhoton.subLeadingView.phoIdMvaWrtChosenVtx : -100",
    "dipho_lead_elveto               := ? diphoIndex!=-1 ? diPhoton.leadingPhoton.passElectronVeto : -100",
    "dipho_sublead_elveto            := ? diphoIndex!=-1 ? diPhoton.subLeadingPhoton.passElectronVeto : -100",
    "dipho_mva                       := ? diphoIndex!=-1 ? diPhotonMVA.result : -100",
    "dipho_index                     := diphoIndex"
    ]


hadronic_variables=[
        "jet_pt1            :=  ? jetVector.size()>0 ? jetVector[0].pt : -100 ",
	"jet_eta1           :=  ? jetVector.size()>0 ? jetVector[0].eta : -100 ",
	"jet_phi1           :=  ? jetVector.size()>0 ? jetVector[0].phi : -100 ",
	"jet_energy1        :=  ? jetVector.size()>0 ? jetVector[0].energy : -100 ",
	"jet_bdiscriminant1 :=  ? jetVector.size()>0 ? jetVector[0].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_hadronFlavour1 :=  ? jetVector.size()>0 ? jetVector[0].hadronFlavour : -100",
	"jet_pt2            :=  ? jetVector.size()>1 ? jetVector[1].pt : -100 ",
	"jet_eta2           :=  ? jetVector.size()>1 ? jetVector[1].eta : -100 ",
	"jet_phi2           :=  ? jetVector.size()>1 ? jetVector[1].phi : -100 ",
	"jet_energy2        :=  ? jetVector.size()>1 ? jetVector[1].energy : -100 ",
	"jet_bdiscriminant2 :=  ? jetVector.size()>1 ? jetVector[1].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_hadronFlavour2 :=  ? jetVector.size()>1 ? jetVector[1].hadronFlavour : -100",
	"jet_pt3            :=  ? jetVector.size()>2 ? jetVector[2].pt : -100 ",
	"jet_eta3           :=  ? jetVector.size()>2 ? jetVector[2].eta : -100 ",
	"jet_phi3           :=  ? jetVector.size()>2 ? jetVector[2].phi : -100 ",
	"jet_energy3        :=  ? jetVector.size()>2 ? jetVector[2].energy : -100 ",
	"jet_bdiscriminant3 :=  ? jetVector.size()>2 ? jetVector[2].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour3 :=  ? jetVector.size()>2 ? jetVector[2].hadronFlavour : -100",
	"jet_pt4            :=  ? jetVector.size()>3 ? jetVector[3].pt : -100 ",
	"jet_eta4           :=  ? jetVector.size()>3 ? jetVector[3].eta : -100 ",
	"jet_phi4           :=  ? jetVector.size()>3 ? jetVector[3].phi : -100 ",
	"jet_energy4        :=  ? jetVector.size()>3 ? jetVector[3].energy : -100 ",
	"jet_bdiscriminant4 :=  ? jetVector.size()>3 ? jetVector[3].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour4 :=  ? jetVector.size()>3 ? jetVector[3].hadronFlavour : -100",
	"jet_pt5            :=  ? jetVector.size()>4 ? jetVector[4].pt : -100 ",
	"jet_eta5           :=  ? jetVector.size()>4 ? jetVector[4].eta : -100 ",
	"jet_phi5           :=  ? jetVector.size()>4 ? jetVector[4].phi : -100 ",
	"jet_energy5        :=  ? jetVector.size()>4 ? jetVector[4].energy : -100 ",
	"jet_bdiscriminant5 :=  ? jetVector.size()>4 ? jetVector[4].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour5 :=  ? jetVector.size()>4 ? jetVector[4].hadronFlavour : -100",
	"jet_pt6            :=  ? jetVector.size()>5 ? jetVector[5].pt : -100 ",
	"jet_eta6           :=  ? jetVector.size()>5 ? jetVector[5].eta : -100 ",
	"jet_phi6           :=  ? jetVector.size()>5 ? jetVector[5].phi : -100",
	"jet_energy6        :=  ? jetVector.size()>5 ? jetVector[5].energy : -100 ",
	"jet_bdiscriminant6 :=  ? jetVector.size()>5 ? jetVector[5].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour6 :=  ? jetVector.size()>5 ? jetVector[5].hadronFlavour : -100",
	"jet_pt7            :=  ? jetVector.size()>6 ? jetVector[6].pt : -100 ",
	"jet_eta7           :=  ? jetVector.size()>6 ? jetVector[6].eta : -100 ",
	"jet_phi7           :=  ? jetVector.size()>6 ? jetVector[6].phi : -100",
	"jet_energy7        :=  ? jetVector.size()>6 ? jetVector[6].energy : -100 ",
	"jet_bdiscriminant7 :=  ? jetVector.size()>6 ? jetVector[6].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour7 :=  ? jetVector.size()>6 ? jetVector[6].hadronFlavour : -100",
	"jet_pt8            :=  ? jetVector.size()>7 ? jetVector[7].pt : -100 ",
	"jet_eta8           :=  ? jetVector.size()>7 ? jetVector[7].eta : -100 ",
	"jet_phi8           :=  ? jetVector.size()>7 ? jetVector[7].phi : -100",
	"jet_energy8        :=  ? jetVector.size()>7 ? jetVector[7].energy : -100 ",
	"jet_bdiscriminant8 :=  ? jetVector.size()>7 ? jetVector[7].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour8 :=  ? jetVector.size()>7 ? jetVector[7].hadronFlavour : -100",
	"jet_pt9            :=  ? jetVector.size()>8 ? jetVector[8].pt : -100 ",
	"jet_eta9           :=  ? jetVector.size()>8 ? jetVector[8].eta : -100 ",
	"jet_phi9           :=  ? jetVector.size()>8 ? jetVector[8].phi : -100",
	"jet_energy9        :=  ? jetVector.size()>8 ? jetVector[8].energy : -100 ",
	"jet_bdiscriminant9 :=  ? jetVector.size()>8 ? jetVector[8].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_hadronFlavour9 :=  ? jetVector.size()>8 ? jetVector[8].hadronFlavour : -100",
	"njet               :=  nJet",
	"nbjetLoose         :=  nBLoose",
	"nbjetMedium        :=  nBMedium",
	"nbjetTight         :=  nBTight",
	"ttHMVA             :=  tthMvaRes",
	"MetPt              :=  MetPt",
	"MetPhi             :=  MetPhi"
]


truth_photon_variables=[
	"GenLeadPhton_pt         :=  tagTruth().LeadingPhotonPt ",
	"GenLeadPhton_eta        :=  tagTruth().LeadingPhotonEta ",
	"GenLeadPhton_phi        :=  tagTruth().LeadingPhotonPhi ",
	"GenLeadPhton_energy     :=  tagTruth().LeadingPhotonEnergy ",
	"GenSubleadPhton_pt      :=  tagTruth().SubleadingPhotonPt ",
	"GenSubleadPhton_eta     :=  tagTruth().SubleadingPhotonEta ",
	"GenSubleadPhton_phi     :=  tagTruth().SubleadingPhotonPhi ",
	"GenSubleadPhton_energy  :=  tagTruth().SubleadingPhotonEnergy "
    ]


truth_hadronic_variables=[
	"Top1_pt1                := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[0].pt() : -100",
	"Top1_energy1            := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[0].energy() : -100",
	"Top1_eta1               := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[0].eta() : -100",
	"Top1_phi1               := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[0].phi() : -100",
	"Top1_pdgId1             := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[0].pdgId() : -100",
	"Top1_pt2                := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[1].pt() : -100",
	"Top1_energy2            := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[1].energy() : -100",
	"Top1_eta2               := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[1].eta() : -100",
	"Top1_phi2               := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[1].phi() : -100",
	"Top1_pdgId2             := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[1].pdgId() : -100",
	"Top1_pt3                := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[2].pt() : -100",
	"Top1_energy3            := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[2].energy() : -100",
	"Top1_eta3               := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[2].eta(): -100",
	"Top1_phi3               := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[2].phi() : -100",
	"Top1_pdgId3             := ? tagTruth().GenTop1.size()==3 ? tagTruth().GenTop1[2].pdgId() : -100",
	"Top2_pt1                := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[0].pt() : -100",
	"Top2_energy1            := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[0].energy(): -100",
	"Top2_eta1               := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[0].eta() : -100",
	"Top2_phi1               := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[0].phi() : -100",
	"Top2_pdgId1             := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[0].pdgId() : -100",
	"Top2_pt2                := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[1].pt() : 100",
	"Top2_energy2            := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[1].energy() : -100",
	"Top2_eta2               := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[1].eta() : -100",
	"Top2_phi2               := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[1].phi() : -100",
	"Top2_pdgId2             := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[1].pdgId() : -100",
	"Top2_pt3                := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[2].pt() : -100",
	"Top2_energy3            := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[2].energy() : -100",
	"Top2_eta3               := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[2].eta() : -100",
	"Top2_phi3               := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[2].phi() : -100",
	"Top2_pdgId3             := ? tagTruth().GenTop2.size()==3 ? tagTruth().GenTop2[2].pdgId() : -100",
    ]





efficiency_variables=[
        "GenJet_pt1            :=  ? tagTruth().GenJets.size() > 0 ? tagTruth().GenJets[0].pt : -100 ",
	"GenJet_eta1           :=  ? tagTruth().GenJets.size() > 0 ? tagTruth().GenJets[0].eta : -100 ",
	"GenJet_phi1           :=  ? tagTruth().GenJets.size() > 0 ? tagTruth().GenJets[0].phi : -100 ",
	"GenJet_energy1        :=  ? tagTruth().GenJets.size() > 0 ? tagTruth().GenJets[0].energy : -100 ",
        "GenJet_pt2            :=  ? tagTruth().GenJets.size() > 1 ? tagTruth().GenJets[1].pt : -100 ",
	"GenJet_eta2           :=  ? tagTruth().GenJets.size() > 1 ? tagTruth().GenJets[1].eta : -100 ",
	"GenJet_phi2           :=  ? tagTruth().GenJets.size() > 1 ? tagTruth().GenJets[1].phi : -100 ",
	"GenJet_energy2        :=  ? tagTruth().GenJets.size() > 1 ? tagTruth().GenJets[1].energy : -100 ",
        "GenJet_pt3            :=  ? tagTruth().GenJets.size() > 2 ? tagTruth().GenJets[2].pt : -100 ",
	"GenJet_eta3           :=  ? tagTruth().GenJets.size() > 2 ? tagTruth().GenJets[2].eta : -100 ",
	"GenJet_phi3           :=  ? tagTruth().GenJets.size() > 2 ? tagTruth().GenJets[2].phi : -100 ",
	"GenJet_energy3        :=  ? tagTruth().GenJets.size() > 2 ? tagTruth().GenJets[2].energy : -100 ",
        "GenJet_pt4            :=  ? tagTruth().GenJets.size() > 3 ? tagTruth().GenJets[3].pt : -100 ",
	"GenJet_eta4           :=  ? tagTruth().GenJets.size() > 3 ? tagTruth().GenJets[3].eta : -100 ",
	"GenJet_phi4           :=  ? tagTruth().GenJets.size() > 3 ? tagTruth().GenJets[3].phi : -100 ",
	"GenJet_energy4        :=  ? tagTruth().GenJets.size() > 3 ? tagTruth().GenJets[3].energy : -100 ",
	"NGenJets              :=  ? tagTruth().GenJets.size() > 0 ? tagTruth().NGenJets : 0",
	"NGenJetsEta2p4        :=  ? tagTruth().GenJets.size() > 0 ? tagTruth().NGenJetsEta2p4 : 0 ",
	"hasRecoLepton         :=  HasRecoLeptons",
	"passHLT               :=  passHLT",
	"passPreselection      :=  passPreselection"
   ]
























