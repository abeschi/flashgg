dipho_variables=[
    "dipho_sumpt            := diPhoton.sumPt",
    "dipho_mass             := diPhoton.mass",
    "dipho_leadPt           := diPhoton.leadingPhoton.pt",
    "dipho_leadEt           := diPhoton.leadingPhoton.et",
    "dipho_leadEta          := diPhoton.leadingPhoton.eta",
    "dipho_leadPhi          := diPhoton.leadingPhoton.phi",
    "dipho_leadEnergy          := diPhoton.leadingPhoton.energy",
    "dipho_lead_sieie       := diPhoton.leadingPhoton.sigmaIetaIeta",
    "dipho_lead_hoe         := diPhoton.leadingPhoton.hadronicOverEm",
    "dipho_lead_sigmaEoE    := diPhoton.leadingPhoton.sigEOverE",
    "dipho_lead_ptoM        := diPhoton.leadingPhoton.pt/diPhoton.mass",
    "dipho_leadR9           := diPhoton.leadingPhoton.full5x5_r9",
    "dipho_subleadPt        := diPhoton.subLeadingPhoton.pt",
    "dipho_subleadEt        := diPhoton.subLeadingPhoton.et",
    "dipho_subleadEta       := diPhoton.subLeadingPhoton.eta",
    "dipho_subleadPhi       := diPhoton.subLeadingPhoton.phi",
    "dipho_subleadEnergy    := diPhoton.subLeadingPhoton.energy",
    "dipho_sublead_sieie    := diPhoton.subLeadingPhoton.sigmaIetaIeta",
    "dipho_sublead_hoe      := diPhoton.subLeadingPhoton.hadronicOverEm",
    "dipho_sublead_sigmaEoE := diPhoton.subLeadingPhoton.sigEOverE",
    "dipho_sublead_ptoM     := diPhoton.subLeadingPhoton.pt/diPhoton.mass",
    "dipho_subleadR9        := diPhoton.subLeadingPhoton.full5x5_r9",
    "dipho_leadIDMVA        := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
    "dipho_subleadIDMVA     := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
    "dipho_lead_elveto      := diPhoton.leadingPhoton.passElectronVeto",
    "dipho_sublead_elveto   := diPhoton.subLeadingPhoton.passElectronVeto",
    "dipho_mva              := diPhotonMVA.result"
    ]


hadronic_variables=[
        "jet_pt1            :=  ? jetVector.size()>0 ? jetVector[0].pt : -100 ",
	"jet_eta1           :=  ? jetVector.size()>0 ? jetVector[0].eta : -100 ",
	"jet_phi1           :=  ? jetVector.size()>0 ? jetVector[0].phi : -100 ",
	"jet_energy1        :=  ? jetVector.size()>0 ? jetVector[0].energy : -100 ",
	"jet_bdiscriminant1 :=  ? jetVector.size()>0 ? jetVector[0].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_pt2            :=  ? jetVector.size()>1 ? jetVector[1].pt : -100 ",
	"jet_eta2           :=  ? jetVector.size()>1 ? jetVector[1].eta : -100 ",
	"jet_phi2           :=  ? jetVector.size()>1 ? jetVector[1].phi : -100 ",
	"jet_energy2        :=  ? jetVector.size()>1 ? jetVector[1].energy : -100 ",
	"jet_bdiscriminant2 :=  ? jetVector.size()>1 ? jetVector[1].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_pt3            :=  ? jetVector.size()>2 ? jetVector[2].pt : -100 ",
	"jet_eta3           :=  ? jetVector.size()>2 ? jetVector[2].eta : -100 ",
	"jet_phi3           :=  ? jetVector.size()>2 ? jetVector[2].phi : -100 ",
	"jet_energy3        :=  ? jetVector.size()>2 ? jetVector[2].energy : -100 ",
	"jet_bdiscriminant3 :=  ? jetVector.size()>2 ? jetVector[2].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt4            :=  ? jetVector.size()>3 ? jetVector[3].pt : -100 ",
	"jet_eta4           :=  ? jetVector.size()>3 ? jetVector[3].eta : -100 ",
	"jet_phi4           :=  ? jetVector.size()>3 ? jetVector[3].phi : -100 ",
	"jet_energy4        :=  ? jetVector.size()>3 ? jetVector[3].energy : -100 ",
	"jet_bdiscriminant4 :=  ? jetVector.size()>3 ? jetVector[3].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt5            :=  ? jetVector.size()>4 ? jetVector[4].pt : -100 ",
	"jet_eta5           :=  ? jetVector.size()>4 ? jetVector[4].eta : -100 ",
	"jet_phi5           :=  ? jetVector.size()>4 ? jetVector[4].phi : -100 ",
	"jet_energy5        :=  ? jetVector.size()>4 ? jetVector[4].energy : -100 ",
	"jet_bdiscriminant5 :=  ? jetVector.size()>4 ? jetVector[4].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt6            :=  ? jetVector.size()>5 ? jetVector[5].pt : -100 ",
	"jet_eta6           :=  ? jetVector.size()>5 ? jetVector[5].eta : -100 ",
	"jet_phi6           :=  ? jetVector.size()>5 ? jetVector[5].phi : -100",
	"jet_energy6        :=  ? jetVector.size()>5 ? jetVector[5].energy : -100 ",
	"jet_bdiscriminant6 :=  ? jetVector.size()>5 ? jetVector[5].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt7            :=  ? jetVector.size()>6 ? jetVector[6].pt : -100 ",
	"jet_eta7           :=  ? jetVector.size()>6 ? jetVector[6].eta : -100 ",
	"jet_phi7           :=  ? jetVector.size()>6 ? jetVector[6].phi : -100",
	"jet_energy7        :=  ? jetVector.size()>6 ? jetVector[6].energy : -100 ",
	"jet_bdiscriminant7 :=  ? jetVector.size()>6 ? jetVector[6].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt8            :=  ? jetVector.size()>7 ? jetVector[7].pt : -100 ",
	"jet_eta8           :=  ? jetVector.size()>7 ? jetVector[7].eta : -100 ",
	"jet_phi8           :=  ? jetVector.size()>7 ? jetVector[7].phi : -100",
	"jet_energy8        :=  ? jetVector.size()>7 ? jetVector[7].energy : -100 ",
	"jet_bdiscriminant8 :=  ? jetVector.size()>7 ? jetVector[7].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt9            :=  ? jetVector.size()>8 ? jetVector[8].pt : -100 ",
	"jet_eta9           :=  ? jetVector.size()>8 ? jetVector[8].eta : -100 ",
	"jet_phi9           :=  ? jetVector.size()>8 ? jetVector[8].phi : -100",
	"jet_energy9        :=  ? jetVector.size()>8 ? jetVector[8].energy : -100 ",
	"jet_bdiscriminant9 :=  ? jetVector.size()>8 ? jetVector[8].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"njet               :=  nJet ",
	"nbjet              :=  nBMedium ",
	"ttHMVA             :=  tthMvaRes ",
	"MetPt              :=  MetPt ",
	"MetPhi             :=  MetPhi "
]


leptonic_variables=[
	"mu_pt1             :=  ? muons.size()>0 ? muons[0].pt() : -100 ",
	"mu_eta1            :=  ? muons.size()>0 ? muons[0].eta() : -100 ",
	"mu_phi1            :=  ? muons.size()>0 ? muons[0].phi() : -100 ",
	"mu_pt2             :=  ? muons.size()>1 ? muons[1].pt() : -100 ",
	"mu_eta2            :=  ? muons.size()>1 ? muons[1].eta() : -100 ",
	"mu_phi2            :=  ? muons.size()>1 ? muons[1].phi() : -100 ",
	"ele_pt1            :=  ? electrons.size()>0 ? electrons[0].pt() : -100 ",
	"ele_eta1           :=  ? electrons.size()>0 ? electrons[0].eta() : -100 ",
	"ele_phi1           :=  ? electrons.size()>0 ? electrons[0].phi() : -100 ",
	"ele_pt2            :=  ? electrons.size()>1 ? electrons[1].pt() : -100 ",
	"ele_eta2           :=  ? electrons.size()>1 ? electrons[1].eta() : -100 ",
	"ele_phi2           :=  ? electrons.size()>1 ? electrons[1].phi() : -100 ",
	"jet_pt1            :=  ? jets.size()>0 ? jets[0].pt() : -100 ",
	"jet_eta1           :=  ? jets.size()>0 ? jets[0].eta() : -100 ",
	"jet_phi1           :=  ? jets.size()>0 ? jets[0].phi() : -100 ",
	"jet_energy1        :=  ? jets.size()>0 ? jets[0].energy() : -100 ",
	"jet_bdiscriminant1 :=  ? jets.size()>0 ? jets[0].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_pt2            :=  ? jets.size()>1 ? jets[1].pt() : -100 ",
	"jet_eta2           :=  ? jets.size()>1 ? jets[1].eta() : -100 ",
	"jet_phi2           :=  ? jets.size()>1 ? jets[1].phi() : -100 ",
	"jet_energy2        :=  ? jets.size()>1 ? jets[1].energy() : -100 ",
	"jet_bdiscriminant2 :=  ? jets.size()>1 ? jets[1].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_pt3            :=  ? jets.size()>2 ? jets[2].pt() : -100 ",
	"jet_eta3           :=  ? jets.size()>2 ? jets[2].eta() : -100 ",
	"jet_phi3           :=  ? jets.size()>2 ? jets[2].phi() : -100 ",
	"jet_energy3        :=  ? jets.size()>2 ? jets[2].energy() : -100 ",
	"jet_bdiscriminant3 :=  ? jets.size()>2 ? jets[2].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt4            :=  ? jets.size()>3 ? jets[3].pt() : -100 ",
	"jet_eta4           :=  ? jets.size()>3 ? jets[3].eta() : -100 ",
	"jet_phi4           :=  ? jets.size()>3 ? jets[3].phi() : -100 ",
	"jet_energy4        :=  ? jets.size()>3 ? jets[3].energy() : -100 ",
	"jet_bdiscriminant4 :=  ? jets.size()>3 ? jets[3].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt5            :=  ? jets.size()>4 ? jets[4].pt() : -100 ",
	"jet_eta5           :=  ? jets.size()>4 ? jets[4].eta() : -100 ",
	"jet_phi5           :=  ? jets.size()>4 ? jets[4].phi() : -100 ",
	"jet_energy5        :=  ? jets.size()>4 ? jets[4].energy() : -100 ",
	"jet_bdiscriminant5 :=  ? jets.size()>4 ? jets[4].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"MetPt              :=  MetPt ",
	"MetPhi             :=  MetPhi "
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



























