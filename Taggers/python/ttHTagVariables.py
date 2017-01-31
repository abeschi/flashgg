dipho_variables=[
    "dipho_sumpt            := diPhoton.sumPt",
    "dipho_mass             := diPhoton.mass",
    "dipho_leadPt           := diPhoton.leadingPhoton.pt",
    "dipho_leadEt           := diPhoton.leadingPhoton.et",
    "dipho_leadEta          := diPhoton.leadingPhoton.eta",
    "dipho_leadPhi          := diPhoton.leadingPhoton.phi",
    "dipho_lead_sieie       := diPhoton.leadingPhoton.sigmaIetaIeta",
    "dipho_lead_hoe         := diPhoton.leadingPhoton.hadronicOverEm",
    "dipho_lead_sigmaEoE    := diPhoton.leadingPhoton.sigEOverE",
    "dipho_lead_ptoM        := diPhoton.leadingPhoton.pt/diPhoton.mass",
    "dipho_leadR9           := diPhoton.leadingPhoton.r9",
    "dipho_subleadPt        := diPhoton.subLeadingPhoton.pt",
    "dipho_subleadEt        := diPhoton.subLeadingPhoton.et",
    "dipho_subleadEta       := diPhoton.subLeadingPhoton.eta",
    "dipho_subleadPhi       := diPhoton.subLeadingPhoton.phi",
    "dipho_sublead_sieie    := diPhoton.subLeadingPhoton.sigmaIetaIeta",
    "dipho_sublead_hoe      := diPhoton.subLeadingPhoton.hadronicOverEm",
    "dipho_sublead_sigmaEoE := diPhoton.subLeadingPhoton.sigEOverE",
    "dipho_sublead_ptoM     := diPhoton.subLeadingPhoton.pt/diPhoton.mass",
    "dipho_subleadR9        := diPhoton.subLeadingPhoton.r9",
    "dipho_leadIDMVA        := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
    "dipho_subleadIDMVA     := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
    "dipho_lead_elveto      := diPhoton.leadingPhoton.passElectronVeto",
    "dipho_sublead_elveto   := diPhoton.subLeadingPhoton.passElectronVeto",
    "dipho_mva              := diPhotonMVA.result"
    ]


hadronic_variables=[
        "jet_pt1            :=  ? jetVector.size()>0 ? jetVector[0].pt : -100 ",
	"jet_phi1           :=  ? jetVector.size()>0 ? jetVector[0].eta : -100 ",
	"jet_eta1           :=  ? jetVector.size()>0 ? jetVector[0].phi : -100 ",
	"jet_bdiscriminant1 :=  ? jetVector.size()>0 ? jetVector[0].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_pt2            :=  ? jetVector.size()>1 ? jetVector[1].pt : -100 ",
	"jet_phi2           :=  ? jetVector.size()>1 ? jetVector[1].eta : -100 ",
	"jet_eta2           :=  ? jetVector.size()>1 ? jetVector[1].phi : -100 ",
	"jet_bdiscriminant2 :=  ? jetVector.size()>1 ? jetVector[1].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_pt3            :=  ? jetVector.size()>2 ? jetVector[2].pt : -100 ",
	"jet_phi3           :=  ? jetVector.size()>2 ? jetVector[2].eta : -100 ",
	"jet_eta3           :=  ? jetVector.size()>2 ? jetVector[2].phi : -100 ",
	"jet_bdiscriminant3 :=  ? jetVector.size()>2 ? jetVector[2].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt4            :=  ? jetVector.size()>3 ? jetVector[3].pt : -100 ",
	"jet_phi4           :=  ? jetVector.size()>3 ? jetVector[3].eta : -100 ",
	"jet_eta4           :=  ? jetVector.size()>3 ? jetVector[3].phi : -100 ",
	"jet_bdiscriminant4 :=  ? jetVector.size()>3 ? jetVector[3].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt5            :=  ? jetVector.size()>4 ? jetVector[4].pt : -100 ",
	"jet_phi5           :=  ? jetVector.size()>4 ? jetVector[4].eta : -100 ",
	"jet_eta5           :=  ? jetVector.size()>4 ? jetVector[4].phi : -100 ",
	"jet_bdiscriminant5 :=  ? jetVector.size()>4 ? jetVector[4].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt6            :=  ? jetVector.size()>5 ? jetVector[5].pt : -100 ",
	"jet_phi6           :=  ? jetVector.size()>5 ? jetVector[5].eta : -100 ",
	"jet_eta6           :=  ? jetVector.size()>5 ? jetVector[5].phi : -100",
	"jet_bdiscriminant6 :=  ? jetVector.size()>5 ? jetVector[5].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"njet               :=  nJet ",
	"nbjet              :=  nBMedium "

]


leptonic_variables=[
	"mu_pt1         :=  ? muons.size()>0 ? muons[0].pt() : -100 ",
	"mu_phi1        :=  ? muons.size()>0 ? muons[0].eta() : -100 ",
	"mu_eta1        :=  ? muons.size()>0 ? muons[0].phi() : -100 ",
	"mu_pt2         :=  ? muons.size()>1 ? muons[1].pt() : -100 ",
	"mu_phi2        :=  ? muons.size()>1 ? muons[1].eta() : -100 ",
	"mu_eta2        :=  ? muons.size()>1 ? muons[1].phi() : -100 ",
	"ele_pt1        :=  ? electrons.size()>0 ? electrons[0].pt() : -100 ",
	"ele_phi1       :=  ? electrons.size()>0 ? electrons[0].eta() : -100 ",
	"ele_eta1       :=  ? electrons.size()>0 ? electrons[0].phi() : -100 ",
	"ele_pt2        :=  ? electrons.size()>1 ? electrons[1].pt() : -100 ",
	"ele_phi2       :=  ? electrons.size()>1 ? electrons[1].eta() : -100 ",
	"ele_eta2       :=  ? electrons.size()>1 ? electrons[1].phi() : -100 ",
	"jet_pt1        :=  ? jets.size()>0 ? jets[0].pt() : -100 ",
	"jet_phi1       :=  ? jets.size()>0 ? jets[0].eta() : -100 ",
	"jet_eta1       :=  ? jets.size()>0 ? jets[0].phi() : -100 ",
	"jet_bdiscriminant1 :=  ? jets.size()>0 ? jets[0].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_pt2        :=  ? jets.size()>1 ? jets[1].pt() : -100 ",
	"jet_phi2       :=  ? jets.size()>1 ? jets[1].eta() : -100 ",
	"jet_eta2       :=  ? jets.size()>1 ? jets[1].phi() : -100 ",
	"jet_bdiscriminant2 :=  ? jets.size()>1 ? jets[1].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100",
	"jet_pt3        :=  ? jets.size()>2 ? jets[2].pt() : -100 ",
	"jet_phi3       :=  ? jets.size()>2 ? jets[2].eta() : -100 ",
	"jet_eta3       :=  ? jets.size()>2 ? jets[2].phi() : -100 ",
	"jet_bdiscriminant3 :=  ? jets.size()>2 ? jets[2].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 ",
	"jet_pt4        :=  ? jets.size()>3 ? jets[3].pt() : -100 ",
	"jet_phi4       :=  ? jets.size()>3 ? jets[3].eta() : -100 ",
	"jet_eta4       :=  ? jets.size()>3 ? jets[3].phi() : -100 ",
	"jet_bdiscriminant4 :=  ? jets.size()>3 ? jets[3].bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags') : -100 "
]
truth_variables=[
    
    ]