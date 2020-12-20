from coffea.analysis_objects import JaggedCandidateArray

pionmass=0.13957018

import awkward as ak
import numpy as np
import pyjet

def makeJets(vectors,R):

    # cluster jets for a single event from track vectors 
    sequence = pyjet.cluster(vectors, R=R, p=-1)
    psuedojets = sequence.inclusive_jets()
    #print(psuedojets)

    jets = []
    for jet in psuedojets: 
        # save a dict
        jet_dict = { "pt":jet.pt, 
                "eta":jet.eta, 
                "phi":jet.phi,
                "mass":jet.mass,
                "ntracks":len(jet) }
        jets.append(jet_dict)

    return jets 

def convertToJagged(jets):
    # takes a list of list of dicts
    # converts to JaggedArrayCandidates

    # make awkward array
    jagged_jets = ak.JaggedArray.fromiter(jets)

    # make jagged candidate array
    fatjets = JaggedCandidateArray.candidatesfromcounts(jagged_jets.counts,
                                                  pt=jagged_jets.pt.flatten(),
                                                  eta=jagged_jets.eta.flatten(),
                                                  phi=jagged_jets.phi.flatten(),
                                                  mass=jagged_jets.mass.flatten(),
                                                  ntracks=jagged_jets.ntracks.flatten())
    return fatjets 

def makeFatJets(tracks):

    print("making fat jets")
    # clusters fat jets from tracks
    # makes jets with three different radii
    # "anti-pattern"

    # list of list of dicts
    jets_20 = []
    jets_15 = []
    jets_10 = []

    for ievt in range(tracks.size):
        
        # for debug
        #if ievt > 10: continue
        if ievt % 1000 == 0: print("Event {}".format(ievt))

        # make a structured array from tracks

        vectors = np.zeros(tracks.pt[ievt].size, np.dtype([('pT', 'f8'), ('eta', 'f8'),
                                                ('phi', 'f8'), ('mass', 'f8')]) )
        vectors['pT'  ] = tracks.pt[ievt]
        vectors['eta' ] = tracks.eta[ievt]
        vectors['phi' ] = tracks.phi[ievt]
        vectors['mass'] = tracks.mass[ievt]
        #print(vectors)

        # make events jets 
        jets_20.append(makeJets(vectors, 2.0))
        jets_15.append(makeJets(vectors, 1.5))
        jets_10.append(makeJets(vectors, 1.0))

    fatjets_20 = convertToJagged(jets_20) 
    fatjets_15 = convertToJagged(jets_15) 
    fatjets_10 = convertToJagged(jets_10) 
    
    return (fatjets_10, fatjets_15, fatjets_20)

class Objects:
    def __init__(self, df):
        #self.df = df
        print("Getting Electrons")
        self.electrons = JaggedCandidateArray.candidatesfromcounts(
            df['Electrons_charge'].counts,
            eta=df['Electrons.fCoordinates.fEta'].flatten(),
            phi=df['Electrons.fCoordinates.fPhi'].flatten(),
            pt=df['Electrons.fCoordinates.fPt'].flatten(),
            energy=df['Electrons.fCoordinates.fE'].flatten(),
            passIso=df['Electrons_passIso'].flatten(),
            charge=df['Electrons_charge'].flatten()
        )    
        print("Getting Tracks")
        print(df['Tracks_charge'].counts)
        self.tracks = JaggedCandidateArray.candidatesfromcounts(
            df['Tracks_charge'].counts,
            px=df['Tracks.fCoordinates.fX'].flatten(),
            py=df['Tracks.fCoordinates.fY'].flatten(),
            pz=df['Tracks.fCoordinates.fZ'].flatten(),
            mass=pionmass*abs(df['Tracks_charge']).flatten(),
            fromPV0=df['Tracks_fromPV0'].flatten(),
            matchedToPFCandidate=df['Tracks_matchedToPFCandidate'].flatten(),
            quality=df['Tracks_quality'].flatten()
        )
        print("Getting Muons")
        self.muons = JaggedCandidateArray.candidatesfromcounts(
            df['Muons_charge'].counts,
            pt=df['Muons.fCoordinates.fPt'].flatten(),
            eta=df['Muons.fCoordinates.fEta'].flatten(),
            phi=df['Muons.fCoordinates.fPhi'].flatten(),
            energy=df['Muons.fCoordinates.fE'].flatten(),
            passIso=df['Muons_passIso'].flatten(),
            charge=df['Muons_charge'].flatten()
        )
        print("Getting Jets")
        self.jets = JaggedCandidateArray.candidatesfromcounts(
            df['Jets_ptD'].counts,
            pt=df['Jets.fCoordinates.fPt'].flatten(),
            eta=df['Jets.fCoordinates.fEta'].flatten(),
            phi=df['Jets.fCoordinates.fPhi'].flatten(),
            energy=df['Jets.fCoordinates.fE'].flatten(),
            ptD=df['Jets_ptD'].flatten(),
            bDeepCSVprob=df['Jets_bJetTagDeepCSVBvsAll'].flatten(),
        )
        
        (self.fatjets10, self.fatjets15, self.fatjets20) = makeFatJets(self.tracks)
        self.trackEtaCut = 2.5
        self.etaCut = 2.4

    def goodElectrons(self):
        # # Good Electrons
        electronQualityCut = (self.electrons.pt > 37) & (abs(self.electrons.eta) < self.etaCut)
        return self.electrons[electronQualityCut]

    def goodTracks(self):
        # # Good Tracks 
        trackQualityCut = (self.tracks.pt > 1) & (abs(self.tracks.eta) < self.trackEtaCut) & (self.tracks.fromPV0 >= 2 ) & ( self.tracks.matchedToPFCandidate > 0 )
        return self.tracks[trackQualityCut]

    def goodMuons(self):
        # # Good Muons
        muonQualityCut = (self.muons.pt > 30) & (abs(self.muons.eta) < self.etaCut)
        return self.muons[muonQualityCut]

    def goodJets(self):
        # # Good AK4 Jets Cut
        ak4QualityCut = (self.jets.pt > 30) & (abs(self.jets.eta) < self.etaCut)
        return self.jets[ak4QualityCut]
