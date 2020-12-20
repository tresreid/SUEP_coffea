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
    
    return (fatjets_20, fatjets_15, fatjets_10)
