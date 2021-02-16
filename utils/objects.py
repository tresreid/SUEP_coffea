from coffea.analysis_objects import JaggedCandidateArray

pionmass=0.13957018

import awkward as ak
import numpy as np
#import pyjet
from pyjet import cluster, DTYPE_PTEPM
import math

def jet_width(jet):
    girth = 0
    pt_avg =0
    count =0
    for i,constituent in enumerate([x for x in jet.constituents_array() if x['pT'] > 1]):
      if constituent['pT'] < 1:
        continue
      #min_dR = 9999
      #is_scalar = True
      #get scalar constituents
      phi = constituent['phi']
      eta = constituent['eta']
      pt = constituent['pT']
      dPhi = abs(jet.phi-phi)
      if dPhi > math.pi:
        dPhi = dPhi - 2*math.pi
      dEta = jet.eta-eta
      dR = math.sqrt(dEta*dEta + dPhi*dPhi)
      girth += pt * dR
      pt_avg += pt
      count += 1
      #if dR < min_dR:
      #  is_scalar = True
      #  min_dR = dR
    return ((girth / jet.pt),((pt_avg/count)/jet.pt))

def jet_constituents(jet,scalars,isrs):
    scalar_part = 0
    isr_part = 0
    for i,constituent in enumerate([x for x in jet.constituents_array() if x['pT'] > 1]):
      if constituent['pT'] < 1:
        continue
      min_dR = 9999
      is_scalar = True
      #get scalar constituents
      phi = constituent['phi']
      eta = constituent['eta']
      for scalar in scalars:
        #if abs(scalar.pt -constituent['pT'])/scalar.pt > 0.3:
        #  continue
        dPhi = abs(scalar['phi']-phi)
        if dPhi > math.pi:
          dPhi = dPhi - 2*math.pi
        dEta = scalar['eta']-eta
        dR = math.sqrt(dEta*dEta + dPhi*dPhi)
        if dR < min_dR:
          is_scalar = True
          min_dR = dR

      #get isr constituents
      for isr in isrs:
        #if abs(scalar.pt -constituent['pT'])/scalar.pt > 0.3:
        #  continue
        dPhi = abs(isr['phi']-phi)
        if dPhi > math.pi:
          dPhi = dPhi - 2*math.pi
        dEta = isr['eta']-eta
        dR = math.sqrt(dEta*dEta + dPhi*dPhi)
        if dR < min_dR:
          is_scalar = False
          min_dR = dR
      if is_scalar:
        scalar_part = scalar_part +1
      else:
        isr_part = isr_part +1
    return (scalar_part, isr_part)



def makeJets_Comparisons(jet_algo,R, particles):
  eta_min, eta_max = -4., 4.
  bins = 200
  cone = R*10
  eta = np.linspace(eta_min, eta_max, bins + 1)[:-1] + (eta_max - eta_min) / (2 * bins)
  phi = np.linspace(-np.pi, np.pi, bins + 1)[:-1] + (np.pi / bins)
  X, Y = np.meshgrid(phi, eta)
  event = np.zeros(len(particles),dtype=DTYPE_PTEPM)
  event['pT'] = particles['pT']#[p.pt for p in particles]
  event['eta'] =particles['eta']# [p.eta for p in particles]
  event['phi'] =particles['phi']# [p.phi for p in particles]
  ghosts = np.zeros(eta.shape[0] * phi.shape[0], dtype=DTYPE_PTEPM)
  ghosts['pT'] = 1e-8
  ghosts['eta'] = Y.ravel()
  ghosts['phi'] = X.ravel()

  # add ghosts to the event
  event = np.concatenate([event, ghosts], axis=0)

  # p = -1 (ak), 0 (CA), 1 (kt)
  sequence = cluster(event,R=R,p=jet_algo)
  jets = sequence.inclusive_jets(ptmin=30)
 # print(jets)
  return jets

def makeJets(vectors,R):

    # cluster jets for a single event from track vectors 
    #sequence = pyjet.cluster(vectors, R=R, p=-1)
    sequence = cluster(vectors, R=R, p=-1)
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

def makeFatJets(tracks,Ht):

    print("making fat jets")
    # clusters fat jets from tracks
    # makes jets with three different radii
    # "anti-pattern"

    # list of list of dicts
    jets_20 = []
    jets_15 = []
    jets_10 = []
    jets_08 = []
    events = []

    for ievt in range(tracks.size):
        
        # for debug
        if ievt < 6960: continue
        #print("Event {}".format(ievt))
        #if ievt % 1000 == 0: print("Event {}".format(ievt))
        if Ht[ievt] < 1000: continue
        print("Event {}/{} : Ht {}".format(ievt,tracks.size,Ht[ievt]))
        #print("Ht {}".format(Ht[ievt]))

        # make a structured array from tracks

        vectors = np.zeros(tracks.pt[ievt].size, np.dtype([('pT', 'f8'), ('eta', 'f8'),
                                                ('phi', 'f8'), ('mass', 'f8')]) )
        vectors['pT'  ] = tracks.pt[ievt]
        vectors['eta' ] = tracks.eta[ievt]
        vectors['phi' ] = tracks.phi[ievt]
        vectors['mass'] = tracks.mass[ievt]
        #print(vectors)

        # make events jets 
        events.append(ievt)
        jets_20.append([makeJets_Comparisons(-1,2.0,vectors),makeJets_Comparisons(0,2.0,vectors),makeJets_Comparisons(1,2.0,vectors)])
        jets_15.append([makeJets_Comparisons(-1,1.5,vectors),makeJets_Comparisons(0,1.5,vectors),makeJets_Comparisons(1,1.5,vectors)])
        jets_10.append([makeJets_Comparisons(-1,1.0,vectors),makeJets_Comparisons(0,1.0,vectors),makeJets_Comparisons(1,1.0,vectors)])
        jets_08.append([makeJets_Comparisons(-1,0.8,vectors),makeJets_Comparisons(0,0.8,vectors),makeJets_Comparisons(1,0.8,vectors)])
        #jets_20.append(makeJets(vectors, 2.0))
        #jets_15.append(makeJets(vectors, 1.5))
        #jets_10.append(makeJets(vectors, 1.0))

    #fatjets_20 = convertToJagged(jets_20) 
    #fatjets_15 = convertToJagged(jets_15) 
    #fatjets_10 = convertToJagged(jets_10) 
    #fatjets_20 = jets_20 
    #fatjets_15 = jets_15 
    #fatjets_10 = jets_10 
    #fatjets_08 = jets_08 
    
    return (events,(jets_08,jets_10, jets_15, jets_20))

def compareJets(tracks,scalarsx,isrsx,Ht,nvtx,numInteractions,fout):
  #new_jets = makeJets_Comparisons(-1,0.8,parts) 
  #(events,new_jets,j1,j2) = makeFatJets(tracks,Ht) 
  (events,all_jets) = makeFatJets(tracks,Ht) 
  #jet_algo_title = "akt"
  #Ri = 1.0

  for R_jets,Ri in zip(all_jets,[0.8,1.0,1.5,2.0]):
    for new_jets,jet_algo_title in zip(R_jets,["akt","ca","kt"]):
      for ievt,new_jet in zip(events,new_jets):
        scalars = np.zeros(scalarsx.pt[ievt].size, np.dtype([('pt', 'f8'), ('eta', 'f8'),
                                                    ('phi', 'f8'), ('mass', 'f8')]) )
        isrs = np.zeros(isrsx.pt[ievt].size, np.dtype([('pt', 'f8'), ('eta', 'f8'),
                                                    ('phi', 'f8'), ('mass', 'f8')]) )
        scalars['pt'] = scalarsx.pt[ievt]
        scalars['eta'] = scalarsx.eta[ievt]
        scalars['phi'] = scalarsx.phi[ievt]
        isrs['pt'] =  isrsx.pt[ievt]
        isrs['eta'] = isrsx.eta[ievt]
        isrs['phi'] = isrsx.phi[ievt]
        for jet_i in range(len(new_jet)):
          #print(new_jets[jet_i])
          (girth_all, trackpt_all) = jet_width(new_jet[jet_i])

          #if("signal" in ofile_name):
          if True:
            (lead_scalars_all, lead_isr_all) = jet_constituents(new_jet[jet_i], scalars, isrs)
            (suep_tracks,isr_tracks) = (len(scalars),len(isrs))
          else:
            (lead_scalars_all, lead_isr_all) = (0,0)
            (suep_tracks,isr_tracks) = (0,0)
          #outfile_all.write("%d %s %.1f %d %f %d %f %f %f %d %d %d %d %d %d\n"%(ievt,jet_algo_title,Ri, jet_i,
          #print("%d %s %.1f %d %f %d %f %f %f %d %d %d %d %d %d\n"%(ievt,jet_algo_title,Ri, jet_i,
          fout.write("%d %s %.1f %d %f %d %f %f %f %d %d %d %d %d %d\n"%(ievt,jet_algo_title,Ri, jet_i,
            new_jet[jet_i].pt,
            len([x for x in new_jet[jet_i].constituents_array() if x['pT'] > 1]),
            girth_all,
            new_jet[jet_i].mass,
            trackpt_all,
            lead_scalars_all, lead_isr_all,suep_tracks,isr_tracks,
            nvtx[ievt],numInteractions[ievt]
          ))

class Objects:
    def __init__(self, df,fout):
        #self.df = df
        print("Getting Ht")
        self.Ht = df["HT"]
        self.nvtx = df["NVtx"]
        self.numInteractions = df["NumInteractions"]
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
            bDeepCSVprob=df['Jets_bJetTagDeepCSVBvsAll'].flatten()
        )
        print("Getting Gen")
        self.gen = JaggedCandidateArray.candidatesfromcounts(
            df['GenParticles_Charge'].counts,
            pt=df['GenParticles.fCoordinates.fPt'].flatten(),
            eta=df['GenParticles.fCoordinates.fEta'].flatten(),
            phi=df['GenParticles.fCoordinates.fPhi'].flatten(),
            energy=df['GenParticles.fCoordinates.fE'].flatten(),
            parentId=df['GenParticles_ParentId'].flatten(),
            charge=df['GenParticles_Charge'].flatten(),
            pdgId=df['GenParticles_PdgId'].flatten(),
            status=df['GenParticles_Status'].flatten()
        )
        finalParticles = (self.gen.status == 1) & (self.gen.pt > 1) & (abs(self.gen.charge) ==1) & (abs(self.gen.eta) < 2.5)
        self.scalars = self.gen[finalParticles & (self.gen.parentId == 999998) & ((abs(self.gen.pdgId) == 11) | (abs(self.gen.pdgId) == 13) | (abs(self.gen.pdgId) == 22) | (abs(self.gen.pdgId) > 100 ))]
        self.isrs = self.gen[finalParticles & (self.gen.parentId != 999998) & ((abs(self.gen.pdgId) == 11) | (abs(self.gen.pdgId) == 13) | (abs(self.gen.pdgId) == 22) | (abs(self.gen.pdgId) > 100 ))]

       #(self.fatjets10, self.fatjets15, self.fatjets20) = makeFatJets(self.tracks)
        #for ievt in range(self.tracks.size):
        #  if ievt % 1000 == 0: print("Event {}".format(ievt))
        #  makeJets_Comparisons(-1,0.8,self.tracks[ievt])
        self.trackEtaCut = 2.5
        self.etaCut = 2.4
        trackQualityCut = (self.tracks.pt > 1) & (abs(self.tracks.eta) < self.trackEtaCut) & (self.tracks.fromPV0 >= 2 ) & ( self.tracks.matchedToPFCandidate > 0 )
        self.goodTracks = self.tracks[trackQualityCut]
        compareJets(self.goodTracks,self.scalars,self.isrs,self.Ht,self.nvtx,self.numInteractions,fout)

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
    #def compareJetAlgos(self):
    #    return compareJets(self.goodTracks,self.scalars,self.isrs,self.Ht,self.nvtx,self.numInteractions)
    
