import coffea

import matplotlib.pyplot as plt


plt.style.use(hep.style.ROOT)


# Get the file and import using uproot
mMed = 1000
mDark = 2
temp = 2
decayMode = 'darkPhoHad'


base = '/Users/kdipetri/Documents/Fermilab/SUEP/SUEP_coffea/inputs/'
#base = 'root://cmseos.fnal.gov//store/user/kdipetri/SUEP/Production_v0.2/2018/NTUP/'
datasets = [base +
            'PrivateSamples.SUEP_2018_mMed-%d_mDark-%d_temp-%d_decay-%s'
            '_13TeV-pythia8_n-100_0_RA2AnalysisTree.root'%(mMed, mDark, temp, decayMode),
           ]
rootfile = datasets[0]
fin = uproot.open(rootfile)


Tracks_x = get_branch('Tracks.fCoordinates.fX')
Tracks_y = get_branch('Tracks.fCoordinates.fY')
Tracks_z = get_branch('Tracks.fCoordinates.fZ')
Tracks_fromPV0 = get_branch('Tracks_fromPV0')
Tracks_matchedToPFCandidate = get_branch('Tracks_matchedToPFCandidate')
HT = get_branch(b'HT')









