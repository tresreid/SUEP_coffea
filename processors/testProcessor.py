from coffea import hist, processor
import numpy as np
import awkward1 as ak
from utils import utility as utl
from utils.objects import Objects
from utils import baseline as bl
from utils import fatjets as fj

class TestProcessor(processor.ProcessorABC):
        def __init__(self):
                # dataset_axis = hist.Cat("dataset", "Primary dataset")
                # pt_axis = hist.Bin("pt", r"$p_{T}$ [GeV]", 40, 0, 3500)
                #ntrack_axis                 = hist.Bin("ntracks",               "Number of Tracks",                                 20,     0.0,    500.0)
                #njet_axis                   = hist.Bin("njets",                 "Number of Jets",                                   20,     0.0,    20.0)
                #ht_axis                     = hist.Bin("ht",                    "H_{T} (GeV)",                                   500,    0.0,    5000.0)
                #st_axis                     = hist.Bin("st",                    "S_{T} (GeV)",                                   500,    0.0,    5000.0)
                #met_axis                    = hist.Bin("MET",                   "E_{T}^{miss} (GeV)",                                   200,    0.0,    2000.0)


                self._accumulator = processor.dict_accumulator({
                        # 'jtpt':hist.Hist("Counts", dataset_axis, pt_axis),
                        # 'jteta':hist.Hist("Counts",dataset_axis,eta_axis),
                        #'h_ntracks':                hist.Hist("h_ntracks",              ntrack_axis),
                        #'h_njets':                  hist.Hist("h_njets",                njet_axis),
                        #'h_ht':                     hist.Hist("h_ht",                   ht_axis),
                        #'h_st':                     hist.Hist("h_st",                   st_axis),
                        #'h_met':                    hist.Hist("h_met",                  met_axis),
                        #'cutflow':                  processor.defaultdict_accumulator(int),
                        #'trigger':                  processor.defaultdict_accumulator(int),
                        'nTracks':                  processor.defaultdict_accumulator(int),
                        #'txt':                  processor.defaultdict_accumulator(int),
                })

        @property
        def accumulator(self):
                return self._accumulator

        def process(self, df):
                output = self.accumulator.identity()

                fout=open("test3.txt","w")
                fout.write("Event Jet_algo R jet_id pt_l ntracks_l girth_l mass_l trackpt_l suep_tracks_l isr_tracks_l total_suep total_isr NVtx NumInteractions\n")
                # important variables
                obj = Objects(df,fout)
                mask = bl.getBaselineMask(df)
                #output["txt"] = obj.compareJetAlgos()
                fout.close()

                return output

        def postprocess(self, accumulator):
                return accumulator
