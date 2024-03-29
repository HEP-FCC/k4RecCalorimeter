from ROOT import gSystem
from EventStore import EventStore

store = EventStore(['output_ecalReco_noiseFromFile_test.root'])

for event in store:
    clusters = event.get('caloClusters')
    assert(len(clusters) == 1)
