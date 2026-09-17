"""Sanity checks on the IDEA o2 reconstruction output.

Run by test_IDEAo2_reco.sh once the reco job has finished.  Covers the behaviour
the seeded-clustering chain adds: optical digitization, the truth links it
produces (including the Cherenkov masking), and that clusters are grown at all.
"""

import sys

import podio.reading


def check(path):
    nClusters = nCells = nSignal = nUnlinked = nCherenLinks = 0
    cherenMasked = True

    for frame in podio.reading.get_reader(path).get("events"):
        nClusters += len(frame.get("TopoGrownClusters"))

        # the digitizer emits exactly one cell per optical hit that carries signal
        nSignal += sum(
            1 for hit in frame.get("SCEPCal_MainScounts") if hit.getEnergy() > 0
        )
        cells = frame.get("SCEPCal_digi_scint")
        nCells += len(cells)

        # every digitized cell must carry a truth link to its energy deposit
        scint = frame.get("SCEPCal_scint_link")
        linked = {link.getFrom().getCellID() for link in scint}
        nUnlinked += sum(1 for cell in cells if cell.getCellID() not in linked)

        # maskCherenkovForTruthLink: Cherenkov cells resolve to the same deposits as the
        # scintillation cells.  Without the mask this collection would come out empty.
        cheren = frame.get("SCEPCal_cheren_link")
        nCherenLinks += len(cheren)
        cherenMasked &= {link.getTo().getCellID() for link in cheren} <= {
            link.getTo().getCellID() for link in scint
        }

    print(
        f"clusters {nClusters}, cells {nCells}/{nSignal}, "
        f"unlinked cells {nUnlinked}, cherenkov links {nCherenLinks}"
    )

    checks = [
        (nClusters > 0, "no topo clusters were grown"),
        (nCells == nSignal, "one cell per optical hit with signal was not produced"),
        (nUnlinked == 0, f"{nUnlinked} digitized cells carry no truth link"),
        (
            nCherenLinks > 0 and cherenMasked,
            "cherenkov cells are not linked to the scintillation deposits",
        ),
    ]
    for ok, msg in checks:
        if not ok:
            print("FAIL:", msg)

    return all(ok for ok, _ in checks)


if __name__ == "__main__":
    sys.exit(0 if check(sys.argv[1]) else 1)
