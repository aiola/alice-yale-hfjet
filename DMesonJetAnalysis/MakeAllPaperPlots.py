#!/usr/bin/env python
# python script to do extract B feed down correction factors

import IPython
import BFeedDown_Paper
import Efficiency_Paper
import JetPtSpectrum_Paper
import JetZSpectrum_JetPt_5_15_Paper
import JetZSpectrum_JetPt_15_30_Paper
import Uncertainties_JetPtSpectrum_Paper


def main():
    BFeedDown_Paper.main()
    Efficiency_Paper.main()
    JetPtSpectrum_Paper.main()
    JetZSpectrum_JetPt_5_15_Paper.main()
    JetZSpectrum_JetPt_15_30_Paper.main()
    Uncertainties_JetPtSpectrum_Paper.main()


if __name__ == '__main__':
    main()

    IPython.embed()
