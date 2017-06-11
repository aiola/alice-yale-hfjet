#!/usr/local/bin/python
#python base classes and utilities for D Meson jet analysis

import enum

class AnalysisType(enum.Enum):
    InvMassFit = 1
    SideBand = 2
    LikeSign = 3
    LikeSignFit = 4
    Truth = 5
