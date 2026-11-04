#pragma once

// There is no residual cut on any plot, entry list or column. AnglesRaw keeps every fitted proton and the Angles view selects all of them (CutFiles.h).
// resDRms (RMS transverse distance of the track's hits to its fitted line, detector block position, mm) stays a branch of AnglesRaw (ReconAngleAll.cpp) and of every
// cut_APandT.root tree (CalcAPandT.cpp), so it can be plotted.
// The 4 mm block residual limit of the 3_ResCuts sample, a cut study of 6_cuts (CutPlots.cpp (samples), Plot3D.cpp).
static const double kRes4Mm = 4.0;
