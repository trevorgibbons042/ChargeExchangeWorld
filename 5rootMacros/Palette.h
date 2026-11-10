#pragma once
#include <TColor.h>

// Color blind friendly colors (Okabe and Ito), the same three in every plot of the analysis:
//   CbCE()      CE recon protons                  vermillion   (213, 94, 0)
//   CbProton()  all recon protons                 blue         (0, 114, 178)
//   CbAll()     all recon charged particles       bluish green (0, 158, 115)
// They stand in for ROOT's kRed, kBlue and kGreen + 2. Functions with a function-local static, not globals: a namespace-scope global with a runtime initializer
// fails on the cluster (5rootMacros/CLAUDE.md, Failed to materialize symbols).
inline int CbCE(){static const int c = TColor::GetColor(213, 94, 0); return c;}
inline int CbProton(){static const int c = TColor::GetColor(0, 114, 178); return c;}
inline int CbAll(){static const int c = TColor::GetColor(0, 158, 115); return c;}
