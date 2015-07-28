// Helper class to set the style of histograms.
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TH1.h>
#include "HistoStyler.h"

ClassImp(HistoStyler);

const Color_t HistoStyler::fgkColors[HistoStyler::fgkNColors] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1, kCyan+3, kOrange+2, kViolet, kYellow+3, kPink+2, kTeal+2};
const Style_t HistoStyler::fgkFullMarkerStyles[HistoStyler::fgkNFullMarkerStyles] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullStar};
const Style_t HistoStyler::fgkOpenMarkerStyles[HistoStyler::fgkNOpenMarkerStyles] = {kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown, kOpenDiamond, kOpenStar};

//____________________________________________________________________________________
HistoStyler::HistoStyler() :
  TAttLine(),
  TAttFill(),
  TAttMarker(),
  fVariableMarkerColor(kMarkerColorFixed),
  fVariableMarkerStyle(kMarkerStyleFixed),
  fVariableLineColor(kLineColorFixed),
  fVariableLineStyle(kLineStyleFixed),
  fVariableFillColor(kFillColorFixed),
  fVariableFillStyle(kFillStyleFixed)
{
  // Default constructor.
}

//____________________________________________________________________________________
void HistoStyler::Apply(Int_t n, TH1** histo) const
{
  // Apply style to histograms.

  for (Int_t i = 0; i < n; i++) {
    if (!histo[i]) continue;
    
    histo[i]->SetMarkerColor(GetMarkerColor(i));
    histo[i]->SetMarkerStyle(GetMarkerStyle(i));
    histo[i]->SetMarkerSize(GetMarkerSize());

    histo[i]->SetLineColor(GetLineColor(i));
    histo[i]->SetLineStyle(GetLineStyle(i));

    histo[i]->SetFillColor(GetFillColor(i));
    histo[i]->SetFillStyle(GetFillStyle(i));
  }
}

//____________________________________________________________________________________
Color_t HistoStyler::GetMarkerColor(UShort_t i) const
{
  switch (fVariableMarkerColor) {
  case kMarkerColorFixed:
    {
      return GetMarkerColor();
    }
  case kMarkerColorVariable:
    {
      if (i < fgkNColors) {
        return fgkColors[i];
      }
      else {
        Printf("Warning: color index %d overflow!", i);
        return fgkColors[fgkNColors-1];
      }
    }
  default:  // should not happen
    {
      return GetMarkerColor();
    }
  }
}

//____________________________________________________________________________________
Style_t HistoStyler::GetMarkerStyle(UShort_t i) const
{
  switch (fVariableMarkerStyle) {
  case kMarkerStyleFixed:
    {
      return GetMarkerStyle();
    }
  case kMarkerStyleVariable:
    {
      
      Int_t p = i % 2;
      if (p == 0) {
        UShort_t j = i / 2;
        if (j < fgkNFullMarkerStyles) {
          return fgkFullMarkerStyles[j];
        }
        else {
          Printf("Warning: marker style index %d overflow!", j);
          return fgkFullMarkerStyles[fgkNFullMarkerStyles-1];
        }
      }
      else {
        UShort_t j = (i-1) / 2;
        if (j < fgkNOpenMarkerStyles) {
          return fgkOpenMarkerStyles[j];
        }
        else {
          Printf("Warning: marker style index %d overflow!", j);
          return fgkOpenMarkerStyles[fgkNOpenMarkerStyles-1];
        }
      }
    }
    case kMarkerStyleOpen:
      {
        if (i < fgkNOpenMarkerStyles) {
          return fgkOpenMarkerStyles[i];
        }
        else {
          Printf("Warning: marker style index %d overflow!", i);
          return fgkOpenMarkerStyles[fgkNOpenMarkerStyles-1];
        }
      }
    case kMarkerStyleFull:
      {
        if (i < fgkNFullMarkerStyles) {
          return fgkFullMarkerStyles[i];
        }
        else {
          Printf("Warning: marker style index %d overflow!", i);
          return fgkFullMarkerStyles[fgkNFullMarkerStyles-1];
        }
      }
  default:  // should not happen
    {
      return GetMarkerStyle();
    }
  }
}

//____________________________________________________________________________________
Color_t HistoStyler::GetLineColor(UShort_t i) const
{
  switch (fVariableLineColor) {
  case kLineColorFixed:
    {
      return GetLineColor();
    }
  case kLineColorVariable:
    {
      if (i < fgkNColors) {
        return fgkColors[i];
      }
      else {
        Printf("Warning: color index %d overflow!", i);
        return fgkColors[fgkNColors-1];
      }
    }
  default:  // should not happen
    {
      return GetLineColor();
    }
  }
}

//____________________________________________________________________________________
Style_t HistoStyler::GetLineStyle(UShort_t /*i*/) const
{
  return GetLineStyle();
}

//____________________________________________________________________________________
Color_t HistoStyler::GetFillColor(UShort_t i) const
{
  switch (fVariableFillColor) {
  case kFillColorFixed:
    {
      return GetLineColor();
    }
  case kFillColorVariable:
    {
      if (i < fgkNColors) {
        return fgkColors[i];
      }
      else {
        Printf("Warning: color index %d overflow!", i);
        return fgkColors[fgkNColors-1];
      }
    }
  default:  // should not happen
    {
      return GetFillColor();
    }
  }
}

//____________________________________________________________________________________
Style_t HistoStyler::GetFillStyle(UShort_t /*i*/) const
{
  return GetFillStyle();
}
