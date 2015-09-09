// Helper class to set the style of histograms.
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TH1.h>
#include <TGraph.h>
#include "HistoStyler.h"

ClassImp(HistoStyler);

const Color_t HistoStyler::fgkColors[HistoStyler::fgkNColors] = {kBlue+1, kRed+1, kGreen+2, kMagenta+1, kCyan+3, kOrange+2, kViolet, kYellow+3, kPink+2, kTeal+2};
const Style_t HistoStyler::fgkFullMarkerStyles[HistoStyler::fgkNFullMarkerStyles] = {kFullCircle, kFullSquare, kFullDiamond, kFullStar, kFullCross, kFullTriangleUp, kFullTriangleDown};
const Style_t HistoStyler::fgkOpenMarkerStyles[HistoStyler::fgkNOpenMarkerStyles] = {kOpenCircle, kOpenSquare, kOpenDiamond, kOpenStar, kOpenCross, kOpenTriangleUp, kOpenTriangleDown};

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
void HistoStyler::Apply(Int_t n, TObject** obj) const
{
  // Apply style to histograms.

  for (Int_t i = 0; i < n; i++) {
    Apply(obj[i], i, i);
  }
}

//____________________________________________________________________________________
void HistoStyler::Apply(Int_t n, TH1** obj) const
{
  // Apply style to histograms.

  for (Int_t i = 0; i < n; i++) {
    Apply(obj[i], i, i);
  }
}

//____________________________________________________________________________________
void HistoStyler::Apply(Int_t n, TGraph** obj) const
{
  // Apply style to histograms.

  for (Int_t i = 0; i < n; i++) {
    Apply(obj[i], i, i);
  }
}

//____________________________________________________________________________________
void HistoStyler::Apply(TObject* obj, Int_t icolor, Int_t imarker) const
{
  // Apply style to histograms.

  ApplyLine(dynamic_cast<TAttLine*>(obj), icolor);
  ApplyFill(dynamic_cast<TAttFill*>(obj), icolor);
  ApplyMarker(dynamic_cast<TAttMarker*>(obj), icolor, imarker);
}

//____________________________________________________________________________________
void HistoStyler::ApplyMarker(TAttMarker* obj, Int_t icolor, Int_t imarker) const
{
  // Apply style to histograms.
  if (!obj) return;

  obj->SetMarkerColor(GetMarkerColor(icolor));

  Style_t markerStyle = GetMarkerStyle(imarker);
  Size_t markerSize = GetMarkerSize();

  if (markerStyle == kFullTriangleUp || markerStyle == kOpenTriangleUp || markerStyle == kFullTriangleDown || markerStyle == kOpenTriangleDown) {
    markerSize *= 1.5;
  }

  if (markerStyle == kFullDiamond || markerStyle == kOpenDiamond || markerStyle == kFullStar || markerStyle == kOpenStar
      || markerStyle == kFullCross || markerStyle == kOpenCross) {
    markerSize *= 2.0;
  }

  obj->SetMarkerSize(markerSize);
  obj->SetMarkerStyle(markerStyle);
}

//____________________________________________________________________________________
void HistoStyler::ApplyLine(TAttLine* obj, Int_t icolor) const
{
  // Apply style to histograms.
  if (!obj) return;

  obj->SetLineColor(GetLineColor(icolor));
  obj->SetLineStyle(GetLineStyle());
}

//____________________________________________________________________________________
void HistoStyler::ApplyFill(TAttFill* obj, Int_t icolor) const
{
  // Apply style to histograms.
  if (!obj) return;

  obj->SetFillColor(GetFillColor(icolor));
  obj->SetFillStyle(GetFillStyle());
}

//____________________________________________________________________________________
Color_t HistoStyler::GetMarkerColor(Int_t i) const
{
  switch (fVariableMarkerColor) {
  case kMarkerColorFixed:
    {
      return GetMarkerColor();
    }
  case kMarkerColorVariable:
    {
      if (i < 0) {
        return GetMarkerColor();
      }
      else if (i < fgkNColors) {
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
Style_t HistoStyler::GetMarkerStyle(Int_t i) const
{
  if (i < 0) return GetMarkerStyle();
  
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
Color_t HistoStyler::GetLineColor(Int_t i) const
{
  switch (fVariableLineColor) {
  case kLineColorFixed:
    {
      return GetLineColor();
    }
  case kLineColorVariable:
    {
      if (i < 0) {
        return GetLineColor();
      }
      else if (i < fgkNColors) {
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
Style_t HistoStyler::GetLineStyle(Int_t /*i*/) const
{
  return GetLineStyle();
}

//____________________________________________________________________________________
Color_t HistoStyler::GetFillColor(Int_t i) const
{
  switch (fVariableFillColor) {
  case kFillColorFixed:
    {
      return GetLineColor();
    }
  case kFillColorVariable:
    {
      if (i < 0) {
        return GetFillColor();
      }
      else if (i < fgkNColors) {
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
Style_t HistoStyler::GetFillStyle(Int_t /*i*/) const
{
  return GetFillStyle();
}
