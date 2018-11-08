// 'ChangeColorPalette.C'
// Derek Anderson
// 11.08.2018
//
// Handy way to change the color
// palette!
//
// Courtesy of: http://ultrahigh.org/2007/08/making-pretty-root-color-palettes/


#include <iostream>
#include "TStyle.h"
#include "TColor.h"

using namespace std;


// no. of colors and no. of contours
static const UInt_t NCols(5);
static const UInt_t NCont(255);



void ChangeColorPalette() {

  Double_t stops[NCols] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NCols]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NCols] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NCols]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NCols, stops, red, green, blue, NCont);
  gStyle -> SetNumberContours(NCont);

}

// End ------------------------------------------------------------------------
