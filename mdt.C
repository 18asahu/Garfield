#include <iostream>
#include <fstream>
#include <cstdlib>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/ViewCell.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/TrackHeed.hh"

using namespace Garfield;

bool readTransferFunction(Sensor& sensor) {

  std::ifstream infile;
  infile.open("mdt_elx_delta.txt", std::ios::in);
  if (!infile) {
    std::cerr << "Could not read delta response function.\n";
    return false;
  }
  std::vector<double> times;
  std::vector<double> values;
  while (!infile.eof()) {
    double t = 0., f = 0.;
    infile >> t >> f;
    if (infile.eof() || infile.fail()) break;
    times.push_back(1.e3 * t);
    values.push_back(f);
  }
  infile.close();
  sensor.SetTransferFunction(times, values);
  return true;
}

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
 
  // Make a gas medium.
  MediumMagboltz gas;
  gas.LoadGasFile("ar_93_co2_7_3bar.gas");
  const std::string path = std::getenv("GARFIELD_INSTALL");
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
  
  std::cout<<"~~~~~~~~~~~~~CHECK 0: Setting gas properties now: ~~~~~~~~~~~~~~" <<std::endl;
  gas.SetComposition("ar", 93., "co2", 7.);
  gas.SetPressure(3 * 760.);
  gas.SetTemperature(293.15);

  std::cout<<"Specifying number of electric fields" <<std::endl;
  gas.SetFieldGrid(100., 100.e3, 20, true);
//20 field points, between 100V/cm and 100kV/cm, logarithmic spacing

  const int ncoll =10;
  gas.GenerateGasTable(ncoll);
  //no. of collisions in multiples of 10^7 over which electron is traced by magboltz

  gas.WriteGasFile("ar_93_co2_7.gas");
  gas.LoadGasFile("ar_93_co2_7.gas");

std::cout << "Check 1" <<std::endl;


  // Make a component with analytic electric field.
  ComponentAnalyticField cmp;
  cmp.SetMedium(&gas);
  // Wire radius [cm]
  const double rWire = 25.e-4;
  // Outer radius of the tube [cm]
  const double rTube = 0.71;
  // Voltages
  const double vWire = 2730.;
  const double vTube = 0.;
  // Add the wire in the centre.
  cmp.AddWire(0, 0, 2 * rWire, vWire, "s");
  // Add the tube.
  cmp.AddTube(rTube, vTube, 0, "t");
  // Request calculation of the weighting field. 
  cmp.AddReadout("s");
  std::cout << "Check 2" <<std::endl;

  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&cmp, "s");
  std::cout <<"############################### Check 3: SENSOR MADE" <<std::endl; 

 // Set the signal time window.
  const double tstep = 0.5;
  const double tmin = -0.5 * tstep;
  const unsigned int nbins = 1000;
  sensor.SetTimeWindow(tmin, tstep, nbins);

  std::cout <<"############################## Check 4: SIGNAL TIME WINDOW SET" <<std::endl;

  // Set the delta reponse function.
  if (!readTransferFunction(sensor)) return 0;
  sensor.ClearSignal();
  std::cout <<"################################## Check 5" <<std::endl;

  // Set up Heed.
  TrackHeed track;
  track.SetParticle("muon");
  track.SetEnergy(170.e9);
  track.SetSensor(&sensor);
  std::cout <<"Check 6" <<std::endl;


  // RKF integration.
  DriftLineRKF drift;
  drift.SetSensor(&sensor);
  drift.SetGainFluctuationsPolya(0., 20000.);
  std::cout <<"Check 7: gas integration complete" <<std::endl;
 
  TCanvas* cD = nullptr;
  ViewCell cellView;
  ViewDrift driftView;
  constexpr bool plotDrift = true;
  if (plotDrift) {
    cD = new TCanvas("cD", "", 600, 600);
    cellView.SetCanvas(cD);
    cellView.SetComponent(&cmp);
    driftView.SetCanvas(cD);
    drift.EnablePlotting(&driftView);
    track.EnablePlotting(&driftView);
  }
 
std::cout <<"~~~~~~~~CHECK 8: TCANVAS STUFF ONGOING~~~~~~~~" <<std::endl;

  TCanvas* cS = nullptr;
  ViewSignal signalView;
  constexpr bool plotSignal = true;
  if (plotSignal) {
    cS = new TCanvas("cS", "", 600, 600);
    signalView.SetCanvas(cS);
    signalView.SetSensor(&sensor);
    signalView.SetLabelY("signal [fC]");
  } 

  const double rTrack = 0.3;
  const double x0 = rTrack;
  const double y0 = -sqrt(rTube * rTube - rTrack * rTrack);
  const unsigned int nTracks = 1;
  for (unsigned int j = 0; j < nTracks; ++j) {
    sensor.ClearSignal();
    track.NewTrack(x0, y0, 0, 0, 0, 1, 0);
    double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
    int nc = 0;
    while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      for (int k = 0; k < nc; ++k) {
        double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
        double dx = 0., dy = 0., dz = 0.;
        track.GetElectron(k, xe, ye, ze, te, ee, dx, dy, dz);
        drift.DriftElectron(xe, ye, ze, te);
      }
    }

std::cout <<"~~~~~~~~~~~~~~~~CHECK 9: TCANVAS STUFF DONE~~~~~~" <<std::endl;

    if (plotDrift) {
      cD->Clear();
      cellView.Plot2d();
      constexpr bool twod = true;
      constexpr bool drawaxis = false;
      driftView.Plot(twod, drawaxis);
    }
    sensor.ConvoluteSignals();
    int nt = 0;
    if (!sensor.ComputeThresholdCrossings(-2., "s", nt)) continue;
    if (plotSignal) signalView.PlotSignal("s");
  }
std::cout <<"~~~~~~~~~~~~~~~~CHECK 10 ~~~~~~~~~~~~~~~~~~~~~~" <<std::endl;
  app.Run(kTRUE);
std::cout <<"~~~~~~~~~~~~~~~~CHECK 11~~~~~~~~~~~~~~~~~~~~~~~" <<std::endl;
}
