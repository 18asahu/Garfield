import ROOT
import os, sys
import math
import ctypes

path = os.getenv('GARFIELD_INSTALL')
if sys.platform == 'darwin':
  ROOT.gSystem.Load(path + '/lib/libmagboltz.dylib')
  ROOT.gSystem.Load(path + '/lib/libGarfield.dylib')
else:
  ROOT.gSystem.Load(path + '/lib/libmagboltz.so')
  ROOT.gSystem.Load(path + '/lib/libGarfield.so')

gas = ROOT.Garfield.MediumMagboltz()
gas.LoadGasFile('ar_93_co2_7_3bar.gas')
gas.LoadIonMobility(path + '/share/Garfield/Data/IonMobility_Ar+_Ar.txt')

cmp = ROOT.Garfield.ComponentAnalyticField()
cmp.SetMedium(gas)
# Wire radius [cm]
rWire = 25.e-4
# Outer radius of the tube [cm]
rTube = 0.71
# Voltages
vWire = 2730.
vTube = 0.
# Add the wire in the centre.
cmp.AddWire(0, 0, 2 * rWire, vWire, 's')
# Add the tube.
cmp.AddTube(rTube, vTube, 0, 't')
# Request calculation of the weighting field. 
cmp.AddReadout('s')

# Make a sensor.
sensor = ROOT.Garfield.Sensor()
sensor.AddComponent(cmp);
sensor.AddElectrode(cmp, 's')
# Set the signal time window.
tstep = 0.5;
tmin = -0.5 * tstep
nbins = 1000
sensor.SetTimeWindow(tmin, tstep, nbins)
# Set the delta reponse function.
infile = open('mdt_elx_delta.txt', 'r')
times = ROOT.std.vector('double')()
values = ROOT.std.vector('double')()
for line in infile:
  line = line.strip()
  line = line.split()
  times.push_back(1.e3 * float(line[0]))
  values.push_back(float(line[1]))
infile.close()
sensor.SetTransferFunction(times, values)

# Set up Heed.
track = ROOT.Garfield.TrackHeed()
track.SetParticle('muon')
track.SetEnergy(170.e9)
track.SetSensor(sensor)

# RKF integration.
drift = ROOT.Garfield.DriftLineRKF()
drift.SetSensor(sensor)
drift.SetGainFluctuationsPolya(0., 20000.)

driftView = ROOT.Garfield.ViewDrift()
cD = ROOT.TCanvas('cD', '', 600, 600)
driftView.SetCanvas(cD)
cellView = ROOT.Garfield.ViewCell()
plotDrift = True
if plotDrift:
  drift.EnablePlotting(driftView)
  track.EnablePlotting(driftView)
  cellView.SetComponent(cmp)
  cellView.SetCanvas(driftView.GetCanvas())

signalView = ROOT.Garfield.ViewSignal()
cS = ROOT.TCanvas('cS', '', 600, 600)
signalView.SetCanvas(cS)
signalView.SetSensor(sensor)
plotSignal = True

rTrack = 0.3
x0 = rTrack
y0 = -math.sqrt(rTube * rTube - rTrack * rTrack)

nTracks = 1
for j in range(nTracks):
  sensor.ClearSignal()
  track.NewTrack(x0, y0, 0, 0, 0, 1, 0)
  xc = ctypes.c_double(0.)
  yc = ctypes.c_double(0.)
  zc = ctypes.c_double(0.)
  tc = ctypes.c_double(0.)
  ec = ctypes.c_double(0.)
  extra = ctypes.c_double(0.)
  nc = ctypes.c_int(0)
  while track.GetCluster(xc, yc, zc, tc, nc, ec, extra):
    for k in range(nc.value):
      xe = ctypes.c_double(0.)
      ye = ctypes.c_double(0.)
      ze = ctypes.c_double(0.)
      te = ctypes.c_double(0.)
      ee = ctypes.c_double(0.)
      dx = ctypes.c_double(0.)
      dy = ctypes.c_double(0.)
      dz = ctypes.c_double(0.)
      track.GetElectron(k, xe, ye, ze, te, ee, dx, dy, dz)
      drift.DriftElectron(xe.value, ye.value, ze.value, te.value)
    if plotDrift:
      driftView.GetCanvas().Clear()
      cellView.Plot2d()
      driftView.Plot(True, False)
      ROOT.gPad.Update()
  sensor.ConvoluteSignals()
  nt = ctypes.c_int(0)
  if sensor.ComputeThresholdCrossings(-2., 's', nt) == False:
    continue
  if plotSignal:
    signalView.PlotSignal('s')

