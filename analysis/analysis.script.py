classLoader = UserClassLoader(path, sys.getClassLoader())
QuakePredictor = classLoader.findClass("QuakePredictor")
qp = QuakePredictor();
dataFolder = "/Volumes/Time Machine Backups/temp/"
seed = 4
h = 0
info = qp.loadEQEventInfo(dataFolder, seed)
h0_data = qp.loadEQHourlyInfo(dataFolder, seed, h)
from jhplot import *;from jhplot.stat import LinReg;from jarray import array;from java.lang import Double;import math
c1 = HPlot(); c1.setAutoRange(); c1.visible(); c1.clear()
pf_0_0 = P1D("site 0 channel 0")
pf_0_0.fill([ x / 1.0 for x in range(h0_data.sitesData.get(0).data[0].__len__())], h0_data.sitesData.get(0).data[0].tolist())
pf_0_0.setSymbolSize(0)
c1.clear(); c1.setRangeX(52000, 103000); c1.setRangeY(8500000.0, 8600000.0); c1.setMarginLeft(100); c1.setMarginTop(100)
c1.draw(pf_0_0)

pf_0_1 = P1D("site 0 channel 1")
pf_0_1.fill([ x / 1.0 for x in range(h0_data.sitesData.get(0).data[1].__len__())], h0_data.sitesData.get(0).data[1].tolist())
pf_0_1.setSymbolSize(0)
c1.clear(); c1.setRangeX(52000, 103000); c1.setRangeY(8500000.0, 8600000.0); c1.setMarginLeft(100); c1.setMarginTop(100)
c1.draw(pf_0_1)

pf_0_2 = P1D("site 0 channel 2")
pf_0_2.fill([ x / 1.0 for x in range(h0_data.sitesData.get(0).data[2].__len__())], h0_data.sitesData.get(0).data[2].tolist())
pf_0_2.setSymbolSize(0)
c1.clear(); c1.setRangeX(52000, 103000); c1.setRangeY(8500000.0, 8600000.0); c1.setMarginLeft(100); c1.setMarginTop(100)
c1.draw(pf_0_1)

# Show hourly magnetic field profiles of 3 channels and relevant properties.
def show3Profiles(seed, site, h):
    
    # Load earthquake relevant information.
    info = qp.loadEQEventInfo(dataFolder, seed)
    h_data = qp.loadEQHourlyInfo(dataFolder, seed, h)
    
    # Print earthquake information.
    print "Seed: " + str(seed)
    print "Earthquake position: " + str(info.EQEventLatLong.latitude) + ", " + str(info.EQEventLatLong.longitude)
    print "Hour: " + str(info.EQEventHour)
    print "Magnitude: " + str(info.EQEventMagnitude)
    print "The closest site: " + str(info.EQEventSite)
    print "Distance from the closest site: " + str(info.EQEventDistToQuake)
    
    # Print site relevant information.
    print " "
    print "Measurement site: " + str(site) + "/" + str(info.sitesInfo.sites.size())
    print "Sampling rate: " + str(info.sitesInfo.sampleRate)
    print "Position: " + str(info.sitesInfo.sites.get(site).latitude) + ", " + str(info.sitesInfo.sites.get(site).longitude)
    print "Distance from the epicenter: " + str(QuakePredictor.calEarthDistance(info.EQEventLatLong, info.sitesInfo.sites.get(site)))
    print "Transit time from the epicenter: " + str(QuakePredictor.calEQTransitTime(info.EQEventLatLong, info.sitesInfo.sites.get(site)))
    
    # Show profiles of 3 channels.
    print " "
    print "Hour: " + str(h)
    
    c1 = HPlot("site " + str(site) + " hour " + str(h)); c1.visible(); c1.clear(); c1.setAutoRange(); c1.setMarginLeft(70); c1.setMarginTop(50)

    pf_0 = P1D("site " + str(site) + " channel 0")
    pf_0.fill([ x / 1.0 for x in range(h_data.sitesData.get(site).data[0].__len__())], h_data.sitesData.get(site).data[0].tolist())
    pf_0.setSymbol(11); pf_0.setSymbolSize(1); pf_0.setDrawLine(True); pf_0.setDrawSymbol(False);pf_0.setColor(Color.RED)

    pf_1 = P1D("site " + str(site) + " channel 1")
    pf_1.fill([ x / 1.0 for x in range(h_data.sitesData.get(site).data[1].__len__())], h_data.sitesData.get(site).data[1].tolist())
    pf_1.setSymbol(11); pf_1.setSymbolSize(1); pf_1.setDrawLine(True); pf_1.setDrawSymbol(False);pf_1.setColor(Color.GREEN)
    
    pf_2 = P1D("site " + str(site) + " channel 2")
    pf_2.fill([ x / 1.0 for x in range(h_data.sitesData.get(site).data[2].__len__())], h_data.sitesData.get(site).data[2].tolist())
    pf_2.setSymbol(11); pf_2.setSymbolSize(1); pf_2.setDrawLine(True); pf_2.setDrawSymbol(False);pf_2.setColor(Color.BLUE)
    
    c1.draw(pf_0); c1.draw(pf_1); c1.draw(pf_2);
    
dSet_4 = qp.loadDataSet(dataFolder, 4)
p_k4 = P1D("Seed 4: air conductivity")
p_k4.fill([x / 1.0 for x in range(dSet_4.EMA.__len__())], dSet_4.EMA.tolist()); p_k4.setDrawLine(True); p_k4.setSymbol(11); p_k4.setSymbolSize(0)
c2 = HPlot("Seed 4: Air conductivity"); c2.clear(); c2.visible(); c2.setAutoRange(); c2.setMarginLeft(70); c2.setMarginTop(50)

dSet_6 = qp.loadDataSet(dataFolder, 6)
p_k6 = P1D("Seed 6: air conductivity")
p_k6.fill([x / 1.0 for x in range(dSet_6.EMA.__len__())], dSet_6.EMA.tolist()); p_k6.setDrawLine(True); p_k6.setSymbol(11); p_k6.setSymbolSize(0)
c6 = HPlot("Seed 6: Air conductivity"); c6.clear(); c6.visible(); c6.setAutoRange(); c6.setMarginLeft(70); c6.setMarginTop(50)
c6.draw(p_k6)

dSet_7 = qp.loadDataSet(dataFolder, 7)
p_k7 = P1D("Seed 7: air conductivity")
p_k7.fill([x / 1.0 for x in range(dSet_7.EMA.__len__())], dSet_7.EMA.tolist()); p_k7.setDrawLine(True); p_k7.setSymbol(11); p_k7.setSymbolSize(0)
c7 = HPlot("Seed 7: Air conductivity"); c7.clear(); c7.visible(); c7.setAutoRange(); c7.setMarginLeft(70); c7.setMarginTop(50)
c7.draw(p_k7)

dSet_8 = qp.loadDataSet(dataFolder, 8)
p_k8 = P1D("Seed 8: air conductivity")
p_k8.fill([x / 1.0 for x in range(dSet_8.EMA.__len__())], dSet_8.EMA.tolist()); p_k8.setDrawLine(True); p_k8.setSymbol(11); p_k8.setSymbolSize(0)
c8 = HPlot("Seed 8: Air conductivity"); c8.clear(); c8.visible(); c8.setAutoRange(); c8.setMarginLeft(70); c8.setMarginTop(50)
c8.draw(p_k8)

      
