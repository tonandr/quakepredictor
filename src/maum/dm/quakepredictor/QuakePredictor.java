package maum.dm.quakepredictor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import maum.dm.Matrix;
import maum.dm.SparkNeuralNetwork;

/**
 * @author Inwoo Chung (gutomitai@gmail.com)
 */
public class QuakePredictor {

	public static boolean DEBUG = true;
		
	// Constants.
	public static int NUM_CHANNELS = 3;
	public static int GLOBAL_QUAKES_LENGTH = 5;
	public static int PREDICTION_HOUR_TIME_LENGTH = 2160; // 90 days.
	public static int PREDICTION_START_HOUR = 768;
	public static int NUM_M_FACTOR = 360;
	public static double K_FACTOR_NORM_NUM = 10.0;
	public static double M_FACTOR_NORM_NUM = 8500000.0;
	
	public static double LEARNING_RATE = 0.005; //?
	public static double NUM_ITERATION = 250; //?
	
	// Magnetic field values.
	public static int NUM_FACTORS = NUM_M_FACTOR;
	
	// Latitude and longitude position.
	public class LatLong {
		public double latitude;
		public double longitude;
		
		public LatLong() {
		}
		
		public LatLong(double latitude, double longitude) {
			this.latitude = latitude;
			this.longitude = longitude;
		}
	}
	
	// Sites information.
	public class SitesInfo {
		public int sampleRate;
		public Map<Integer, LatLong> sites = new HashMap<Integer, LatLong>();
		
		// Parse sites information from input data.
		public void parse(int sampleRate, int numOfSites, double[] sitesData) {
			
			// Input arguments are assumed to be valid.
			
			// Initialize sites information.
			sites.clear();
			
			// Parse.
			this.sampleRate = sampleRate;
			
			for (int i = 0; i < numOfSites; i++) {
				sites.put(i, new LatLong(sitesData[i * 2], sitesData[i * 2 + 1]));
			}
		}
	}
	
	// Earthquake event information.
	public class EQEventInfo {
		public SitesInfo sitesInfo = new SitesInfo();
		public int EQEventSite;
		public int EQEventHour;
		public LatLong EQEventLatLong = new LatLong();
		public double EQEventMagnitude;
		public double EQEventDistToQuake; // km unit.
	}
	
	// Global earthquake event information.
	public class GlobalEQEventInfo {
		public LatLong EQEventLatLong = new LatLong();
		public double EQEventDepth;
		public double EQEventMagnitude;
		public double EQEventTime; // Time in seconds from the start of a test case when an earthquake happened.
	}
	
	// HourlyMagnetometer data of a site.
	public class HourlyMagnetData {
		public int[][] data = new int[NUM_CHANNELS][];
	}
	
	// Earthquake hourly information.
	public class EQHourlyInfo {
		public int numOfSites;
		public int sampleRates;
		public int hour; // Zero-indexed hour value from the start time of magnetometer data measurement.
		public double K; // Planetary magnetic activity from 0 to 10.
		public Map<Integer, HourlyMagnetData> sitesData = new HashMap<Integer, HourlyMagnetData>();
		public int numGlobalQuakes;
		public Map<Integer, GlobalEQEventInfo> globalEQs = new HashMap<Integer, GlobalEQEventInfo>();
		
		// Parse.
		public void parse(int numOfSites, int sampleRates, int hour, int[] data, double K, double[] globalQuakes) {
			
			// Input arguments are assumed to be valid.
			
			// Initialize.
			this.numOfSites = numOfSites;
			this.sampleRates = sampleRates;
			this.hour = hour;
			this.K = K;
			sitesData.clear();
			globalEQs.clear();
			
			// Parse.
			// Magnetometer data.
			int channelDataLength = this.sampleRates * 3600; 
			
			for (int j = 0; j < numOfSites; j++) {
				HourlyMagnetData hourlyMagnetData = new HourlyMagnetData();
				
				for (int c = 0; c < NUM_CHANNELS; c++) {
					hourlyMagnetData.data[c] = new int[channelDataLength];
					
					for (int i = 0; i < channelDataLength; i++) {
						
						// If a magnetometer data value is invalid, -1 is assigned.
						hourlyMagnetData.data[c][i] 
								= data[j * channelDataLength * NUM_CHANNELS + c * channelDataLength + i]; 
					}
				}
				
				sitesData.put(j, hourlyMagnetData);
			}
			
			// Global quakes.
			numGlobalQuakes = globalQuakes.length / GLOBAL_QUAKES_LENGTH;
			
			for (int i = 0; i < numGlobalQuakes; i++) {
				GlobalEQEventInfo globalEQEventInfo = new GlobalEQEventInfo();
				
				globalEQEventInfo.EQEventLatLong.latitude = globalQuakes[i * 5];
				globalEQEventInfo.EQEventLatLong.longitude = globalQuakes[i * 5 + 1];
				globalEQEventInfo.EQEventDepth = globalQuakes[i * 5 + 2];
				globalEQEventInfo.EQEventMagnitude = globalQuakes[i * 5 + 3];
				globalEQEventInfo.EQEventTime = globalQuakes[i * 5 + 4];
				
				globalEQs.put(i, globalEQEventInfo);
			}
		}
	}
	
	// Calculate distance on the Earth surface in km.
	public static double calEarthDistance(LatLong p1, LatLong p2) {
		
		// Input arguments are assumed to be valid.
		
		// Calculate.
		final double earthRadius = 6371.01;
		
		double deltaLong = Math.abs(p1.longitude - p2.longitude);
		
		if (deltaLong > 180.0) 
			deltaLong = 360.0 - deltaLong;
		
		// Vincenty formula.
		double dist = earthRadius * Math.atan2(
				Math.sqrt( Math.pow((Math.cos(Math.toRadians(p1.latitude)) * Math.sin(Math.toRadians(deltaLong))), 2.0) 
						+ Math.pow((Math.cos(Math.toRadians(p2.latitude)) * Math.sin(Math.toRadians(p1.latitude)) 
								- Math.sin(Math.toRadians(p2.latitude)) * Math.cos(Math.toRadians(p1.latitude)) 
								* Math.cos(Math.toRadians(deltaLong))), 2.0)), 
								Math.sin(Math.toRadians(p2.latitude)) * Math.sin(Math.toRadians(p1.latitude)) 
								+ Math.cos(Math.toRadians(p2.latitude)) * Math.cos(Math.toRadians(p1.latitude)) 
								* Math.cos(Math.toRadians(deltaLong)));
		
		return dist;
	}
	
	// Calculate earthquake transit time.
	public static double calEQTransitTime(LatLong pEQ, LatLong pTarget) {
		
		// Input arguments are assumed to be valid.
		
		// Calculate.
		final double earthRadius = 6372.795;
		final double highSlope = 60.0 / 122.0;
		final double lowSlope = 53.3 / 120.0;
		double meanSlope = 60.0 * (highSlope + lowSlope) / 2.0;
		double earthCircumference = 2.0 * Math.PI * earthRadius;
		
		double earthDistance = calEarthDistance(pEQ, pTarget);
		double transitTime = (earthDistance / earthCircumference) * 360.0 * meanSlope; // Second?
		
		return transitTime;
	}
	
	// For analysis.
	// Load EQ event information.
	public EQEventInfo loadEQEventInfo(String dataFolder, int seed) throws Exception {
		
		// Input arguments are assumed to be valid.
		
		// Load.
		QuakeTester qt = new QuakeTester();
		qt.seed = seed;
		
		QuakeTester.DataSet dset = qt.new DataSet(dataFolder+seed+"/"); //?
        double[] sitesData = new double[dset.numOfSites*2];
        
        for (int i=0;i<dset.numOfSites;i++) {
           sitesData[i*2] = dset.sites[i].latitude;
           sitesData[i*2+1] = dset.sites[i].longitude;
        }
		
        // load ground truth data
        dset.readGTF();
        
        // Create EQ event information.
        EQEventInfo info = new EQEventInfo();
        
        info.sitesInfo.parse(dset.sampleRate, dset.numOfSites, sitesData);
        info.EQEventLatLong.latitude = dset.gtfLatitude;
        info.EQEventLatLong.longitude = dset.gtfLongitude;
        info.EQEventHour = dset.gtfEQHour;
        info.EQEventSite = dset.gtfSite;
        info.EQEventMagnitude = dset.gtfMagnitude;
        info.EQEventDistToQuake = dset.gtfDistToEQ;
        		
        return info;
	}
	
	// Load EQ hourly information.
	public EQHourlyInfo loadEQHourlyInfo(String dataFolder, int seed, int h) 
			throws Exception {
		
		// Input arguments are assumed to be valid.
		
		// Load.
		QuakeTester qt = new QuakeTester();
		qt.seed = seed;
		
		QuakeTester.DataSet dset = qt.new DataSet(dataFolder+seed+"/");
		
        int[] hourlyData = dset.loadHour(dataFolder+seed+"/", h);
        double[] otherQuakes = dset.getOtherQuakes(h);
        
        // Create EQ hourly information.
        EQHourlyInfo info = new EQHourlyInfo();
        
        info.parse(dset.numOfSites, dset.sampleRate, h, hourlyData, dset.EMA[h], otherQuakes); // EMA?
        
		return info;
	}
	
	// Load EQ relevant data.
	public QuakeTester.DataSet loadDataSet(String dataFolder, int seed) throws Exception {
		
		// Input arguments are assumed to be valid.
		
		// Load.
		QuakeTester qt = new QuakeTester();
		qt.seed = seed;
		
		QuakeTester.DataSet dset = qt.new DataSet(dataFolder+seed+"/");
		
		return dset;
	}
	
	// Practical training data unit.
	public class TrainingDataUnit {
		public double[] f = new double[NUM_FACTORS];
		public int tEQ;
	}
		
	// Earthquake prediction model.
	public class EQPredModel extends SparkNeuralNetwork {
		
		public EQPredModel(int numLayers, int[] numActs) {
			super(numLayers, numActs);
			// TODO Auto-generated constructor stub
		}
		
		// Earthquake event information.
		public EQEventInfo eInfo = new EQEventInfo();
		
		// Valid training information.
		public List<TrainingDataUnit> trainingInfo 
			= new ArrayList<TrainingDataUnit>();
			
		// Extract valid training information.
		public void extractValidTrainingInfo(EQEventInfo eInfo, Map<Integer, EQHourlyInfo> hInfos) {
			
			// Input arguments are assumed to be valid.
			
			// Make valid training information for the closest site from an earthquake.
			int site = eInfo.EQEventSite;
			int startHourIndex = PREDICTION_START_HOUR - 1;
			int stopHourIndex = eInfo.EQEventHour - 1;
			
			for (int hIndex = startHourIndex; hIndex < stopHourIndex; hIndex++) {
				
				// Create training data.
				TrainingDataUnit trainingData = new TrainingDataUnit();
				
				// Extract factors.

				// Make normalized magnetic field values of the closest site.
				// Make data via 3 channel averaging + sampling rate decreasing.
				double[][] sampleData = new double[NUM_CHANNELS][NUM_M_FACTOR];
				double[] aData = new double[NUM_M_FACTOR];
				
				int numAvgSamples = eInfo.sitesInfo.sampleRate * 3600 / NUM_M_FACTOR; // The assigned value is assumed to be an integer value.
				
				for (int i = 0; i < NUM_CHANNELS; i++) {
					sampleData[i] = new double[NUM_M_FACTOR];
					
					for (int j = 0; j < NUM_M_FACTOR; j++) {
						double sum = 0.0;
						
						for (int k = 0; k < numAvgSamples; k++) 
							sum += hInfos.get(hIndex).sitesData.get(site).data[i][j * numAvgSamples + k];
						
						sum /= numAvgSamples;
						sampleData[i][j] = sum;
					}
				}
			
				for (int i = 0; i < NUM_M_FACTOR; i++) {
					for(int j = 0; j < NUM_CHANNELS; j++) {
						aData[i] += sampleData[j][i];
					}
					
					aData[i] = aData[i] / NUM_CHANNELS / M_FACTOR_NORM_NUM;
				}
				
				// Noise reduction. Magnemeter noise?
				// Median filtering for high frequency filtering.
				aData = medianFiltering(aData);
				
				for (int i = 0; i < NUM_M_FACTOR; i++) 
					trainingData.f[i] = aData[i];
											
				// Earthquake event time.
				trainingData.tEQ = eInfo.EQEventHour - (hIndex + 1);
				
				trainingInfo.add(trainingData);
				printMsgSameLine("Tranining data " + trainingInfo.size() + " is added.");
			}
		}
		
		// Median filtering.
		private double[] medianFiltering(double[] raw) {
			double[] result = new double[raw.length];
			double[] samples = new double[3];
			
			samples[0] = raw[0];
			samples[1] = raw[0];
			samples[2] = raw[1];
			
			result[0] = calMedian(samples);
			
			for (int i = 1; i < raw.length - 1; i++) {
				samples[0] = raw[i - 1];
				samples[1] = raw[i];
				samples[2] = raw[i + 1];
				
				result[i] = calMedian(samples);
			}

			samples[0] = raw[raw.length - 2];
			samples[1] = raw[raw.length - 1];
			samples[2] = raw[raw.length - 1];
			
			result[raw.length - 1] = calMedian(samples);
			
			return result;
		}
		
		// Calculate a median value.
		private double calMedian(double[] data) {
			double[] raw = data.clone();
			
			// Sort.
			for (int i = 0; i < raw.length - 1; i++) {
				for (int j = i + 1; j < raw.length; j++) {
					if (raw[i] > raw[j]) {
						double temp = raw[i];
						raw[i] = raw[j];
						raw[j] = temp;
					}
				}
			}
			
			return raw.length % 2 == 0 
					? (raw[raw.length / 2 - 1] + raw[raw.length / 2]) / 2.0 : raw[raw.length / 2]; //?   
		}
		
		// Create the neural network parameter matrix.
		public void createNNParameterMatrix() {
			
			// It is assumed that the training data is obtained.
			
			// Create matrixes for both input and output features.
			int numSamples = trainingInfo.size();
			Matrix X = new Matrix(NUM_FACTORS, 1, trainingInfo.get(0).f);
			Matrix Y = new Matrix(PREDICTION_HOUR_TIME_LENGTH, 1, 0.0);
			Y.setVal(trainingInfo.get(0).tEQ, 1, 1.0);
			
			for (int i = 1; i < numSamples; i++) { // Memory limitation?
				X = X.horizontalAdd(new Matrix(NUM_FACTORS, 1, trainingInfo.get(i).f));
				
				Matrix tempY = new Matrix(PREDICTION_HOUR_TIME_LENGTH, 1, 0.0);
				tempY.setVal(trainingInfo.get(i).tEQ, 1, 1.0);
				
				Y = Y.horizontalAdd(tempY);
			}
			
			// Train.
			train(X, Y);
		}
		
		// Forecast.
		public double[] forecast(EQHourlyInfo hInfo) {
			
			// An input argument is assumed to be valid.
			
			// Calculate earthquake probability for each site.
			double[] result = new double[eInfo.sitesInfo.sites.size() * PREDICTION_HOUR_TIME_LENGTH];
			double[][] fs = new double[eInfo.sitesInfo.sites.size()][NUM_FACTORS];
			
			for (int siteIndex = 0; siteIndex < eInfo.sitesInfo.sites.size(); siteIndex++) {
				
				// Get a factor vector.
				fs[siteIndex] = new double[NUM_FACTORS];
				
				// Extract factors.
				fs[siteIndex][0] = 1.0;
				
				// Make normalized magnetic field values of the closest site.
				// Make data via 3 channel averaging + sampling rate decreasing.
				double[][] sampleData = new double[NUM_CHANNELS][NUM_M_FACTOR];
				double[] aData = new double[NUM_M_FACTOR];
				
				int numAvgSamples = eInfo.sitesInfo.sampleRate * 3600 / NUM_M_FACTOR; // The assigned value is assumed to be an integer value.
				
				for (int i = 0; i < NUM_CHANNELS; i++) {
					sampleData[i] = new double[NUM_M_FACTOR];
					
					for (int j = 0; j < NUM_M_FACTOR; j++) {
						double sum = 0.0;
						
						for (int k = 0; k < numAvgSamples; k++) 
							sum += hInfo.sitesData.get(siteIndex).data[i][j * numAvgSamples + k];
						
						sum /= numAvgSamples;
						sampleData[i][j] = sum;
					}
				}
			
				for (int i = 0; i < NUM_M_FACTOR; i++) {
					for(int j = 0; j < NUM_CHANNELS; j++) {
						aData[i] += sampleData[j][i];
					}
					
					aData[i] = aData[i] / NUM_CHANNELS / M_FACTOR_NORM_NUM;
				}
				
				// Noise reduction. Magnemeter noise?
				// Median filtering for high frequency filtering.
				aData = medianFiltering(aData);
				
				for (int i = 0; i < NUM_M_FACTOR; i++) 
					fs[siteIndex][i + 1] = aData[i];
			}
			
			// Calculate probability values of each hour.
			Matrix X = new Matrix(NUM_FACTORS, 1, fs[0]);
			
			for (int i = 1; i <= fs.length - 1; i++) {
				X = X.horizontalAdd(new Matrix(NUM_FACTORS, 1, fs[i]));
			}
			
			Matrix Y = judge(X);
			
			for (int siteIndex = 0; siteIndex < eInfo.sitesInfo.sites.size(); siteIndex++) { //?
				for (int i = 0; i < PREDICTION_HOUR_TIME_LENGTH; i++) {				
					result[PREDICTION_HOUR_TIME_LENGTH * eInfo.sitesInfo.sites.size() + i] 
							= Y.getVal(i + 1, siteIndex + 1);
				}
			}
			
			return result;
		}

		@Override
		public Matrix judge(Matrix X) {
			return feedForward(X);
		}
	}
		
	/** Earth quake prediction model. */
	public EQPredModel predModel;
	
	// Train. Caution: Memory limitation.
	public int train(String dataFolder, int[] seeds) throws Exception {
		
		// Input arguments are assumed to be valid.
		
		printMessage("Start training...");
		
		// Construct the earthquake prediction model.
		int[] acts = {360, 30, 60, 2160}; //?
		predModel = new EQPredModel(4, acts);
		
		// Load earthquake relevant information and make valid training information.	
		for (int seed : seeds) {
			printMessage("Load seed " + seed + " data.");
			EQEventInfo eInfo = loadEQEventInfo(dataFolder, seed);
			Map<Integer, EQHourlyInfo> hInfos = new HashMap<Integer, EQHourlyInfo>();
			
			for (int h = PREDICTION_START_HOUR - 1; h < eInfo.EQEventHour; h++) {
				printMsgSameLine("Load " + h + " data...");
				hInfos.put(h, loadEQHourlyInfo(dataFolder, seed, h));
			}
			
			// Extract valid training information.
			printMessage("Extract valid training information from seed " + seed + " data...");	
			predModel.extractValidTrainingInfo(eInfo, hInfos);
			hInfos.clear();
			hInfos = null;
		}
		
		// Create the neural network parameter matrix.
		printMessage("Create the neural network parameter matrix...");
		predModel.createNNParameterMatrix();
		predModel.trainingInfo.clear();
		
		return 0;
	}
			
	// Test.
	// Initialize.
	public int init(int sampleRate, int numOfSites, double[] sitesData) throws IOException {
		
		// Initialize the earthquake prediction model.
		predModel.trainingInfo.clear();
		
		// Get earthquake event information.
		predModel.eInfo.sitesInfo.parse(sampleRate, numOfSites, sitesData);
		
		return 0;
	}
	
	// Forecast earth quake.
	public double[] forecast(int hour, int[] data, double K, double[] globalQuakes) {
		
		// Input arguments are assumed to be valid.
				
		// Create earthquake hourly information.
		EQHourlyInfo hInfo = new EQHourlyInfo();
		
		hInfo.parse(predModel.eInfo.sitesInfo.sites.size()
				, predModel.eInfo.sitesInfo.sampleRate, hour, data, K, globalQuakes);
		
		// Forecast.
		double[] result = predModel.forecast(hInfo);
		
		return result;
	}
	
    public void printMessage(String s) {
        if (DEBUG) {
            System.out.println(s);
        }
    }

	// Print a message at a same line.
	public static void printMsgSameLine(String msg) {
		if (DEBUG) {
			System.out.print("\r" + msg);
		}
	}
}
