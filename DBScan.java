
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.TreeMap;

public class DBScan {

	private static final int UNCLASSIFIED = 0;
	private static final int NOISE = -1;
	static int clusterID = 1;
	static double[][] distanceMatrix;
	static Map<Integer,Integer> clusterMap=new TreeMap<Integer, Integer>();
	
	public static void main(String[] args) throws IOException {
		File file=new File("C:/Users/Rahul/OneDrive/Fall 2015/Data Mining/Project2/iyer.txt");
		HashMap<Integer, List<Double>> geneMap = new HashMap<Integer, List<Double>>();
		HashMap<Integer, Integer> groundTruth = new HashMap<Integer, Integer>();
		
		//INPUT PARAMETERS:
		int minpts = 3;
		Double eps = 2.78;
		
		Scanner input = new Scanner(file);
		
		while(input.hasNext())
		{
			String inputRecord = input.nextLine();
			String[] geneRecord=inputRecord.split("\t");
			
			List<Double> recordPointList = new ArrayList<Double>();
			
			int geneId = Integer.parseInt(geneRecord[0]);
			int geneTruth = Integer.parseInt(geneRecord[1]);
			
			for(int i=2; i<geneRecord.length; i++)
			{
				recordPointList.add(new Double(geneRecord[i]));
			}
			
			geneMap.put(geneId, recordPointList);
			groundTruth.put(geneId, geneTruth);
			clusterMap.put(geneId,UNCLASSIFIED);
		}
		input.close();

		distanceMatrix = new double[geneMap.size()][geneMap.size()];
	
		for(Map.Entry<Integer,Integer> entry:clusterMap.entrySet())
		{
			int geneIndex=entry.getKey();
			if(clusterMap.get(geneIndex) == UNCLASSIFIED)
			{
				assignCluster(geneMap, geneIndex, eps, minpts);
			}
			
		}
		
		System.out.println("\nAfter Clustering:");
		
		for(Map.Entry<Integer,Integer> entry:clusterMap.entrySet())
		{
			System.out.println("GeneID: " + entry.getKey() + "\tClusterID: " + entry.getValue());
		}
		
		// VALIDATION:
		System.out.println("\nValidation:\n");
		
		//1. External Index using Jaccard Coefficient
		ValidateCoefficients validate = new ValidateCoefficients();
		double externalIndex = validate.JaccardCoefficient(clusterMap, groundTruth);
		System.out.println("Jaccard Coefficient: " + externalIndex + "\n");
		
		//2. Internal Index using Silhouette Coefficient
		validate.silhouetteCoefficient(clusterMap, distanceMatrix,clusterID);
		
		FileWriter fstream = new FileWriter("DBScanResults.txt");
        BufferedWriter out = new BufferedWriter(fstream);
        
        for(Map.Entry<Integer,Integer> entry:clusterMap.entrySet())
		{
        	out.write(entry.getValue() + "\n");
		}
        
        out.close();
	}
	
	public static void assignCluster(HashMap<Integer, List<Double>> geneMap, int geneIndex, Double eps, int minPts)
	{
		List<Integer> neighbourList = neighbours(geneMap, geneIndex, eps);
		
		if(neighbourList.size() >= minPts-1)
		{
			clusterMap.put(geneIndex,clusterID);
			for(int neighbourIndex = 0; neighbourIndex<neighbourList.size(); neighbourIndex++)
			{
				clusterMap.put(neighbourList.get(neighbourIndex),clusterID);
				extendForCorePoints(geneMap, eps, neighbourList.get(neighbourIndex), minPts);
				//assignCluster(geneMap, neighbourList.get(neighbourIndex), eps, minPts);
			}
			clusterID++;
			
		}
		else //if(clusterMap.get(geneIndex)==UNCLASSIFIED)
		{
			clusterMap.put(geneIndex,NOISE);
		}
	}

	private static void extendForCorePoints(HashMap<Integer, List<Double>> geneMap, Double eps, Integer corePointIndex,int minPts) 
	{
		List<Integer> associatedNeighbours = neighbours(geneMap, corePointIndex, eps);
		if(associatedNeighbours.size() >= minPts-1)
		{
			for(int index = 0; index<associatedNeighbours.size(); index++)
			{
				if(clusterMap.get(associatedNeighbours.get(index)) == UNCLASSIFIED || clusterMap.get(associatedNeighbours.get(index)) == NOISE)
				{
					clusterMap.put(associatedNeighbours.get(index),clusterID);				
				}
			}
		}
	}
	
	public static List<Integer> neighbours(HashMap<Integer, List<Double>> geneMap, int geneIndex, Double eps)
	{
		List<Integer> neighbourList = new ArrayList<Integer>();
		for(Entry<Integer, List<Double>> entry : geneMap.entrySet()) {
		{	
			int compareGeneIndex=entry.getKey();
			
			if(compareGeneIndex==geneIndex)
			{
				distanceMatrix[compareGeneIndex-1][geneIndex-1] = 0.0;
				continue;
			}
			
			Double distance = euclideanDistance(geneMap.get(compareGeneIndex), geneMap.get(geneIndex));
			distanceMatrix[compareGeneIndex-1][geneIndex-1] = distance;
			distanceMatrix[geneIndex-1][compareGeneIndex-1] = distance;
			
			if(distance<eps)
				neighbourList.add(compareGeneIndex);
			}
		}
		return neighbourList;
	}
	
	public static Double euclideanDistance(List<Double> p, List<Double> q)
	{
		int listSize = p.size();
		Double distanceSum = 0.0;
		for(int listIndex = 0; listIndex<listSize; listIndex++)
		{
			Double difference = 0.0;
			difference = q.get(listIndex) - p.get(listIndex);
			distanceSum += Math.pow(difference,2);
		}
		return Math.sqrt(distanceSum);
	}
	
	
}
