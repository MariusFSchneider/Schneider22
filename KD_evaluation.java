// this plugin implements Particle Analyzer, the watersheding and tresholding tool, similar to Nucleus counter
// By stating threshold values the plugin compares whether a pixels is green and red, green and blue, red and blue or all and generates overlay pictures.
// plugin works on series 
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;
import java.awt.List;
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.lang.*;

import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.filter.Analyzer;
import ij.plugin.frame.PlugInFrame;
import ij.*;
import ij.process.*;
import ij.io.FileSaver;
import java.awt.image.*;
import java.awt.*;
import java.io.*;
import ij.io.OpenDialog;
import java.util.Iterator;
import javax.imageio.*;
import javax.imageio.stream.*;

public class KD_evaluation extends java.lang.Object implements PlugIn {

// imp1 = green channel
// imp2 = red channel
// imp3 = overlay: green/red = transduced neurons
// imp4 = blue channel
// imp5 =  overlay: blue and green = nuclei of transduced cells. Cells are counted by particle analyzer 
// imp6 = overlay: red and blue = nuclei of neurons. Cells are counted by particle analyzer 
// imp7 = overlay: red/green and blue = nuclei of transduced neurons. Cells are counted by particle analyzer. 
// imp8 = modified blue channel



    private ImagePlus imp1;
	
	// imp1 = DAPI channels
    private ImagePlus imp2;
	private ImagePlus imp3;
	private ImagePlus imp4;
	private ImagePlus imp5;
	private ImagePlus imp6;
	private ImagePlus imp7;
	public  ImagePlus imp8;
	public float[] fl;
	public double [] areas;
	public double meanA = 0;
	public double STDa = 0;
	public String Dapi;
	public String GFP;
	public String TujI;
	public String directory = (String)Prefs.get("", "C:/");
	public ResultsTable rt; 
	public ResultsTable summaryTable = new ResultsTable(); 
	
	public String path = "";
	
	public int maxG = (int)Prefs.get("NC_max.intA", 500);
	public int maxR = (int)Prefs.get("NC_max.intB", 500);
	public int maxB = (int)Prefs.get("NC_max.intC", 900);
	// area ints count green, red and both colored pixels, respectively , but are not used further
	public int area_green;
	public int area_red;
	public int area_overlay; 
	public String Titel;
	public String Titel1;
	public String Titel2;
	public String Titel3;
	public String Titel4;
// nucleus counter 	
// int max = largest particle size
// int minP = smalles particle size
	public int max = (int)Prefs.get("NC_max.int", 5000);
	public int minP = (int)Prefs.get("NC_min.intP", 150);

	
//	public int threshIndex=(int)Prefs.get("NC_threshIndex.int",0);
	public int smoothIndex=(int)Prefs.get("NC_smoothIndex.int",0);
	public boolean watershed=Prefs.get("NC_watershed.boolean",true);
	public boolean summarize=Prefs.get("NC_summarize.boolean",true);
	public boolean add=Prefs.get("NC_record.boolean",true);
	public boolean fullStats=Prefs.get("NC_fullStats.boolean",true);
// int minthreshold2 = 12bit value	
	public int minthreshold2;
// int minthreshold = not used anymore
	public int minthreshold;
// int minThreshold = 8bit value	
	public int minThreshold;
// int maxThreshold = 8bit value
	public int maxThreshold = 255;
// int maxhreshold = not used anymore
	public int maxthreshold;

	public double ln1;
	public double ln2;
	public double ln3;
	public double ratio;
	public int choice;	
public boolean overlay(ImagePlus imp1, ImagePlus imp2, ImagePlus imp4, int maxG, int maxR, int maxB, String Titel, String Titel2, String Titel3, String Titel4) {
	

	summaryTable.incrementCounter();
	int area_green = 0;
    int area_red = 0;
	
	int w = imp1.getWidth();
	int h = imp1.getHeight();
	
	imp3 = NewImage.createShortImage(Titel, w, h,
				1, NewImage.FILL_BLACK);
	imp5 = NewImage.createShortImage(Titel2, w, h,
				1, NewImage.FILL_BLACK);
	imp6 = NewImage.createShortImage(Titel3, w, h,
				1, NewImage.FILL_BLACK);	

	imp7 = NewImage.createShortImage(Titel4, w, h,
				1, NewImage.FILL_BLACK);	
				
	imp8 = NewImage.createShortImage("DAPI modified", w, h,
				1, NewImage.FILL_BLACK);	
				
	ImageProcessor IP1 = imp1.getProcessor();
	ImageProcessor IP2 = imp2.getProcessor();
	ImageProcessor IP3 = imp3.getProcessor();
	
	ImageProcessor IP4 = imp4.getProcessor();
	ImageProcessor IP5 = imp5.getProcessor();
	ImageProcessor IP6 = imp6.getProcessor();	
	ImageProcessor IP7 = imp7.getProcessor();
	ImageProcessor IP8 = imp8.getProcessor();
	
	short[] pixels1 = (short[])IP1.getPixels();
	short[] pixels2 = (short[])IP2.getPixels();
	short[] pixels3 = (short[])IP3.getPixels();

	short[] pixels4 = (short[])IP4.getPixels();
	short[] pixels5 = (short[])IP5.getPixels();
	short[] pixels6 = (short[])IP6.getPixels();	
	short[] pixels7 = (short[])IP7.getPixels();	
	short[] pixels8 = (short[])IP8.getPixels();	

// calculation of overlay images
	
	for (int i =0; i < pixels1.length; i++){
		
			int p1 = pixels1[i];
			int p2 = pixels2[i];
			int p4 = pixels4[i];
			int p3 = 0;
			int p5 = 0;
			int p6 = 0;
			int p7 = 0;
			int p8 = 0;
			
			
			
			if (p1 >= maxG){
				area_green = area_green +1;
			}
			if (p2 >= maxR){
				area_red = area_red + 1;
			}			
			if (p1 >= maxG && p2 >= maxR) {
				p3 = (p1+p2)/2;
			}
			if (p1 >= maxG && p2 >= maxR && p4 >= maxB) {
	
				p7 = p4;
			}
			if (p1 >= maxG && p4 >= maxB) {
				p5 = p4;
			}	
			if (p2 >= maxR && p4 >= maxB) {
				p6 = p4;
			}		
			if (p4 >= maxB){
				p8 = 4095;
			}
			pixels3[i] = (short)(p3);
			pixels5[i] = (short)(p5);
			pixels6[i] = (short)(p6);
			pixels7[i] = (short)(p7);	
			pixels8[i] = (short)(p8);
		}	

	
// saving overlay pictures as jpg files
	imp3.show();
	imp3.updateAndDraw();
	path = directory + "/" + Titel4 + "2" + ".jpg";
	path = path.replace("," , "-" );
	saveAsJpeg(imp3, path, 75);
	this.imp3.close();
	
	imp5.show();
	imp5.updateAndDraw();	
	path = directory+ "/"  + Titel2 + ".jpg";
	path = path.replace("," , "-" );	
	saveAsJpeg(imp5, path, 75);

	
	imp6.show();
	imp6.updateAndDraw();	
	path = directory+ "/"  + Titel3+ ".jpg";
	path = path.replace("," , "-" );	
	saveAsJpeg(imp6, path, 75);

	
	
	imp7.show();
	imp7.updateAndDraw();
	path = directory + "/" + Titel4+ ".jpg";
	path = path.replace("," , "-" );
	saveAsJpeg(imp7, path, 75);


	imp8.show();
	imp8.updateAndDraw();
	path = directory + "/" + Titel1 + ".jpg";
	path = path.replace("," , "-" );
	saveAsJpeg(imp8, path, 75);
	
// running of Particle Analyzer method for each overlay: 	

// running of Particle Analyzer method for each overlay: 	



// total cells, DAPI modified
	int[] minU1 = new int [1];
	minU1[0] = maxB;
	count(minU1, this.imp8, Titel1, true);
	this.imp8.close();
	
// transduced cells
	int[] minU2 = new int [2];
	minU2[0] = maxB;	
	minU2[1] = maxG;	
	count(minU2, this.imp5, Titel2, true);
	this.imp5.close();
	
// neurons	
	int[] minU4 = new int [2];
	minU4[0] = maxB;	
	minU4[1] = maxR;	
	count(minU4,  this.imp6, Titel3, true);
	this.imp6.close();
	
// KD neurons
	int[] minU5 = new int [3];
	minU5[0] = maxB;	
	minU5[1] = maxG;
	minU5[2] = maxR;
	count(minU5, this.imp7, Titel4, false);
	this.imp7.close();
	



	
	summaryTable.show("summary KD evaluation");	
//	IJ.showMessage("area of green pixels = " + area_green + "\n" + "area of red pixels = " + area_red);
	return true;
	}
	

public void ParticlAnalyzer(ImagePlus imp, String Titel, int minThreshold, int minthreshold2){
// many parts copied from nucleus counter, default threshold values were changed, respectively.... Default Thresholding method = "OtsuThresholding 16Bit"
	
	ImageProcessor ip1 = imp.getProcessor();
	ImageProcessor ip2 = ip1.duplicate();
	new ImagePlus("Analysis", ip2).show();
	IJ.run("Grays");

	ImagePlus imp2 = WindowManager.getCurrentImage();
	ImageWindow winimp2=imp2.getWindow();

				

//duplicate this for image thresholding
	ImageProcessor ip3 = ip2.duplicate();
	ip3.setThreshold(minThreshold,maxThreshold,4);
	new ImagePlus("Threshold", ip3).show();
	ImagePlus imp3 = WindowManager.getCurrentImage();
	ImageWindow winimp3=imp3.getWindow();
	ImagePlus imp5=null;
	ImageWindow winimp5=null;
	//IJ.showMessage("thres="+threshIndex);
		
//set threshold

		if (smoothIndex==1) IJ.run("Mean...", "radius=2 separable");
		else if(smoothIndex==2) IJ.run("Mean...", "radius=3 separable");
		else if(smoothIndex==3) IJ.run("Median...", "radius=2");
		else if(smoothIndex==4) IJ.run("Median...", "radius=3");
	


		WindowManager.setCurrentWindow(winimp3);
		
//		if (threshIndex==1&&imp1.getType()==imp1.GRAY8)	IJ.run("OtsuThresholding 8Bit");
//		if (threshIndex==1&&imp1.getType()==imp1.GRAY16) 
		IJ.run("OtsuThresholding 16Bit");
//		if (threshIndex==2){IJ.run("8-bit");	IJ.run("Entropy Threshold");}
//		if (threshIndex==3){IJ.run("8-bit");	IJ.run("Mixture Modeling threshold");}

//		minthreshold =(int)ip3.getMinThreshold();
//		maxthreshold =(int)ip3.getMaxThreshold();	
//		WindowManager.setCurrentWindow(winimp3);
//		IJ.setThreshold(minThreshold, maxThreshold);

//		if (threshIndex==5) 
//			{WindowManager.setCurrentWindow(winimp3);
//			IJ.run("Adapative3DThreshold ");
//			imp5 = WindowManager.getCurrentImage();
//			ImageProcessor ip5=imp5.getProcessor();
//			winimp5=imp5.getWindow();
//			WindowManager.setCurrentWindow(winimp5);
//			IJ.run("Rename...", "title=Threshold");
//			IJ.setThreshold(1, 255); imp5.changes=false;
//			}

//		if (threshIndex==4) 
//			{WindowManager.setCurrentWindow(winimp3);
//			IJ.run("k-means Clustering", "number=2 cluster=0.00010000 randomization=48");
//			imp5 = WindowManager.getCurrentImage();
//			ImageProcessor ip5=imp5.getProcessor();
//			winimp5=imp5.getWindow();
//			WindowManager.setCurrentWindow(winimp5);
//			IJ.run("Rename...", "title=Threshold");
//			IJ.run("Invert");
//			IJ.setThreshold(2, 2); imp5.changes=false;
//			}


		
		//IJ.showMessage("Min="+minthreshold+"   Max = "+maxthreshold);

		IJ.run("Convert to Mask");
//watershed
		if (watershed) IJ.run("Watershed");
		IJ.setThreshold(minThreshold,maxThreshold);

		ImagePlus mask = WindowManager.getCurrentImage();


//IJ.run("Analyze Particles...", "size=0.65-647.67 circularity=0.00-1.00 show=Outlines display exclude clear summarize");

//analyze particles
		String analyseStr = "size=" + minP +"-" + max +" circularity=0.5-1.0 show=";
		analyseStr+="Outlines";
		//analyseStr+=" display";
		analyseStr+=" exclude clear";

		if (summarize) analyseStr+=" summarize";
		if (add) analyseStr+=" add";
		//if (fullStats) analyseStr+=" size";

//IJ.showMessage(analyseStr);

		IJ.run("Analyze Particles...",   analyseStr+ " add");


//get Outline image
		ImagePlus imp4 = WindowManager.getCurrentImage();
		ImageProcessor ip4 = imp4.getProcessor();
		ImageWindow winimp4 = imp4.getWindow();
		IJ.run("Rename...", "title=boundaries");
		path = directory + "/" + Titel+ "-outlines.jpg";
		path = path.replace("," , "-" );
		saveAsJpeg(imp4, path, 75);
		
		
	//	IJ.run("Invert");		
	//	WindowManager.setCurrentWindow(winimp2);
	//
	//	IJ.run("8-bit");
	//	IJ.run("Image Calculator...", "image1='Analysis' operation=Subtract image2='boundaries' ");
	//	WindowManager.setCurrentWindow(winimp4);
	//	IJ.run("Red");
	//	IJ.run("RGB Color");

//get duplicate image		
	//	WindowManager.setCurrentWindow(winimp2);
	//	IJ.run("RGB Color");
	//	IJ.run("Image Calculator...", "image1='Analysis' operation=Add image2='boundaries' create ");

		imp2.changes = false; winimp2.close();
//imp3.changes = false;
	imp4.changes = false;
		
	winimp4.close();	

//		if (!add||threshIndex==4) {winimp3.close();}
//		if (threshIndex==4||threshIndex==5&&!add) {imp5.changes=false; winimp5.close();}
		
//		if(add)
//			{
//			if(threshIndex==4) 
///				{WindowManager.setCurrentWindow(winimp5);
//				IJ.setThreshold(128,255);}
//			else 
//				{WindowManager.setCurrentWindow(winimp3);
//				IJ.setThreshold(minthreshold, maxthreshold);}
			//IJ.run("Multi Measure");
			//IJ.showMessage("Click on MultiMeasure 'Add Particles' button\nnext after clearing existing list if required. \nThen close 'Threshold' image.");
					
//			}
		mask.changes=false;

		mask.close();	
}
public void run(String arg) {

	this.showDialog();


		
	
	}
// Generic Dialog Box	
public boolean showDialog() {
	/** 
 * Generate Dialog box
 * Get set image processor
 * overlay:
 */

	String[] channel = new String[4];
	channel[0] = "C=0";
	channel[1] = "C=1";
	channel[2] = "C=2";
	channel[3] = "C=3";
	
	int index1 = 2;
	int index2 = 3;
	int index3 = 1;
	
// Nucleus counter
//	String [] 	thresholds  = {"Current", "Otsu", "Max. Entropy", "Mixture Modelling", "k-means Clustering", "Adaptive"};	
	String [] 	smooths= {"None", "Mean 3x3", "Mean 5x5", "Median 3x3", "Median 5x5"};	
//  overlay	
 	GenericDialog gd = new GenericDialog("KD Evaluation plKO.1-GFP vectors");
	gd.addStringField("Directory of output files: ", directory);
	gd.addNumericField("Threshold GFP Intensity", maxG,0);
	gd.addNumericField("Threshold TujI staining", maxR,0);
	gd.addNumericField("Threshold DAPI", maxB,0);
    gd.addChoice("GFP channel", channel, channel[index1]);
    gd.addChoice("TujI staining", channel, channel[index2]);
	gd.addChoice("DAPI staining", channel, channel[index3]);
	
	
//  Nucleus counter
	gd.addNumericField("Smallest Particle Size", minP,0);
	gd.addNumericField("Largest Particle Size", max,0);
//	gd.addChoice("Threshold Method", thresholds,thresholds[threshIndex]);
	gd.addChoice("Smooth Method", smooths,smooths[smoothIndex]);

//	gd.addCheckbox("Smooth prior to segmentation",smooth);
	gd.addCheckbox("Watershed filter",watershed);
	gd.addCheckbox("Add Particles to ROI manager",add);
//	gd.addCheckbox("Show full statistics",fullStats);
	gd.addCheckbox("Show summary",summarize);
	
	
	gd.showDialog();        
	if (gd.wasCanceled()) {
	         return false;
	}
	
// overlay
    index1 = gd.getNextChoiceIndex();
    index2 = gd.getNextChoiceIndex();
	index3 = gd.getNextChoiceIndex();
	directory = gd.getNextString();
	GFP = channel[index1];
    TujI = channel[index2];
	Dapi = channel[index3];
	maxG=(int)gd.getNextNumber();
	maxR=(int)gd.getNextNumber();
	maxB=(int)gd.getNextNumber();
	
// Nucleus counter
	minP=(int)gd.getNextNumber();
	max=(int)gd.getNextNumber();
//	threshIndex = 1;
	smoothIndex=gd.getNextChoiceIndex();


	watershed=gd.getNextBoolean();
	add=gd.getNextBoolean();
//	fullStats=gd.getNextBoolean();
	summarize=gd.getNextBoolean();
   ImagePlus imagePlus;	
   int[] arrn = WindowManager.getIDList();
	String[] arrstring = new String[arrn.length];
    int n = 0;
    while (n < arrn.length) {
		
        imagePlus = WindowManager.getImage((int)arrn[n]);
        if (imagePlus != null) {
             arrstring[n] = imagePlus.getTitle();
	
        } else {
        arrstring[n] = "";
        }
        ++n;
        }
	String imA = "";
	String imB = "";
	String imC = "";
	String imD = "";
	
// looping through opened lif file	
		
	for (int i= 0; i < (arrstring.length)-3; i++){
//	IJ.showMessage("click");
		if (arrstring[i].contains("=") && arrstring[i+1].contains("=") && arrstring[i+2].contains("=") && arrstring[i+3].contains("=")){
//	IJ.showMessage("clickagain");		
			imA = arrstring[i];
			imB = arrstring[i+1];
			imC = arrstring[i+2];
			imD = arrstring[i+3];
			String[] ImA = imA.split("-\\s+");
			String[] ImB = imB.split("-\\s+");
			String[] ImC = imC.split("-\\s+");
			String[] ImD = imD.split("-\\s+");	
				if(ImA[1].equals(ImB[1]) &&  ImC[1].equals(ImD[1]) && ImD[1].equals(ImA[1])){
	
//				IJ.showMessage("toomanyClicks");
				this.imp1 = WindowManager.getImage((int)arrn[i+ index1]);
				IJ.run("Subtract Background...", "rolling=50");
				this.imp2 = WindowManager.getImage((int)arrn[i+ index2]);
				IJ.run("Subtract Background...", "rolling=50");
				this.imp4 = WindowManager.getImage((int)arrn[i+ index3]);
				IJ.run("Subtract Background...", "rolling=50");				
				Titel1 = ImA[1] + "," + "total cells";
				Titel2 = ImA[1] + "," + "transduced cells";
				Titel3 = ImA[1] + "," + "neurons";
				Titel4 = ImA[1] + "," + "transduced neurons";
				this.overlay(this.imp1, this.imp2, this.imp4, maxG, maxR,maxB, Titel, Titel2, Titel3, Titel4);

				
				}

		}
	}
	return true;

}
String saveAsJpeg(ImagePlus imp, String path, int quality) {
        int width = imp.getWidth();
        int height = imp.getHeight();
        int biType = BufferedImage.TYPE_INT_RGB;
        boolean overlay = imp.getOverlay()!=null && !imp.getHideOverlay();
        if (imp.getProcessor().isDefaultLut() && !imp.isComposite() && !overlay)
            biType = BufferedImage.TYPE_BYTE_GRAY;
        BufferedImage bi = new BufferedImage(width, height, biType);
        String error = null;
        try {
            Graphics g = bi.createGraphics();
            Image img = imp.getImage();
            if (overlay)
                img = imp.flatten().getImage();
            g.drawImage(img, 0, 0, null);
            g.dispose();            
            Iterator iter = ImageIO.getImageWritersByFormatName("jpeg");
            ImageWriter writer = (ImageWriter)iter.next();
            File f = new File(path);
            String originalPath = null;
            boolean replacing = f.exists();
            if (replacing) {
                originalPath = path;
                path += ".temp";
                f = new File(path);
            }
            ImageOutputStream ios = ImageIO.createImageOutputStream(f);
            writer.setOutput(ios);
            ImageWriteParam param = writer.getDefaultWriteParam();
            param.setCompressionMode(param.MODE_EXPLICIT);
            param.setCompressionQuality(quality/100f);
            if (quality == 100)
                param.setSourceSubsampling(1, 1, 0, 0);
            IIOImage iioImage = new IIOImage(bi, null, null);
            writer.write(null, iioImage, param);
            ios.close();
            writer.dispose();
            if (replacing) {
                File f2 = new File(originalPath);
                boolean ok = f2.delete();
                if (ok) f.renameTo(f2);
            }
        } catch (Exception e) {
            error = ""+e;
            IJ.error("Jpeg Writer", ""+error);
        }
        return error;
    }
public void count (int[] min, ImagePlus impA, String Titel, boolean incrementCol){
	

		try {
		minthreshold = getMinValue(min);
		ln1 = Math.log(choice);
		ratio = ln1 / 12;
		ln2 = 8*ln2;
		ln3 = Math.exp(ln2);
		minThreshold = (int)ln3;
		ParticlAnalyzer(impA, Titel, minThreshold, minthreshold);
		rt = Analyzer.getResultsTable();
		fl = rt.getColumn(0);
		int l4 = fl.length;
		areas = rt.getColumnAsDoubles(0);
		meanA = mean_area(areas);
		STDa = std_area(meanA, areas);
		summaryTable.addValue("sample", Titel);
		summaryTable.addValue("counts", l4);	
		summaryTable.addValue("mean_area", meanA);
		summaryTable.addValue("std_area", STDa);
		if (incrementCol == true){
			summaryTable.incrementCounter();
			}
		
		impA.close();
	}   catch (Exception e) {  
      //      e.printStackTrace();  
		summaryTable.addValue("sample", Titel);	
		summaryTable.addValue("counts", 0);	
		summaryTable.addValue("mean_area", 0);
		summaryTable.addValue("std_area", 0);		
		if (incrementCol == true){
			summaryTable.incrementCounter();
			}
		impA.close();
 //     return false;  
        } 

}	
public double mean_area (double [] areas){
	int length = areas.length;
	double sum = 0;
	for (int i =0; i < length; i ++){
		sum = sum + areas[i];
	}
	return sum/length;
	
}	

public double std_area (double mean, double[] areas){
	
	int length = areas.length;
	double sum = 0;
	for (int i =0; i < length; i ++){
		sum = sum + ((areas[i] - mean)*(areas[i] - mean));
	}
	double var = sum/(length-1);
	double std = Math.sqrt(var);
	return std;
	
}
public int getMinValue(int[] min){
	int Min = 0;
	Min = min[0];
	for (int i = 0; i < min.length; i++){
		if (min[i] <= Min){
			Min = min[i];
		}
	}
	return Min;
}
}

