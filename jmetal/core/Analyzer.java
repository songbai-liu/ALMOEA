package jmetal.core;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

import jmetal.problems.Fonseca;
import jmetal.problems.Kursawe;
import jmetal.problems.ProblemFactory;
import jmetal.problems.Schaffer;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.ZDT1;
import jmetal.problems.ZDT.ZDT2;
import jmetal.problems.ZDT.ZDT3;
import jmetal.problems.ZDT.ZDT4;
import jmetal.problems.ZDT.ZDT6;
import jmetal.problems.cec2009Competition.UF1;
import jmetal.problems.cec2009Competition.UF2;
import jmetal.problems.cec2009Competition.UF3;
import jmetal.problems.cec2009Competition.UF4;
import jmetal.problems.cec2009Competition.UF5;
import jmetal.problems.cec2009Competition.UF6;
import jmetal.problems.cec2009Competition.UF7;
import jmetal.problems.cec2009Competition.UF8;
import jmetal.problems.cec2009Competition.UF9;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator1;
import jmetal.qualityIndicator.fastHypervolume.wfg.wfghvCalculator2;
import jmetal.util.Configuration;
import jmetal.util.JMException;

public class Analyzer {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;//java文件输出流，创建文件流
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;//OutputStreamWriter是字符流通向字节流的桥梁 
	      BufferedWriter bw      = new BufferedWriter(osw)        ;//缓冲区               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");//写到缓冲区
	        bw.newLine(); //换行       
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	
	
	public static void main(String[] args) throws JMException,
	SecurityException, IOException, ClassNotFoundException {
		jmetal.qualityIndicator.util.MetricsUtil utilities_;
		utilities_ = new jmetal.qualityIndicator.util.MetricsUtil();
		for(int fun=6;fun<=9;fun++){
			int runtimes=30;
			double[] IGDarray=new double[runtimes];
			double[] Hypervolume=new double[runtimes];
			
			QualityIndicator indicators; // Object to get quality indicators
			Problem problem=null; // The problem to solve
			
			indicators = null;
			if (args.length == 1) {
				Object[] params = { "Real" };
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else if (args.length == 2) {
				Object[] params = { "Real" };
				problem = (new ProblemFactory()).getProblem(args[0], params);
				indicators = new QualityIndicator(problem, args[1]);
			} // if
			else { // Default problem
				if(fun==1){
			  	      problem = new ZDT1("Real");
			  	      
			  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT1_501.txt" ) ;
			  	    	}//problem = new WFG1("Real");
			  	if(fun==2){
			  	      problem = new ZDT2("Real");
			  	      
			  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT2_501.txt" ) ;
			  	    	}//problem = new WFG1("Real");
			  	if(fun==3){
			  	      problem = new ZDT3("Real");
			  	      
			  	      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT3_269.txt" ) ;
			  	    	}
			  	if(fun==4){
				      problem = new ZDT4("Real",10);
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT4_501.txt" ) ;
				    	}//problem = new WFG1("Real");
				if(fun==5){
				      problem = new ZDT6("Real",10);
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\ZDT6_774.txt" ) ;
				    	}//problem = new WFG1("Real");
				if(fun==6){
					//problem = new DTLZ1("Real",10,2);
				      problem = new DTLZ1("Real",14,10);
				      
				      indicators = new QualityIndicator(problem,"E:\\truePF\\DTLZ1_10D.txt" ) ;
				      //indicators = new QualityIndicator(problem ) ;
				    	}
			  	if(fun==7){
				      //problem = new DTLZ2("Real",10,2);
			  		problem = new DTLZ2("Real",19,10);
				      
			  		indicators = new QualityIndicator(problem,"E:\\truePF\\DTLZ2_10D.txt" ) ;
			  		//indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\DTLZ2.500000p.10D.pf.txt" ) ;
				     //indicators = new QualityIndicator(problem);
				    	}//problem = new WFG1("Real");
				if(fun==8){
				      problem = new DTLZ3("Real",19,10);
				      
				      indicators = new QualityIndicator(problem,"E:\\truePF\\DTLZ3_10D.txt" ) ;
				      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\DTLZ3.500000p.10D.pf.txt" ) ;
				      //indicators = new QualityIndicator(problem);
				    	}//problem = new WFG1("Real");
				if(fun==9){
				      problem = new DTLZ4("Real",19,10);
				      
				      indicators = new QualityIndicator(problem,"E:\\truePF\\DTLZ4_10D.txt" ) ;
				      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\DTLZ4.500000p.10D.pf.txt" ) ;
				      //indicators = new QualityIndicator(problem);
				    	}
			  	if(fun==10){
				      problem = new DTLZ5("Real",19,10);
				      
				      indicators = new QualityIndicator(problem,"E:\\truePF\\DTLZ5_10D.txt" ) ;
				      //indicators = new QualityIndicator(problem);
				    	}//problem = new WFG1("Real");
				if(fun==11){
				      problem = new DTLZ6("Real",19,10);
				      
				      indicators = new QualityIndicator(problem,"E:\\truePF\\DTLZ6_10D.txt" ) ;
				      //indicators = new QualityIndicator(problem);
				    	}//problem = new WFG1("Real");
				if(fun==12){
				      problem = new DTLZ7("Real",27,8);
				      
				      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\DTLZ7.500000p.10D.pf.txt" ) ;
				      indicators = new QualityIndicator(problem);
				    	}
				if(fun==13){
			  	      problem = new WFG1("Real",10,20,6);
			  	      
			  	      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG2.500000p.5D.pf.txt" ) ;
			  	      indicators = new QualityIndicator(problem);
			  	    	}//problem = new WFG1("Real");
			  	if(fun==14){
			  	      problem = new WFG2("Real",10,20,6);
			  	      
			  	      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG2.500000p.5D.pf.txt" ) ;
			  	      indicators = new QualityIndicator(problem);
			  	    	}//problem = new WFG1("Real");
			  	if(fun==15){
			  	      problem = new WFG3("Real",10,20,6);
			  	      
			  	      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG3.500000p.5D.pf.txt" ) ;
			  	      indicators = new QualityIndicator(problem);
			  	    	}
				if(fun==16){
				      problem = new WFG4("Real",10,20,6);
				      
				      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" ) ;
				      indicators = new QualityIndicator(problem);
				    	}//problem = new WFG1("Real");
				if(fun==17){
				      problem = new WFG5("Real",10,20,6);
				      
				      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" ) ;
				      indicators = new QualityIndicator(problem);
				    	}
			  	if(fun==18){
			  	      problem = new WFG6("Real",10,20,6);
			  	      
			  	      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" ) ;
			  	      indicators = new QualityIndicator(problem);
			  	    	}//problem = new WFG1("Real");
			  	if(fun==19){
				      problem = new WFG7("Real",10,20,6);
				      
				      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" ) ;
				      indicators = new QualityIndicator(problem);
				    	}//problem = new WFG1("Real");
				if(fun==20){
				      problem = new WFG8("Real",10,20,6);
				      
				      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" );
				      indicators = new QualityIndicator(problem);
				    	}
				if(fun==21){
				      problem = new WFG9("Real",10,20,6);
				      
				      //indicators = new QualityIndicator(problem,"E:\\SRA\\SRA_2016\\truePF\\500000PF\\WFG4.500000p.5D.pf.txt" ) ;
				      indicators = new QualityIndicator(problem);
				    	}//problem = new WFG1("Real");
			  	if(fun==22){
				      problem = new Fonseca("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\Fonseca.pf" ) ;
				    	}//problem = new WFG1("Real");
				if(fun==23){
				      problem = new Kursawe("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\Kursawe.pf" ) ;
				    	}
				if(fun==24){
				      problem = new Schaffer("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\Schaffer.pf" ) ;
				    	}//problem = new WFG1("Real");
				if(fun==25){
				      problem = new UF1("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF1_500.txt" ) ;//.txt
				    	}
				if(fun==26){
				      problem = new UF2("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF2_500.txt" ) ;
				    	}
				if(fun==27){
				      problem = new UF3("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF3_500.txt" ) ;
				    	}
				if(fun==28){
				      problem = new UF4("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF4_500.txt" ) ;
				    	}
				if(fun==29){
				      problem = new UF5("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF5_21.txt" ) ;
				    	}
				if(fun==30){
				      problem = new UF6("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF6_668.txt" ) ;
				    	}
				if(fun==31){
				      problem = new UF7("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF7_500.txt" ) ;
				    	}
				if(fun==32){
				      problem = new UF8("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF8.DAT" ) ;
				    	}
				if(fun==33){
				      problem = new UF9("Real");
				      
				      indicators = new QualityIndicator(problem,"E:\\sbMaOP\\SbMaOP\\truePF\\UF9.DAT" ) ;
				    	}
			} // else
			
			for(int i=0;i<runtimes;i++){
				String path = "E:\\"+problem.getName()+"\\RVEA_"+problem.getName()+"_"+problem.getNumberOfObjectives()+"_run"+(i+1)+".txt";
				SolutionSet population = utilities_.readNonDominatedSolutionSet(path);
				IGDarray[i] = indicators.getIGD1(population);
				wfghvCalculator wfg = new wfghvCalculator(population,fun);
				Hypervolume[i] = wfg.calculatewfghv();
			}
			
			printGD("RVEA_10Obj_"+problem.getName()+"_HV.txt",Hypervolume);
			printGD("RVEA_10Obj_"+problem.getName()+"_IGD.txt",IGDarray);
		}
	}

}
