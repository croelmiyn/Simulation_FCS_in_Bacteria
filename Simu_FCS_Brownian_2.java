import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;
import mpi.rc.IJ.IJutilities.MersenneTwister;

import java.io.FileWriter;

public class Simu_FCS_Brownian_2 implements PlugIn {

    /* Simulation of FCS in bacteria
     *  Boundary conditions: Redrawing the step if it falls out of the bacterium.
     *  Code that was used for the published results
     */

    private double L,d; // length diameter of the cell
    private double w0, z0, w0r, z0r; // width Gaussian beam
    private double xc,yc,zc; // coordinates of the center of the illumination beam

    private double D, Dr; // diffusion coefficient, sqrt(2D)

    private double dt; // time step
    private double Tmax; // duration

    private int nFrames;
    private int nObjects;

    private double[][] pos;

    private double a,b; // half length-half diameter / half diameter

    private double[] I, CI; //FCS intensity + correlation function
    private double Im;

    private DoubleFFT_1D fFT1D;

    private MersenneTwister rd;

    private String dir,saveNameCorr,saveNameParam,lineSep;

    public void run(String arg){

        initialize();

        simulate();

        compute();

        save();

    }

    private void initialize(){

        // insert dialog here

        lineSep = System.getProperty("line.separator");

        Tmax = 10;
        dt = 1e-6;

        // particle
        nObjects = 100;
        D = 10;

        // cell
        L = 5;
        d = 1;

        // laser beam
        w0 = 0.2;
        z0 = 0.6;

        xc = 0.5;
        yc = 0;
        zc = 0;

        String parm = "test";

        GenericDialog gd = new GenericDialog("nom du fichier de sortie");
        gd.addNumericField("Simu_duration (s) ", Tmax, 1);
        gd.addNumericField("Time_step (s) ", dt, 6);
        gd.addNumericField("Number_of_Objects ", nObjects, 0);
        gd.addNumericField("Diffusion_coeff (um2/s) ", D, 1);
        gd.addMessage("Cell Parameters");
        gd.addNumericField("Cell_length (um)", L, 1);
        gd.addNumericField("Cell_diameter (um)", d, 1);
        gd.addMessage("Laser beam Parameters");
        gd.addNumericField("w0 (um) ", w0, 2);
        gd.addNumericField("z0 (um) ", z0, 2);
        gd.addMessage("Laser beam position relative to cell tip");
        gd.addNumericField("xc (um) ", xc, 1);
        gd.addNumericField("yc (um) ", yc, 1);
        gd.addNumericField("yc (um) ", zc, 1);

        gd.addMessage("Save option");
        gd.addStringField("label:", parm, 0);

        gd.showDialog();
        Tmax = gd.getNextNumber();
        dt = gd.getNextNumber();
        nObjects = (int) gd.getNextNumber();
        D = gd.getNextNumber();
        L = gd.getNextNumber();
        d = gd.getNextNumber();

        w0 = gd.getNextNumber();
        z0 = gd.getNextNumber();
        xc = gd.getNextNumber();
        yc = gd.getNextNumber();
        zc = gd.getNextNumber();

        parm = gd.getNextString();

        // save dialogs

        saveNameCorr = "CorrFct_SimuFCS_"+parm+"";
        saveNameParam = "Params_SimuFCS_"+parm+"";


        SaveDialog sd = new SaveDialog("Output_File",saveNameCorr,".txt");
        saveNameCorr = sd.getDirectory()+sd.getFileName();

        sd = new SaveDialog("Log_File","",saveNameParam,".txt");
        saveNameParam = sd.getDirectory() + sd.getFileName();


        // initialize variables

        rd = new MersenneTwister();

        a = L/2-d/2;
        b = d/2;

        Dr = Math.sqrt(2*D*dt);
        w0r = w0*w0/2;
        z0r = z0*z0/2;

        pos = new double[nObjects][3];

        nFrames = (int)(Tmax/dt);
        I = new double[nFrames];

        fFT1D = new DoubleFFT_1D(2*nFrames); // zero padding    
    }

    private void simulate(){

        double x,y,z,X,r,F,dx,dy,dz;
        Im=0;
        boolean notFoundLegalStep;

        // initialize positions
        for(int k=0;k<nObjects;k++){
            int cnt =0;
            while(!correctStartPosition(k)&&cnt<(1e4*nObjects)){cnt++;}
        }

        // run simulation

        for(int t=0; t<nFrames; t++){

            for(int k=0;k<nObjects;k++){
                x = 0.0;
                y = 0.0;
                z = 0.0;
                notFoundLegalStep = true;

                while(notFoundLegalStep) {
                    dx = Dr * rd.nextGaussian();
                    dy = Dr * rd.nextGaussian();
                    dz = Dr * rd.nextGaussian();

                    // boundary check
                    x = pos[k][0] + dx;
                    y = pos[k][1] + dy;
                    z = pos[k][2] + dz;

                    X = Math.abs(x) - a;
                    X = (X > 0) ? X : 0;

                    r = Math.sqrt(X * X + y * y + z * z);

                    notFoundLegalStep=(r > b); // hard wall boundary condition implementation
                }

                pos[k][0] = x;
                pos[k][1] = y;
                pos[k][2] = z;


                // intensity computation

                I[t] += intensityFCS(pos[k]);

            }
            Im+=I[t];

        }
        Im/=nFrames;

    }

    private void compute(){

        // compute correlation function

        double[] ic;

        ic = toComplex(I);

        fFT1D.complexForward(ic);

        double[] psd = new double[4*nFrames];

        for(int k=0; k<2*nFrames; k++){
            psd[2*k] = ic[2*k]*ic[2*k] + ic[2*k+1]*ic[2*k+1];
        }

        fFT1D.complexInverse(psd,true);

        CI = new double[nFrames];

        for(int k=0;k<nFrames;k++){
            CI[k] = psd[2*k]/(double)(nFrames-k) - Im*Im;
        }

    }

    private void save(){
        saveCorrFct();
        saveParams();
    }

    private void saveCorrFct(){
        String buffer = "tau\t<I(t)I(0)>-<I>^2\t(<I(t)I(0)>-<I>^2)/(<I^2>-<I>^2)\t(<I(t)I(0)>-<I>^2)/<I>^2"+lineSep;

        try {
            FileWriter file = new FileWriter(saveNameCorr);
            file.write(buffer);

            // logarithmically spaced lag times
            int dec = (int) (Math.log(nFrames)/Math.log(10));
            int t0=0,cnt=1;
            for(int pow=0; pow<=20*dec;pow++){
                int t1 = (int) Math.exp(pow*Math.log(10)/20.0);
                if(t1!=t0 && t1<nFrames){
                    cnt++;
                    t0=t1;
                }
            }

            int[] t_log = new int[cnt];

            t0=0;
            cnt=1;
            for(int pow=0; pow<=20*dec;pow++){
                int t1 = (int) Math.exp(pow*Math.log(10)/20);
                if(t1!=t0 && t1<nFrames){
                    t_log[cnt]=t1;
                    cnt++;
                    t0=t1;
                }
            }


            for(int t:t_log){
                buffer = ""+(t*dt)+"\t"+CI[t]+"\t"+CI[t]/CI[0]+"\t"+CI[t]/(Im*Im)+lineSep;
                file.write(buffer);
            }

            file.close();
        } catch (Exception e){
            IJ.log("Erreur doSave --> "+e.getMessage());
            IJ.log("Erreur doSave --> "+e.getCause());
            IJ.log("Erreur doSave --> "+e.getLocalizedMessage());
        }

        IJ.showStatus("Done");
    }

    private void saveParams(){
        String buffer = "Parameters Brownian FCS"+lineSep;

        buffer +="Simu_duration (s)\t"+Tmax+lineSep;
        buffer +="Time_step (s)\t"+dt+lineSep;
        buffer +="Number_of_Objects\t"+nObjects+lineSep;
        buffer +="Diffusion_coeff (um2/s)\t"+D+lineSep;

        buffer +="Cell_length (um)\t"+L+lineSep;
        buffer +="Cell_diameter (um)\t"+d+lineSep;


        buffer +="w0 (um)\t"+w0+lineSep;
        buffer +="z0 (um)\t"+z0+lineSep;
        buffer +="xc (um)\t"+xc+lineSep;
        buffer +="yc (um)\t"+yc+lineSep;
        buffer +="zc (um)\t"+zc+lineSep;


        try {
            FileWriter file = new FileWriter(saveNameParam);
            file.write(buffer);

            file.close();
        } catch (Exception e){
            IJ.log("Erreur doSave --> "+e.getMessage());
            IJ.log("Erreur doSave --> "+e.getCause());
            IJ.log("Erreur doSave --> "+e.getLocalizedMessage());
        }

        IJ.showStatus("Done");
    }

    private double intensityFCS(double[] r){

        double dr = ( (r[0]-xc)*(r[0]-xc) + (r[1]-yc)*(r[1]-yc) )/w0r + ( (r[2]-zc)*(r[2]-zc)  )/z0r ;

        return Math.exp( - dr ); // (dr<9)?Math.exp( - dr ):0

    }

    private boolean correctStartPosition(int k){
        double x,y,z,X;

        pos[k][0]=L*rd.nextDouble()-L/2;
        pos[k][1]=d*rd.nextDouble()-d/2;
        pos[k][2]=d*rd.nextDouble()-d/2;

        // boundary check
        x = pos[k][0];
        y = pos[k][1];
        z = pos[k][2];

        X = Math.abs(x)-a;
        X = (X>0)?X:0;

        return (X*X+y*y+z*z)<=b*b;

    }

    private double[] toComplex(double[] x){

        double[] out = new double[4*x.length]; // zero padding

        for(int k=0;k<x.length; k++){
            out[2*k]=x[k];

        }

        return out;
    }

}
