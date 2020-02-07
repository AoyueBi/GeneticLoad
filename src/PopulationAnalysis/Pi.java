package PopulationAnalysis;

import AoUtils.SplitScript;
import pgl.utils.IOUtils;

import java.io.File;
import java.util.List;

public class Pi {
    public Pi(){
//        this.mkPiCommandbasedwinndow();
        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/001_script_based2Mwindow_1Mstep/pi_based2Mwindow_1Mstep_20200207.sh",23,4); //91cmd

    }

    public void mkPiCommandbasedwinndow(){
        String windowsize = "2000000";
        String windowstep = "1000000";

        //HPC path
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/021_popGen/102_Pi/002_pi_based2Mwindow_1Mstep";
        String groupDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/102_Pi/000_group";

        //local path
        String grouplocalDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/000_group";

        List<File> fs = IOUtils.getVisibleFileListInDir(grouplocalDirS);
        for (int i = 0; i < fs.size(); i++) {
            String groupname = fs.get(i).getName().split(".txt")[0];
            if (groupname.equals("Cultivar") || groupname.equals("Landrace")){
                String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B", "1D", "2D", "3D",  "4D", "5D", "6D", "7D"};
                for (int j = 0; j < chrArr.length; j++) {
                    String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + windowsize+ "Window_" + windowstep + "step").getAbsolutePath();
                    String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --keep " + groupS +
                            " --window-pi " + windowsize +  " --window-pi-step 1000000  --out " + outfileS );
                }

            }
            else if (groupname.equals("Domesticated_emmer") || groupname.equals("Free_threshing_tetraploid") || groupname.equals("Wild_emmer")){
                String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                for (int j = 0; j < chrArr.length; j++) {
                    String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + windowsize+ "Window_" + windowstep + "step").getAbsolutePath();
                    String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --keep " + groupS +
                            " --window-pi " + windowsize +  " --window-pi-step 1000000  --out " + outfileS );
                }

            }
            else if (groupname.equals("Ae.tauschii")){
                String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
                for (int j = 0; j < chrArr.length; j++) {
                    String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + windowsize+ "Window_" + windowstep + "step").getAbsolutePath();
                    String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --keep " + groupS +
                            " --window-pi " + windowsize +  " --window-pi-step 1000000  --out " + outfileS );
                }
            }

        }
    }
}
