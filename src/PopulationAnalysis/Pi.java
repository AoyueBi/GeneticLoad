package PopulationAnalysis;

import AoUtils.AoFile;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Pi {
    public Pi(){
//        this.mkPiCommandbasedwinndow();
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/001_script_based2Mwindow_1Mstep/pi_based2Mwindow_1Mstep_20200207.sh",23,4); //91cmd
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/004_script_based100kbwindow_50kbstep/pi_based100kbwindow_50kbstep_20200213.sh",23,4); //91cmd

//        this.getMeanPIvalue();
//        this.mkPiCommandbasedwinndow_LandrcaeSub();
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/007_script_based2Mwindow_1Mstep_landraceWAEUEA/pi_based2Mwindow_1Mstep.sh",34,5); //168
        this.getMeanPIvaluefromLandraceWAEAEU();
    }

    /**
     * 获取 Landrace细分类后的PI值
     *
     */
    public void getMeanPIvaluefromLandraceWAEAEU(){
        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/groupInfo.txt";
        HashMap<String,String> hm = new AoFile().getHashMapStringKey(taxaList,0,1);

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/008_pi_based2Mwindow_1Mstep";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/009_meanPI/Pi_bySubspecies_20200208.txt";


        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tSub\tGroup\tPloidy\tMEAN_PI");
            bw.newLine();
            for (int i = 0; i < fsList.size(); i++) {
                String infileS = fsList.get(i).getAbsolutePath();
                String name = new File(infileS).getName(); //Cultivar_chr1D_based2000000Window_1000000step.windowed.pi
                String chr = name.substring(name.indexOf("chr")+3,name.indexOf("chr")+5);
                String sub = chr.substring(1);
                String group = name.substring(0,name.indexOf("_chr"));
                String ploidy = hm.get(group);
                String value = this.getMean(infileS);
                bw.write(chr + "\t" + sub + "\t" + group + "\t" + ploidy + "\t" + value);
                bw.newLine();

            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     *  calculate the pi of Landrace_WA EU EA
     *   I am so tired
     *
     */
    public void mkPiCommandbasedwinndow_LandrcaeSub(){
        String windowsize = "2000000";
        String windowstep = "1000000";

//        String windowsize = "100000";
//        String windowstep = "50000";

        //HPC path
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/021_popGen/102_Pi/006_pi_based2Mwindow_1Mstep";
        String groupDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/102_Pi/000_group2";

        //local path
        String grouplocalDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/000_group2";

        List<File> fs = IOUtils.getVisibleFileListInDir(grouplocalDirS);
        for (int i = 0; i < fs.size(); i++) {
            String groupname = fs.get(i).getName().split(".txt")[0];
            if (groupname.contains("Landrace")){
                String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B", "1D", "2D", "3D",  "4D", "5D", "6D", "7D"};
                for (int j = 0; j < chrArr.length; j++) {
                    String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + windowsize+ "Window_" + windowstep + "step").getAbsolutePath();
                    String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --keep " + groupS +
                            " --window-pi " + windowsize +  " --window-pi-step " + windowstep + "  --out " + outfileS );
                }

            }
            else if (groupname.contains("Wild_emmer")){
                String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                for (int j = 0; j < chrArr.length; j++) {
                    String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + windowsize+ "Window_" + windowstep + "step").getAbsolutePath();
                    String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --keep " + groupS +
                            " --window-pi " + windowsize +  " --window-pi-step " + windowstep + "  --out " + outfileS );
                }

            }
//            else if (groupname.equals("Ae.tauschii")){
//                String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
//                for (int j = 0; j < chrArr.length; j++) {
//                    String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
//                    String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + windowsize+ "Window_" + windowstep + "step").getAbsolutePath();
//                    String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
//                    System.out.println("vcftools --vcf " + infileS + " --keep " + groupS +
//                            " --window-pi " + windowsize +  " --window-pi-step " + windowstep +  "  --out " + outfileS );
//                }
//            }

        }
    }


    /**
     * get mean from a file， is the value equals -nan, continue
     *
     * @param infileS
     * @return
     */
    public String getMean(String infileS){
        String out = null;
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            TDoubleArrayList vList = new TDoubleArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String t = l.get(4);
                if (t.startsWith("-n") || t.startsWith("n")) { //pos中 fst值为 -nan 的全部去掉
                    continue;
                }
                double v = Double.parseDouble(t);
                vList.add(v);
                cnt++; //有 fst值的位点数目
            }
            br.close();
            double[] v = vList.toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(v);
            double m = d.getMean();
            out = String.format("%.6f", m);
            System.out.println(infileS + " with mean pi " + out + " is completed");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return out;
    }

    /**
     *
     */
    public void getMeanPIvalue(){
        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
        HashMap<String,String> hm = new AoFile().getHashMapStringKey(taxaList,11,8);

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/002_pi_based2Mwindow_1Mstep";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/003_meadPI/Pi_bySubspecies_20200208.txt";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/008_pi_based2Mwindow_1Mstep";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/009_meanPI/Pi_bySubspecies_20200208.txt";


        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tSub\tGroup\tPloidy\tMEAN_PI");
            bw.newLine();
            for (int i = 0; i < fsList.size(); i++) {
                String infileS = fsList.get(i).getAbsolutePath();
                String name = new File(infileS).getName(); //Cultivar_chr1D_based2000000Window_1000000step.windowed.pi
                String chr = name.substring(name.indexOf("chr")+3,name.indexOf("chr")+5);
                String sub = chr.substring(1);
                String group = name.substring(0,name.indexOf("_chr"));
                String ploidy = hm.get(group);
                String value = this.getMean(infileS);
                bw.write(chr + "\t" + sub + "\t" + group + "\t" + ploidy + "\t" + value);
                bw.newLine();

            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkPiCommandbasedwinndow(){
//        String windowsize = "2000000";
//        String windowstep = "1000000";

        String windowsize = "100000";
        String windowstep = "50000";

        //HPC path
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
//        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/021_popGen/102_Pi/002_pi_based2Mwindow_1Mstep";
        String groupDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/102_Pi/000_group";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/021_popGen/102_Pi/003_pi_based100kbwindow_50kbstep";

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
                            " --window-pi " + windowsize +  " --window-pi-step " + windowstep + "  --out " + outfileS );
                }

            }
            else if (groupname.equals("Domesticated_emmer") || groupname.equals("Free_threshing_tetraploid") || groupname.equals("Wild_emmer")){
                String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                for (int j = 0; j < chrArr.length; j++) {
                    String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + windowsize+ "Window_" + windowstep + "step").getAbsolutePath();
                    String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --keep " + groupS +
                            " --window-pi " + windowsize +  " --window-pi-step " + windowstep + "  --out " + outfileS );
                }

            }
            else if (groupname.equals("Ae.tauschii")){
                String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
                for (int j = 0; j < chrArr.length; j++) {
                    String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + windowsize+ "Window_" + windowstep + "step").getAbsolutePath();
                    String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --keep " + groupS +
                            " --window-pi " + windowsize +  " --window-pi-step " + windowstep +  "  --out " + outfileS );
                }
            }

        }
    }
}
