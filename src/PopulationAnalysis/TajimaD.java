package PopulationAnalysis;

import AoUtils.AoFile;
import AoUtils.AoWinScan;
import AoUtils.SplitScript;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class TajimaD {
    public TajimaD(){
//        this.mkTajimaDCommandbasedwinndow();
//        this.window();
//        this.addGroupToTajimaDwindow();
//        this.getMeanTajimaDvalue();

//        this.changePostoRef();

        this.window();


    }

    public void changePostoRef(){
//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/011_populationPara/004_thetaW/002_merge001/angsd_subspecies26_geneRegion.txt.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/011_populationPara/004_thetaW/002_merge001/angsd_subspecies26_geneRegion_RefChr.txt.gz";

        String infileS = "/Users/Aoyue/Documents/df4.txt";
        String outfileS = "/Users/Aoyue/Documents/df4_RefChr.txt";

        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //read header
            bw.write( br.readLine() + "\tRefChr\tRefPos");
            bw.newLine();
            String temp = null; //read header
            int cnt = 0;
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                int chrID = Integer.parseInt(l.get(3));
                int pos = Integer.parseInt(l.get(4));
                String chr = RefV1Utils.getChromosome(chrID,pos);
                int posOnChrosome = RefV1Utils.getPosOnChromosome(chrID,pos);
                //先找到 chr 所在的列
                StringBuilder sb = new StringBuilder();
                sb.append(temp).append("\t").append(chr).append("\t").append(posOnChrosome);
                bw.write(sb.toString());
                bw.newLine();
            }

            br.close();
            bw.flush();
            bw.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     *
     */
    public void getMeanTajimaDvalue(){
        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        HashMap<String,String> hm = new AoFile().getHashMapStringKey(taxaList,15,3);//亚种和倍性之间的hashmap
//        String infileDirS = "";
//        String outfileDirS = "";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/003_TajimaD/002_TajimaD_based50000window/001";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/003_TajimaD/002_TajimaD_based50000window/004/001_TajimaD_bySubspecies_20200907.txt";


        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tSub\tGroup\tPloidy\tMEAN_TajimaD");
            bw.newLine();
            for (int i = 0; i < fsList.size(); i++) {
                String infileS = fsList.get(i).getAbsolutePath();
                String name = new File(infileS).getName(); //Cultivar_chr1D_based2000000Window_1000000step.windowed.pi
                String chr = name.substring(name.indexOf("chr")+3,name.indexOf("chr")+5);
                String sub = chr.substring(1);
                String group = name.substring(0,name.indexOf("_chr"));
                String ploidy = hm.get(group);
                String value = this.getMean(infileS); //调用函数
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
     * get mean from a file， if the value equals -nan, continue
     *
     * @param infileS
     * @return
     */
    private String getMean(String infileS){
        String out = null;
        try {
            BufferedReader br = AoFile.readFile(infileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            TDoubleArrayList vList = new TDoubleArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String t = l.get(3);
                if (t.startsWith("-n") || t.startsWith("n")) { //pos中 pi值为 -nan 的全部去掉
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


    public void addGroupToTajimaDwindow(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/005_TajimaD_windowbyJava";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/006_TajimaD/chr4A.windowed.weir.TajimaD.txt.gz";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/003_TajimaD/002_TajimaD_based50000window/002_windowbyJava";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/003_TajimaD/002_TajimaD_based50000window/003_merge/001_TajimaD_2Mwindow_1Mstep.txt";

        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        Arrays.sort(fs);
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            /**
             * 需要改动,header的名字可以自定义
             */
            //read header
            bw.write(br.readLine() + "\tGroup");
            bw.newLine();

            int cnttotal = 0;
            //read context
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                /**
                 * 需要改动
                 */
                String name = fs[i].getName();
                String group = name.split("_chr")[0];

                br = AoFile.readFile(infileS);
                br.readLine(); //read header
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp).append("\t").append(group);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(fs[i].getName() + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void window(){

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/004_TajimaD";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/005_TajimaD_windowbyJava";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/003_TajimaD/002_TajimaD_based50000window/001";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/003_TajimaD/002_TajimaD_based50000window/002_windowbyJava";

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            AoFile.readheader(infileS);
            int chrColumn = 0;
            int posIndex = 1;
            int valueIndex = 3;
//            double window = 2000000;
//            double step = 1000000;
            double window = 10000000;
            double step = 1000000;
            String name = new File(infileS).getName().split(".TajimaD")[0] + "_" + window + "window_" + step + "step.txt.gz";
            String parent = new File(infileS).getParent();
            String outfileS = new File(outfileDirS,name).getAbsolutePath();

            new AoWinScan().getwindowDistrbution_general(infileS,chrColumn,posIndex,valueIndex,window,step,outfileS);
        });
    }


    public void mkTajimaDCommandbasedwinndow(){
// ******************************************* 参数设置模板 *********************************************************************
//        //local path
//        String grouplocalDirS = "";
//        //HPC path
//        String infileDirS = "";
//        String groupDirS = "";
//        String outfileDirS = "";
//        //para
//        String window = "";
//        String step = "";
        //        String scriptS = "";

// ******************************************* 参数设置 *********************************************************************
        /**
         * first run
         */
        //HPC path
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
//        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/021_popGen/102_Pi/002_pi_based2Mwindow_1Mstep";
//        String groupDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/102_Pi/000_group";
//        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/021_popGen/102_Pi/003_pi_based100kbwindow_50kbstep";
        //local path
//        String grouplocalDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/000_group";
        //        String window = "2000000";
//        String step = "1000000";
//        String window = "100000";
//        String step = "50000";

// ******************************************* 参数设置 *********************************************************************
        /**
         * VMap2.0 after May 1st
         */
        //local path
//        String grouplocalDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/003_TajimaD/000_group";
//        //HPC path
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/105_VMap2.1ByRef";
//        String groupDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/003_TajimaD/000_group";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/003_TajimaD/002_TajimaD_based50000window/001";
//        //para
//        String window = "50000";
//        int numfile = 10;
//        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/003_TajimaD/001_srcipt_based50000window/tajimaD_based" + window + "window" + "_20200904.sh";


// ******************************************* 参数设置 *********************************************************************
//        //****** local path
//        String grouplocalDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/011_populationPara/003_TajimaD/000_group";
//        //****** HPC path
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/203_VMap2.1ByRef";
//        String groupDirS = "/data4/home/aoyue/vmap2/analysis/049_PopGenetics/003_TajimaD/000_group";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/049_PopGenetics/003_TajimaD/002_TajimaD_based50000window";
//        //****** para
//        String window = "50000";
//        String step = "";
//        String scriptS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/011_populationPara/003_TajimaD/001_script/001/tajimaD_based" + window + "window" + "_2021_10_22_1.sh";

// ******************************************* 参数设置 *********************************************************************
//        //****** local path
//        String grouplocalDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/011_populationPara/003_TajimaD/000_group2";
//        //****** HPC path
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/203_VMap2.1ByRef";
//        String groupDirS = "/data4/home/aoyue/vmap2/analysis/049_PopGenetics/003_TajimaD/000_group2";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/049_PopGenetics/003_TajimaD/004_TajimaD_based50000window";
//        //****** para
//        String window = "50000";
//        String step = "";
//        String scriptS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/011_populationPara/003_TajimaD/001_script/001/tajimaD_based" + window + "window" + "_2021_10_22_2.sh";

//// ******************************************* 参数设置 *********************************************************************
////        //****** local path
//        String grouplocalDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/011_populationPara/003_TajimaD/000_group3";
//        //****** HPC path
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/203_VMap2.1ByRef";
//        String groupDirS = "/data4/home/aoyue/vmap2/analysis/049_PopGenetics/003_TajimaD/000_group3";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/049_PopGenetics/003_TajimaD/006_TajimaD_based50000window";
//        //****** para
//        String window = "50000";
//        String step = "";
//        String scriptS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/011_populationPara/003_TajimaD/001_script/001/tajimaD_based" + window + "window" + "_2021_10_22_3.sh";


        // ******************************************* 参数设置 *********************************************************************
//        //****** local path
        String grouplocalDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/016_1B1R_translocations/001_popPara/000_group";
        //****** HPC path
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/203_VMap2.1ByRef";
        String groupDirS = "/data4/home/aoyue/vmap2/analysis/051_1B1R_translocations/001_popPara/000_group";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/051_1B1R_translocations/001_popPara/003_tajimaD";
        //****** para
        String window = "50000";
        String step = "";
        String scriptS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/016_1B1R_translocations/001_popPara/script/tajimaD_based" + window + "window" + "_2021_10_22_3.sh";

// ******************************************* 参数设置完毕后的 section *********************************************************************
        System.out.println("mkdir 001_srcipt_based" + window + "window" );
        System.out.println("mkdir 002_TajimaD_based" + window + "window");

        List<File> fs = IOUtils.getVisibleFileListInDir(grouplocalDirS);
        try{
            BufferedWriter bw = AoFile.writeFile(scriptS);
            for (int i = 0; i < fs.size(); i++) {
                String groupname = fs.get(i).getName().split(".txt")[0];
                if (groupname.equals("With1B1R") || groupname.equals("No1B1R")){
//                if (groupname.equals("compactum") || groupname.equals("LR_AF") || groupname.equals("LR_AM") || groupname.equals("LR_CSA") || groupname.equals("LR_WA") || groupname.equals("LR_EA") || groupname.equals("LR_EU")
//                        || groupname.equals("Cultivar") || groupname.equals("macha") || groupname.equals("OtherHexaploids") || groupname.equals("petropavlovskyi") || groupname.equals("spelta") || groupname.equals("sphaerococcum")
//                        || groupname.equals("tibeticum") || groupname.equals("vavilovii") || groupname.equals("yunna-nense")){
//                if (groupname.equals("LR_AF") || groupname.equals("LR_AM") || groupname.equals("LR_CSA") || groupname.equals("LR_WA") || groupname.equals("LR_EA") || groupname.equals("LR_EU") || groupname.equals("CL")){
//                if (groupname.equals("CL") || groupname.equals("LR_EU") || groupname.equals("LR_EA")){
//                if (groupname.equals("Cultivar") || groupname.equals("Landrace")){
//                    String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B", "1D", "2D", "3D",  "4D", "5D", "6D", "7D"};
                    String[] chrArr = {"1A", "1B","1D"};

                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + window+ "Window_").getAbsolutePath();
                        String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --keep ").append(groupS);
                        sb.append(" --TajimaD ").append(window);
                        sb.append("  --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }

                }
                else if (groupname.equals("carthlicum") || groupname.equals("dicoccoides") || groupname.equals("dicoccum") || groupname.equals("durum") || groupname.equals("ispahanicum")
                        || groupname.equals("karamyschevii") || groupname.equals("polonicum") || groupname.equals("turanicum") || groupname.equals("turgidum")){
//                else if (groupname.equals("DE") || groupname.equals("FTT") || groupname.equals("WE")){
//                else if (groupname.equals("Domesticated_emmer") || groupname.equals("Free_threshing_tetraploid") || groupname.equals("Wild_emmer")){
                    String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + window+ "Window_").getAbsolutePath();
                        String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --keep ").append(groupS);
                        sb.append(" --TajimaD ").append(window);
                        sb.append("  --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }

                }
                else if (groupname.equals("strangulata")){
//                else if (groupname.equals("AT")){
//                else if (groupname.equals("Ae.tauschii")){
                    String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + window+ "Window_").getAbsolutePath();
                        String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --keep ").append(groupS);
                        sb.append(" --TajimaD ").append(window);
                        sb.append("  --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }
                }

            }
            bw.flush();bw.close();
//            SplitScript.splitScript4(scriptS,numfile); //脚本拆分

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

}
