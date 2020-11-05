package PopulationAnalysis;

import AoUtils.AoFile;
import AoUtils.AoWinScan;
import AoUtils.SplitScript;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class Pi {
    public Pi(){
//        this.mkPiCommandbasedwinndow();
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/001_script_based2Mwindow_1Mstep/pi_based2Mwindow_1Mstep_20200207.sh",23,4); //91cmd
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/004_script_based100kbwindow_50kbstep/pi_based100kbwindow_50kbstep_20200213.sh",23,4); //91cmd

//        this.getMeanPIvalue();
//        this.mkPiCommandbasedwinndow_LandraceSub();
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/007_script_based2Mwindow_1Mstep_landraceWAEUEA/pi_based2Mwindow_1Mstep.sh",34,5); //168
//        this.getMeanPIvaluefromLandraceWAEAEU();


//                this.mkPiCommandbasedwinndow();
//                this.getMeanPIvalue();
//        this.window();
//        this.addGroupToPiwindow();

//        this.mkPiCMD_by2File();
//                this.window();
                this.piRatio();


    }

    /**
     * 求2条染色体的pi ratio
     */
    public void piRatio(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/007_tetraploid_Dm_DE/010_test3_pifst/002_pi/003_pi_window_byaoCode/Domesticated_emmer_chr1A_basedSNP.sites_2000.0window_2000.0step.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/007_tetraploid_Dm_DE/010_test3_pifst/002_pi/003_pi_window_byaoCode/Durum_chr1A_basedSNP.sites_2000.0window_2000.0step.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/007_tetraploid_Dm_DE/010_test3_pifst/002_pi/004_piRatio/chr1A_piRatio.txt";

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedReader br2 = AoFile.readFile(infileS2);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("CHROM\tBIN_START\tBIN_END\tRatio"); bw.newLine();
            String temp = br.readLine();
            String temp2 = br2.readLine();
            int cnt = 0;
            List<String> l = new ArrayList<>();
            List<String> l2 = new ArrayList<>();
            double ratio = Double.NaN;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                temp2 = br2.readLine();
                l = PStringUtils.fastSplit(temp);
                l2 = PStringUtils.fastSplit(temp2);
                String pi1 = l.get(4);
                String pi2 = l2.get(4);
                if (pi1.startsWith("N"))continue;
                if (pi2.startsWith("N"))continue;
                ratio = (double) Double.parseDouble(pi1)/Double.parseDouble(pi2);
                System.out.println(ratio);
                sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(ratio);
                bw.write(sb.toString());
                bw.newLine();
                cnt++;
//                if (cnt > 20) break;
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(cnt + " lines in the file.");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    public void mkPiCMD_by2File(){
        //local path
        String grouplocalDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/002_Pi/000_group";

        String pop1FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_bySubspecies/Durum.txt";
        String pop2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_bySubspecies/Domesticated_emmer.txt";


        //HPC path
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/007_tetraploid_Dm_DE/000_exonVCF_ref";
        String groupFileDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/000_pop_bySubspecies";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/007_tetraploid_Dm_DE/013_fstTest/002_pi";
        //para
//        String window = "100000";
//        String step = "50000";
        int numcmd = 28;
//        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/007_tetraploid_Dm_DE/010_test3_pifst/002_script_Pi_basedSNP/pi_based" + window + "window_" + step + "step_20200904.sh";
        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/007_tetraploid_Dm_DE/010_test3_pifst/002_script_Pi_basedSNP/pi_basedSNP" + "_20201105.sh";

        System.out.println("mkdir 001_srcipt_basedSNP");
        System.out.println("mkdir 002_pi_basedSNP");

//        System.out.println("mkdir 001_srcipt_based" + window + "window_" + step + "step");
//        System.out.println("mkdir 002_pi_based" + window + "window_" + step + "step");

//        List<File> fs = IOUtils.getVisibleFileListInDir(grouplocalDirS);

        List<File> fs = new ArrayList<>();
        fs.add(new File(pop1FileS)); fs.add(new File(pop2FileS));

        try{
            BufferedWriter bw = AoFile.writeFile(scriptS);
            for (int i = 0; i < fs.size(); i++) {
                String groupname = fs.get(i).getName().split(".txt")[0];
                    String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1.vcf.gz").getAbsolutePath();
//                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + window+ "Window_" + step + "step").getAbsolutePath();
                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_basedSNP").getAbsolutePath();
                        String groupS = new File(groupFileDirS,groupname+".txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --keep ").append(groupS);
                        sb.append(" --site-pi");
//                        sb.append(" --window-pi ").append(window).append(" --window-pi-step ").append(step);
                        sb.append("  --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }
            }
            bw.flush();bw.close();
            SplitScript.splitScript4(scriptS,numcmd); //脚本拆分

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


    }


    /**
     * 向window滑窗的结果添加group信息，区分不同分组的比较
     */
    public void addGroupToPiwindow(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/002_Pi/002_pi_based100000window_50000step/003_window";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/002_Pi/002_pi_based100000window_50000step/004_merge/001_Pi_2Mwindow_1Mstep.txt";
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
        // model
//        String infileDirS = "";
//        String outfileDirS = "";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/007_pi";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/008_pi/";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/002_Pi/002_pi_based100000window_50000step/001";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/002_Pi/002_pi_based100000window_50000step/003_window";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/007_tetraploid_Dm_DE/010_test3_pifst/002_pi/002_pi_output";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/007_tetraploid_Dm_DE/010_test3_pifst/002_pi/003_pi_window_byaoCode";

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            AoFile.readheader(infileS);
//            int chrColumn = 0;
//            int posIndex = 1;
//            int valueIndex = 4;
//            double window = 2000000;
//            double step = 1000000;

            int chrColumn = 0;
            int posIndex = 1;
            int valueIndex = 2;
            double window = 2000;
            double step = 2000;
            String name = new File(infileS).getName().split(".pi")[0] + "_" + window + "window_" + step + "step.txt";
            String parent = new File(infileS).getParent();
            String outfileS = new File(outfileDirS,name).getAbsolutePath();

            new AoWinScan().getwindowDistrbution_general(infileS,chrColumn,posIndex,valueIndex,window,step,outfileS);
        });

    }


    /**
     * 根据4A染色体的几个亚群，做一张综合的表
     */
    public void getPiRatio(){
        String infileDirS = "";
        String outfileS = "";

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
     *
     *
     */
    public void mkPiCommandbasedwinndow_LandraceSub(){
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
     * get mean from a file， if the value equals -nan, continue
     *
     * @param infileS
     * @return
     */
    private String getMean(String infileS){
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

    /**
     *
     */
    public void getMeanPIvalue(){
        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        HashMap<String,String> hm = new AoFile().getHashMapStringKey(taxaList,15,3);//亚种和倍性之间的hashmap
//        String infileDirS = "";
//        String outfileS = "";

//        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
//        HashMap<String,String> hm = new AoFile().getHashMapStringKey(taxaList,11,8);//亚种和倍性之间的hashmap

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/002_pi_based2Mwindow_1Mstep";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/003_meadPI/Pi_bySubspecies_20200208.txt";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/008_pi_based2Mwindow_1Mstep";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/009_meanPI/Pi_bySubspecies_20200208.txt";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/002_Pi/002_pi_based100000window_50000step/001";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/002_Pi/002_pi_based100000window_50000step/002/001_Pi_bySubspecies_20200907.txt";

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

        /**
         * VMap2.0 after May 1st
         */
        //local path
        String grouplocalDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/002_Pi/000_group";
        //HPC path
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/105_VMap2.1ByRef";
        String groupDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/002_Pi/000_group";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/002_Pi/002_pi_based100000window_50000step/001";
        //para
        String window = "100000";
        String step = "50000";
        int numcmd = 5;
        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/002_Pi/001_srcipt_based100000window_50000step/pi_based" + window + "window_" + step + "step_20200904.sh";



        System.out.println("mkdir 001_srcipt_based" + window + "window_" + step + "step");
        System.out.println("mkdir 002_pi_based" + window + "window_" + step + "step");

        List<File> fs = IOUtils.getVisibleFileListInDir(grouplocalDirS);
        try{
            BufferedWriter bw = AoFile.writeFile(scriptS);
            for (int i = 0; i < fs.size(); i++) {
                String groupname = fs.get(i).getName().split(".txt")[0];
                if (groupname.equals("Cultivar") || groupname.equals("Landrace")){
                    String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B", "1D", "2D", "3D",  "4D", "5D", "6D", "7D"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + window+ "Window_" + step + "step").getAbsolutePath();
                        String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --keep ").append(groupS);
                        sb.append(" --window-pi ").append(window).append(" --window-pi-step ").append(step);
                        sb.append("  --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }

                }
                else if (groupname.equals("Domesticated_emmer") || groupname.equals("Free_threshing_tetraploid") || groupname.equals("Wild_emmer")){
                    String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + window+ "Window_" + step + "step").getAbsolutePath();
                        String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --keep ").append(groupS);
                        sb.append(" --window-pi ").append(window).append(" --window-pi-step ").append(step);
                        sb.append("  --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }

                }
                else if (groupname.equals("Ae.tauschii")){
                    String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
                    for (int j = 0; j < chrArr.length; j++) {
                        String infileS = new File(infileDirS,"chr" + chrArr[j] + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS,groupname + "_chr" + chrArr[j] + "_based" + window+ "Window_" + step + "step").getAbsolutePath();
                        String groupS = new File(groupDirS,groupname+".txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --keep ").append(groupS);
                        sb.append(" --window-pi ").append(window).append(" --window-pi-step ").append(step);
                        sb.append("  --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }
                }

            }
            bw.flush();bw.close();
            SplitScript.splitScript4(scriptS,numcmd); //脚本拆分

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
