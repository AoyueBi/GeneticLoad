package PopulationAnalysis;

import AoUtils.AoMath;
import AoUtils.CountSites;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Fst {
    public Fst() { //计算不消耗内存
//        this.mkFstCommandbasedSNP();
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/001_scriptSNPbased/fst_basedSNP_20200205.sh",7,23);

//        this.mkFstTable("/Users/Aoyue/Documents/po1_VS_pop2_chr1A.txt","/Users/Aoyue/Documents/out.txt");
//        this.scriptMkFstTable();

//        this.mkFstCommandbasedwinndow();
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/003_scriptbased2Mwindow1Mstep/sh_fst_based2Mwindow_1Mstep_20200205.sh",21,8);
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/005_script_based100kwindow_50kstep/fst_based100kwindow_50kstep_20200213.sh",21,8);

//        this.extractVCFlog();
//        this.mergeTxt();

//        this.mkFstCommandbasedwinndow_newGroup();
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/006_scriptbased2Mwindow1Mstep/fst_based2Mwindow_1Mstep_20200216.sh",40,7); //273

//        this.extractVCFlog();
//        this.addSubandGroup();

        this.mergeFSTwindow();
    }

    /**
     * 将fst window scan 计算的FST,多个文件合并起来成为一个文件，并添加一列分组Group信息
     * 暂时
     *
     */
    public void mergeFSTwindow(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/004_fst_based2Mwindow_1Mstep/001";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/004_fst_based2Mwindow_1Mstep/002_merge/Pi_bySubspecies_20200208.txt";

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
                String value = this.getMean(infileS);
                bw.write(chr + "\t" + sub + "\t" + group + "\t" + "jj" + "\t" + value);
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
     * 将landrace细化分组，为欧洲的landrace和东亚的landrace，比较两两之间的关系
     */
    public void mkFstCommandbasedwinndow_newGroup() {

        //local file： one: group
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group2/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group2/hexaandDi";

        // HPC file: group fileDirS
        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group2/hexaandTetra";
        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group2/hexaandDi";

        // HPC file: output fileDirS
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/008_fst_based2Mwindow_1Mstep";

        List<File> fs = IOUtils.getVisibleFileListInDir(groupHexaandTetraDirS);
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        List<File> fs2 = IOUtils.getVisibleFileListInDir(groupHexaandDiDirS);
        File[] group2FileS = fs2.toArray(new File[fs2.size()]);

        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。

        for (int i = 0; i < group1FileS.length - 1; i++) {
            String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
            for (int j = i + 1; j < group1FileS.length; j++) {
                String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group1FileDirS, group1FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group1FileDirS, group1FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S +
                            " --weir-fst-pop " + group2S + " --fst-window-size 2000000 --fst-window-step 1000000 " +  " --out " + outfileS);
//                            " --weir-fst-pop " + group2S + " --fst-window-size 100000 --fst-window-step 50000 " +  " --out " + outfileS);

                }
            }
        }

        for (int i = 0; i < group2FileS.length - 1; i++) {
            String pop1 = group2FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group2FileS.length; j++) {
                String pop2 = group2FileS[j].getName().replace(".txt", "");
                String[] chrArr = {"1D", "2D", "3D", "4D", "5D", "6D", "7D"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group2FileDirS, group2FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group2FileDirS, group2FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S +
                            " --weir-fst-pop " + group2S + " --fst-window-size 2000000 --fst-window-step 1000000 " +  " --out " + outfileS);
//                            " --weir-fst-pop " + group2S + " --fst-window-size 100000 --fst-window-step 50000 " +  " --out " + outfileS);

                }
            }
        }
    }


//    public void mergeTxt(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/004_fst_based2Mwindow_1Mstep/test/001";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/004_fst_based2Mwindow_1Mstep/test/002_merge/Free_threshing_tetraploid_VS_Landrace_windowed.weir.fst.txt";
//        new CountSites().mergeTxt(infileDirS,outfileS);
//    }

    public void mergeTxt(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/test/003_merge";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/102_Pi/test/004/pi.txt";
        new CountSites().mergeTxt(infileDirS,outfileS);
    }

    /**
     * 为两两群体之间的FST添加①亚基因组的信息和②分组的信息
     *
     */
    public void addSubandGroup(){
        String infileS="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/log/004_merge/001_Fst_bySubspecies_LandraceEUEA_20200217_ori.txt";
        String outfileS="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/log/004_merge/002_Fst_bySubspecies_LandraceEUEA_addSub.txt";
        HashMap<String,Integer> hm = new AoMath().setGrouptoNumber(infileS,1);

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);

            String temp = br.readLine();
            bw.write(temp + "\tGroup\tSub");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String chr = l.get(0);
                String sub = chr.substring(1);
                String name = l.get(1);
                int group = hm.get(name);
                for (int i = 0; i < l.size() ; i++) {
                    bw.write(l.get(i) + "\t");
                }
                bw.write(group + "\t" + sub);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void extractVCFlog(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/log/001_log";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/log/002_merge/Fst_bySubspecies_20200208.txt";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/log/003";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/log/004_merge/001_Fst_bySubspecies_LandraceEUEA_20200217_ori.txt";
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tPOP1_VS_POP2\tWEIGHTED_FST\tMEAN_FST");
            bw.newLine();
            for (int i = 0; i < fsList.size(); i++) {
                String infileS = fsList.get(i).getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(infileS);

                String temp = null;
                while ((temp = br.readLine()) != null) {
                    /**
                     * start to deal with: extract chr pop1 pop2
                     */

                    if (temp.equals("Parameters as interpreted:")){
                        for (int j = 0; j < 16; j++) {
                            temp=br.readLine();
                            if (j==7){
                                String outPath = PStringUtils.fastSplit(temp).get(1).replaceFirst("--out ","");
//                                System.out.println(outPath);
                                String name = new File(outPath).getName(); //Cultivar_VS_Domesticated_emmer_chr1A
                                String chr = name.substring(name.indexOf("chr")+3,name.indexOf("chr")+5);
                                String pop1pop2 = name.substring(0,name.indexOf("_chr"));
                                bw.write(chr + "\t" + pop1pop2 + "\t");
                            }
                            if (j==12){ //Weir and Cockerham mean Fst estimate: 0.14323
                                System.out.println(temp);
                                String meanValue = temp.replaceFirst("Weir and Cockerham mean Fst estimate: ","");
                                temp = br.readLine();
                                String weightedValue = temp.replaceFirst("Weir and Cockerham weighted Fst estimate: ","");
                                bw.write(weightedValue + "\t" + meanValue);
                                bw.newLine();
                            }
                            
                        } //16

                    }

                }
                br.close();
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkFstCommandbasedwinndow() {

        //local file： one: group
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandDi";

        // HPC file: group fileDirS
        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandTetra";
        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandDi";

        // HPC file: output fileDirS
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/004_fst_based2Mwindow_1Mstep/001";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/006_fst_based100kwindow_50kstep/001";

        List<File> fs = IOUtils.getVisibleFileListInDir(groupHexaandTetraDirS);
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        List<File> fs2 = IOUtils.getVisibleFileListInDir(groupHexaandDiDirS);
        File[] group2FileS = fs2.toArray(new File[fs2.size()]);

        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。

        for (int i = 0; i < group1FileS.length - 1; i++) {
            String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
            for (int j = i + 1; j < group1FileS.length; j++) {
                String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group1FileDirS, group1FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group1FileDirS, group1FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S +
//                            " --weir-fst-pop " + group2S + " --fst-window-size 2000000 --fst-window-step 1000000 " +  " --out " + outfileS);
                    " --weir-fst-pop " + group2S + " --fst-window-size 100000 --fst-window-step 50000 " +  " --out " + outfileS);

                }
            }
        }

        for (int i = 0; i < group2FileS.length - 1; i++) {
            String pop1 = group2FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group2FileS.length; j++) {
                String pop2 = group2FileS[j].getName().replace(".txt", "");
                String[] chrArr = {"1D", "2D", "3D", "4D", "5D", "6D", "7D"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group2FileDirS, group2FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group2FileDirS, group2FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S +
//                            " --weir-fst-pop " + group2S + " --fst-window-size 2000000 --fst-window-step 1000000 " +  " --out " + outfileS);
                    " --weir-fst-pop " + group2S + " --fst-window-size 100000 --fst-window-step 50000 " +  " --out " + outfileS);

                }
            }
        }
    }


    /**
     * 进行批量计算
     *
     */
    public void scriptMkFstTable() {

        //local file： one: group two: script
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandDi";

        // HPC file: output fileDirS
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/002_fst_basedSNP";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/003_fst_table";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/log/001";

        List<File> fs = IOUtils.getVisibleFileListInDir(groupHexaandTetraDirS);
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        List<File> fs2 = IOUtils.getVisibleFileListInDir(groupHexaandDiDirS);
        File[] group2FileS = fs2.toArray(new File[fs2.size()]);

        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。

        for (int i = 0; i < group1FileS.length - 1; i++) {
            String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
            for (int j = i + 1; j < group1FileS.length; j++) {
                String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".weir.fst").getAbsolutePath();
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".txt").getAbsolutePath();
                    String logS = new File(logDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".txt").getAbsolutePath();
                    System.out.println("java -jar 036_mkFstTable.jar " + infileS + " " + outfileS + " > " + logS + " 2>&1");
                }
            }
        }

        for (int i = 0; i < group2FileS.length - 1; i++) {
            String pop1 = group2FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group2FileS.length; j++) {
                String pop2 = group2FileS[j].getName().replace(".txt", "");
                String[] chrArr = {"1D", "2D", "3D", "4D", "5D", "6D", "7D"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".weir.fst").getAbsolutePath();
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".txt").getAbsolutePath();
                    String logS = new File(logDirS, pop1 + "_VS_" + pop2 + "_chr" + chr + ".txt").getAbsolutePath();
                    System.out.println("java -jar 036_mkFstTable.jar " + infileS + " " + outfileS + " > " + logS + " 2>&1");
                }
            }
        }
    }


    /**
     *
     *
     */
    public void mkFstTable(String infileS, String outfileS){
        //pseudo-code
        /**
         * 四倍体：WE DE FT
         * 六倍体：LR CU
         * 二倍体
         * 合计12种组合，两两之间。每种组合分为A B subgenome,每条Asub 分为 1-7条染色体
         *
         * 先跑每一种组合的7条，写成1A 2A 3A 4A 5A 6A 7A
         * 只要知道infileS,就能得出 chr group value , 就能输出到 outfileS 中；
         */

        try {
            String name = new File(infileS).getName();
            String chr = this.getChr(infileS);
            String pop1pop2 = name.substring(0,name.indexOf("_chr"));
            String value = this.getMean(infileS);

            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tPOP1_VS_POP2\tMEAN_FST");
            bw.newLine();
            bw.write(chr+ "\t" + pop1pop2 + "\t" + value);
            bw.newLine();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据文件的名字获取chr信息，不进行里面的读取
     *
     * @param infileS
     * @return
     */
    public String getChr(String infileS){
        String out = null;

        String name = new File(infileS).getName();
        int index = name.indexOf("chr");
        int index1 = index+3;
        int index2 = index+5;
        out = name.substring(index1,index2);
        System.out.println(out);

        return out;
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
                String t = l.get(2);
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
            out = String.format("%.4f", m);
            System.out.println(infileS + " with mean fst " + out + " is completed");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return out;
    }

    /**
     *
     */
    public void mkFstCommandbasedSNP() {

        //local file： one: group two: script
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/000_group/hexaandDi";

        // HPC file: group fileDirS
        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandTetra";
        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandDi";

        // HPC file: output fileDirS
        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/013_VMapIIbyRef";
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/002_fst_basedSNP";

        List<File> fs = IOUtils.getVisibleFileListInDir(groupHexaandTetraDirS);
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        List<File> fs2 = IOUtils.getVisibleFileListInDir(groupHexaandDiDirS);
        File[] group2FileS = fs2.toArray(new File[fs2.size()]);

        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。

        for (int i = 0; i < group1FileS.length - 1; i++) {
            String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
            for (int j = i + 1; j < group1FileS.length; j++) {
                String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                //vcftools --gzvcf test.vcf.gz --weir-fst-pop ../groups/Teosinte.txt --weir-fst-pop ../groups/Stiff_stalk.txt --weir-fst-pop ../groups/Non_stiff_stalk.txt --out out.txt
                String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outputDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group1FileDirS, group1FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group1FileDirS, group1FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S + " --weir-fst-pop " + group2S + " --out " + outfileS);
                }
            }
        }

        for (int i = 0; i < group2FileS.length - 1; i++) {
            String pop1 = group2FileS[i].getName().replace(".txt", "");
            for (int j = i + 1; j < group2FileS.length; j++) {
                String pop2 = group2FileS[j].getName().replace(".txt", "");
                String[] chrArr = {"1D", "2D", "3D", "4D", "5D", "6D", "7D"};
                for (int k = 0; k < chrArr.length; k++) {
                    String chr = chrArr[k];
                    String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                    String outfileS = new File(outputDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                    String group1S = new File(group2FileDirS, group2FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group2FileDirS, group2FileS[j].getName()).getAbsolutePath();
                    System.out.println("vcftools --vcf " + infileS + " --weir-fst-pop " + group1S + " --weir-fst-pop " + group2S + " --out " + outfileS);
                }
            }
        }
    }

}



