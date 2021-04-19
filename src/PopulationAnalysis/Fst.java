package PopulationAnalysis;

import AoUtils.*;
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

//        this.mergeFSTwindow();

        //********************************* VMap2.0 after -- new Version ********************//
//        this.mkFstCommandbasedwinndow2();
//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/001_script_based100kwindow_50kstep/fst_baesd100kwindow_50kstep_20200907.sh",21,1);
//        this.window();
//        this.addGroupToFstwindow();

//        this.getMeanFstValue();
//        this.getMatrixFst();

//        this.mkFstCMD_by2File();
//        this.mergeExonVCF();


//        this.changeChrPos();
//        this.filterMAF();
//        this.mkJAVAscript(); //script for filtering MAF of VCF
//        this.mkFstCMD_single();

    }

    /**
     * 2件事情：1.将chr pos 转换为参考基因组的1A形式，2.将1和2文件合并，并命名为1A
     */
    public void mergeExonVCF(){

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/007_tetraploid_Dm_DE/000_exonVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/007_tetraploid_Dm_DE/000_exonVCF_ref";
        new CountSites().mergeVCFfileandChangeChrPos_chr1and2(infileDirS,outfileDirS);
    }

    /**
     * 适用于只计算2个群体的 fst
     */
    public void mkFstCMD_by2File(){
        int window = 2000;
        int step = 2000;
        int numcmd = 1;
        //其他需要修改参数： 输入文件名称
//        String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();

        System.out.println("mkdir 001_srcipt_based" + window + "window_" + step + "step");
        System.out.println("mkdir 002_fst_based" + window + "window_" + step + "step");

        //local file： one: group
        String pop1FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_bySubspecies/Durum.txt";
        String pop2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_bySubspecies/Domesticated_emmer.txt";

        // HPC file: group fileDirS
        String groupFileDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/000_pop_bySubspecies";

        // HPC file: output fileDirS
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/007_tetraploid_Dm_DE/000_exonVCF_ref";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/007_tetraploid_Dm_DE/013_fstTest/001";
        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/007_tetraploid_Dm_DE/010_test3_pifst/001_script_based100kwindow_2kbstep/fst_based2kwindow_2kstep_20201105.sh"; //local script file
        new File(scriptS).getParentFile().mkdirs();

        List<File> fs = new ArrayList<>();
        fs.add(new File(pop1FileS)); fs.add(new File(pop2FileS));
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        try{
            BufferedWriter bw = AoFile.writeFile(scriptS);
            for (int i = 0; i < group1FileS.length - 1; i++) {
                String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
                for (int j = i + 1; j < group1FileS.length; j++) {
                    String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                    String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                    for (int k = 0; k < chrArr.length; k++) {
                        String chr = chrArr[k];
                        String infileS = new File(infileDirS, "chr" + chr + "_exon_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                        String group1S = new File(groupFileDirS, group1FileS[i].getName()).getAbsolutePath();
                        String group2S = new File(groupFileDirS, group1FileS[j].getName()).getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --weir-fst-pop ").append(group1S).append(" --weir-fst-pop ").append(group2S);
                        sb.append(" --fst-window-size ").append(window).append(" --fst-window-step ").append(step);
                        sb.append(" --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();

                    }
                }
            }

            bw.flush();bw.close();
            SplitScript.splitScript3(scriptS,numcmd); //脚本拆分

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }



    /**
     * 将分组文件转换成矩阵，并用来画PCA图,分成 AB D 两个文件
     */
    public void getMatrixFst(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/002_fst_based100000window_50000step/002/002_fst_mean.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/002_fst_based100000window_50000step/002/003_fst_matrix.txt";
        String outfileS2 = new File(outfileS).getAbsolutePath().replaceFirst(".txt","_Dsub.txt");
        //建立含有二维数组的List

        String[] subspecies = {"Wild_emmer","Domesticated_emmer","Free_threshing_tetraploid","Landrace","Cultivar"};
        String[] subspecies2 = {"Landrace","Cultivar","Ae.tauschii"};
        Arrays.sort(subspecies);Arrays.sort(subspecies2);
        double[][] aArray = new double[subspecies.length][subspecies.length];
        double[][] bArray = new double[subspecies.length][subspecies.length];
        double[][] dArray = new double[subspecies2.length][subspecies2.length];


        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            BufferedWriter bw2 = AoFile.writeFile(outfileS2);

            ///////////////////////////////////////////// A sub
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                if (!l.get(0).equals("A"))continue;
                cnt++;
                String pop1 = l.get(1).split("_VS_")[0];
                String pop2 = l.get(1).split("_VS_")[1];
                double fst = Double.parseDouble(l.get(2));
                int index1 = Arrays.binarySearch(subspecies,pop1);
                int index2 = Arrays.binarySearch(subspecies,pop2);
                aArray[index1][index2] = fst;
                aArray[index2][index1] = fst;
            }
            System.out.println(cnt + " chr A");
            br.close();

            bw.write("Sub\tGroup\t" + subspecies[0]);
            for (int i = 1; i < subspecies.length; i++) {
                bw.write("\t");
                bw.write(subspecies[i]);
            }
            bw.newLine();
            for (int i = 0; i < aArray.length; i++) {
                bw.write("A\t" + subspecies[i]);
                for (int j = 0; j < aArray[0].length; j++) {
                    String v = String.valueOf(aArray[i][j]);
                    bw.write("\t");
                    bw.write(v);
                }
                bw.newLine();
            }

            ///////////////////////////////////////////// B sub
            br = AoFile.readFile(infileS);
            header = br.readLine();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                if (!l.get(0).equals("B"))continue;
                cnt++;
                String pop1 = l.get(1).split("_VS_")[0];
                String pop2 = l.get(1).split("_VS_")[1];
                double fst = Double.parseDouble(l.get(2));
                int index1 = Arrays.binarySearch(subspecies,pop1);
                int index2 = Arrays.binarySearch(subspecies,pop2);
                bArray[index1][index2] = fst;
                bArray[index2][index1] = fst;

            }
            System.out.println(cnt + " chr B");
            br.close();

            for (int i = 0; i < bArray.length; i++) {
                bw.write("B\t" + subspecies[i]);
                for (int j = 0; j < bArray[0].length; j++) {
                    String v = String.valueOf(bArray[i][j]);
                    bw.write("\t");
                    bw.write(v);
                }
                bw.newLine();
            }

            ///////////////////////////////////////////// D sub
            br = AoFile.readFile(infileS);
            header = br.readLine();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                if (!l.get(0).equals("D"))continue;
                cnt++;
                String pop1 = l.get(1).split("_VS_")[0];
                String pop2 = l.get(1).split("_VS_")[1];
                double fst = Double.parseDouble(l.get(2));
                int index1 = Arrays.binarySearch(subspecies2,pop1);
                int index2 = Arrays.binarySearch(subspecies2,pop2);
                dArray[index1][index2] = fst;
                dArray[index2][index1] = fst;
            }
            System.out.println(cnt + " chr D");
            br.close();


            bw2.write("Sub\tGroup\t" + subspecies2[0]);
            for (int i = 1; i < subspecies2.length; i++) {
                bw2.write("\t");
                bw2.write(subspecies2[i]);
            }
            bw2.newLine();
            for (int i = 0; i < dArray.length; i++) {
                bw2.write("D\t" + subspecies2[i]);
                for (int j = 0; j < dArray[0].length; j++) {
                    String v = String.valueOf(dArray[i][j]);
                    bw2.write("\t");
                    bw2.write(v);
                }
                bw2.newLine();
            }
            bw.flush();bw2.flush();
            bw.close();bw2.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void getMeanFstValue(){
        String taxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        HashMap<String,String> hm = new AoFile().getHashMapStringKey(taxaList,15,3);//亚种和倍性之间的hashmap
//        String infileDirS = "";
//        String outfileS = "";


        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/002_fst_based100000window_50000step/001";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/002_fst_based100000window_50000step/002/001_Fst_bySubspecies_20200907.txt";

        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tSub\tGroup\tPloidy\tMEAN_FST");
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
     * 向window滑窗的结果添加group信息，区分不同分组的比较
     */
    public void addGroupToFstwindow(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/003_fst_windowbyJava";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/001_fst";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/002_fst/chr4A.windowed.weir.fst";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/002_fst_based100000window_50000step/003_windowbyJava";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/002_fst_based100000window_50000step/004_merge/001_Fst_2Mwindow_1Mstep.txt.gz";

//        String infileDirS = "";
//        String outfileDirS = "";

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

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/001_fst";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/004_mix_4A/003_fst_windowbyJava";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/002_fst_based100000window_50000step/001";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/002_fst_based100000window_50000step/003_windowbyJava";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/003_fst_based100000window_50000step_testforMAF/001";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/003_fst_based100000window_50000step_testforMAF/003_windowbyJava";

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            AoFile.readheader(infileS);
            int chrColumn = 0;
            int posIndex = 1;
            int valueIndex = 4;
            double window = 2000000;
            double step = 1000000;
            String name = new File(infileS).getName().split(".fst")[0] + "_" + window + "window_" + step + "step.txt.gz";
            String outfileS = new File(outfileDirS,name).getAbsolutePath();

            new AoWinScan().getwindowDistrbution_general(infileS,chrColumn,posIndex,valueIndex,window,step,outfileS);
        });

    }


    public void mkJAVAscript(){
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/005_changeChrPos_from002";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/005_changeChrPos_from002/filterMAF0.1";
        String[] input = {"Asubgenome_RefChrPos.vcf.gz","Bsubgenome_RefChrPos.vcf.gz","Dsubgenome_RefChrPos.vcf.gz"};
        String ratio = "0.1";
        for (int i = 0; i < input.length; i++) {
            String name = input[i];
            String infileS = new File(infileDirS,name).getAbsolutePath();
            String outfileS = new File(outfileDirS, name.replaceFirst(".vcf.gz","_MAF0.1.vcf.gz")).getAbsolutePath();
            System.out.println("java -jar 053_filterVCFbyMAF.jar " + infileS + " " + ratio + " "+ outfileS + " > log_053_filterVCFbyMAF_2021-04-11.txt 2>&1 &" );
        }

    }

    public void filterMAF(){
        String infileS = "/Users/Aoyue/Documents/testout/Dsubgenome_diploid_RefChrPos.vcf.gz";
        String outfileS ="/Users/Aoyue/Documents/testout/Dsubgenome_diploid_RefChrPos_maf0.1.vcf.gz";
        double ratio = 0.1;
        CalVCF.filterMAFinVCF(infileS,ratio,outfileS);
    }

    /**
     * 将抽样的VCF文件的位置转换成RefChrPos
     */
    public void changeChrPos(){

        String infileDirS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/002_mergeVCFtoSub";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/005_changeChrPos_from002";
//        String infileDirS = "/Users/Aoyue/Documents/test";
//        String outfileDirS = "/Users/Aoyue/Documents/testout";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().split(".vcf.gz")[0] + "_RefChrPos.vcf.gz").getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                List<String> l = new ArrayList<>();
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    if (temp.startsWith("#")){
                        bw.write(temp);
                        bw.newLine();
                    }else{
                        int chr = Integer.parseInt(l.get(0));
                        int pos = Integer.parseInt(l.get(1));
                        String RefChr = RefV1Utils.getChromosome(chr,pos);
                        int RefPos = RefV1Utils.getPosOnChromosome(chr, pos);
                        sb.setLength(0);
                        sb.append(RefChr).append("\t").append(RefPos);
                        for (int i = 2; i < l.size(); i++) {
                            sb.append("\t").append(l.get(i));
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

        // java -jar GeneticLoad.jar > log_changeChrPos_2021-04-10.txt 2>&1 &
    }


    /**
     * 输出运行Fst的单行命令
     * vcftools --gzvcf /data4/home/aoyue/vmap2/genotype/mergedVCF/105_VMap2.1ByRef/chr7D_vmap2.1.vcf.gz --weir-fst-pop /data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandDi/Cultivar.txt --weir-fst-pop /data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandDi/Landrace.txt --fst-window-size 100000 --fst-window-step 50000 --out /data4/home/aoyue/vmap2/analysis/031_popGen/002_fst_based100000window_50000step/001/Cultivar_VS_Landrace_chr7D
     *
     * 当VCF文件只有一个时，多个组的循环
     */
    public void mkFstCMD_single(){

        // ******************************************* 参数设置 *********************************************************************
// 先做一个测试，MAF过滤前后，FST是否有变化
        //        //local file： one: group
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/000_group/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/000_group/hexaandDi";
        // HPC file: group fileDirS
        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/000_group/hexaandTetra";
        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/000_group/hexaandDi";

//        String infileDirS = ""; //HPC
//        String outfileDirS = ""; //HPC
//        String infileS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/002_mergeVCFtoSub/Asubgenome.vcf.gz"; //未转换成RefChrPos
//        String infileS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/005_changeChrPos_from002/beforeFilterMAF0.1/Dsubgenome_RefChrPos.vcf.gz";
        String infileS = "/data4/home/aoyue/vmap2/analysis/025_subsetVCF/005_changeChrPos_from002/filterMAF0.1/Bsubgenome_RefChrPos_MAF0.1.vcf.gz";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/005_testFst_MAFfilter/002_fst_based100000window_50000step";
        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_Fst/001_script_based100kwindow_50kstep/fst_based100kwindow_50kstep_2021-04-10.sh";

        //        //参数
        int window = 100000;
        int step = 50000;
        int numcmd = 1;
        //其他需要修改参数： 输入文件名称
//        String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
        System.out.println("mkdir 001_script_based" + window + "window_" + step + "step");
        System.out.println("mkdir 002_fst_based" + window + "window_" + step + "step");


        List<File> fs = AoFile.getFileListInDir(groupHexaandTetraDirS);
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        List<File> fs2 = AoFile.getFileListInDir(groupHexaandDiDirS);
        File[] group2FileS = fs2.toArray(new File[fs2.size()]);


        try{
            BufferedWriter bw = AoFile.writeFile(scriptS);
            for (int i = 0; i < group1FileS.length - 1; i++) {
                String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
                for (int j = i + 1; j < group1FileS.length; j++) {
                    String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_" + new File(infileS).getName().replaceFirst(".vcf.gz","")).getAbsolutePath();
                    String group1S = new File(group1FileDirS, group1FileS[i].getName()).getAbsolutePath();
                    String group2S = new File(group1FileDirS, group1FileS[j].getName()).getAbsolutePath();
                    StringBuilder sb = new StringBuilder();
                    sb.append("vcftools --gzvcf ").append(infileS).append(" --weir-fst-pop ").append(group1S).append(" --weir-fst-pop ").append(group2S);
                    sb.append(" --fst-window-size ").append(window).append(" --fst-window-step ").append(step);
                    sb.append(" --out ").append(outfileS);
                    System.out.println(sb.toString());
                    bw.write(sb.toString());bw.newLine();
                }
            }

//            for (int i = 0; i < group2FileS.length - 1; i++) {
//                String pop1 = group2FileS[i].getName().replace(".txt", "");
//                for (int j = i + 1; j < group2FileS.length; j++) {
//                    String pop2 = group2FileS[j].getName().replace(".txt", "");
//                    String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_" + new File(infileS).getName().replaceFirst(".vcf.gz","")).getAbsolutePath();
//                    String group1S = new File(group2FileDirS, group2FileS[i].getName()).getAbsolutePath();
//                    String group2S = new File(group2FileDirS, group2FileS[j].getName()).getAbsolutePath();
//                    StringBuilder sb = new StringBuilder();
//                    sb.append("vcftools --gzvcf ").append(infileS).append(" --weir-fst-pop ").append(group1S).append(" --weir-fst-pop ").append(group2S);
//                    sb.append(" --fst-window-size ").append(window).append(" --fst-window-step ").append(step);
//                    sb.append(" --out ").append(outfileS);
//                    System.out.println(sb.toString());
//                    bw.write(sb.toString());bw.newLine();
//
//                }
//            }

            bw.flush();bw.close();
            SplitScript.splitScript3(scriptS,numcmd); //脚本拆分

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 终极版本！运行fst命令的脚本
     */
    public void mkFstCommandbasedwinndow2() {
        /**
         * input: including local path and HPC path
         */
        // 分组文件
        //String groupHexaandTetraDirS = ""; //本地
        //String groupHexaandDiDirS = ""; //本地
        //String group1FileDirS = ""; //HPC
//        String group2FileDirS = ""; //HPC
        //输入输出文件
//        String infileDirS = ""; //HPC
//        String outfileDirS = ""; //HPC
//          String scriptS = ""; //local script file
//        //参数
//        int window = 100000;
//        int step = 50000;
//        int numcmd = 8;

//        //其他需要修改参数： 输入文件名称
////        String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
//
//        System.out.println("mkdir 001_srcipt_based" + window + "window_" + step + "step");
//        System.out.println("mkdir 002_fst_based" + window + "window_" + step + "step");
//
//

// ******************************************* 参数设置 *********************************************************************
//        //local file： one: group
//        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/000_group/hexaandTetra";
//        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/000_group/hexaandDi";
//
//        // HPC file: group fileDirS
//        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandTetra";
//        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandDi";
//
//        // HPC file: output fileDirS
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/105_VMap2.1ByRef";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/002_fst_based100000window_50000step/001";
//        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_script_based100kwindow_50kstep/fst_based100kwindow_50kstep_20200213.sh"; //local script file
//        new File(scriptS).getParentFile().mkdirs();

//        int window = 2000;
//        int step = 2000;
//        int numcmd = 21;

//        //其他需要修改参数： 输入文件名称
////        String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
//        System.out.println("mkdir 001_script_based" + window + "window_" + step + "step");
//        System.out.println("mkdir 002_fst_based" + window + "window_" + step + "step");

// ******************************************* 参数设置 *********************************************************************
//        //local file： one: group
//        String groupHexaandTetraDirS = "";
//        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/000_group/hexaandDi";
//
//        // HPC file: group fileDirS
//        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandTetra";
//        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/021_popGen/101_Fst/000_group/hexaandDi";
//
//        // HPC file: output fileDirS
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/105_VMap2.1ByRef";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/002_fst_based100000window_50000step/001";
//        String scriptS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/001_script_based100kwindow_50kstep/fst_based100kwindow_50kstep_20200213.sh"; //local script file
//        new File(scriptS).getParentFile().mkdirs();

// ******************************************* 参数设置 *********************************************************************
// 先做一个测试，MAF过滤前后，FST是否有变化
        //        //local file： one: group
        String groupHexaandTetraDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/000_group/hexaandTetra";
        String groupHexaandDiDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/039_popGen/000_group/hexaandDi";
        // HPC file: group fileDirS
        String group1FileDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/000_group/hexaandTetra";
        String group2FileDirS = "/data4/home/aoyue/vmap2/analysis/031_popGen/000_group/hexaandDi";

        String infileDirS = ""; //HPC
        String outfileDirS = ""; //HPC
        String scriptS = ""; //local script file

        //        //参数
        int window = 100000;
        int step = 50000;
        int numcmd = 8;
        //其他需要修改参数： 输入文件名称
//        String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
        System.out.println("mkdir 001_srcipt_based" + window + "window_" + step + "step");
        System.out.println("mkdir 002_fst_based" + window + "window_" + step + "step");


        List<File> fs = AoFile.getFileListInDir(groupHexaandTetraDirS);
        File[] group1FileS = fs.toArray(new File[fs.size()]);

        List<File> fs2 = AoFile.getFileListInDir(groupHexaandDiDirS);
        File[] group2FileS = fs2.toArray(new File[fs2.size()]);


        try{
            BufferedWriter bw = AoFile.writeFile(scriptS);
            for (int i = 0; i < group1FileS.length - 1; i++) {
                String pop1 = group1FileS[i].getName().replace(".txt", ""); //第一组的名字
                for (int j = i + 1; j < group1FileS.length; j++) {
                    String pop2 = group1FileS[j].getName().replace(".txt", ""); //第二组的名字
                    String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};
                    for (int k = 0; k < chrArr.length; k++) {
                        String chr = chrArr[k];
                        String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                        String group1S = new File(group1FileDirS, group1FileS[i].getName()).getAbsolutePath();
                        String group2S = new File(group1FileDirS, group1FileS[j].getName()).getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --weir-fst-pop ").append(group1S).append(" --weir-fst-pop ").append(group2S);
                        sb.append(" --fst-window-size ").append(window).append(" --fst-window-step ").append(step);
                        sb.append(" --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();

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
                        String infileS = new File(infileDirS, "chr" + chr + "_vmap2.1.vcf.gz").getAbsolutePath();
                        String outfileS = new File(outfileDirS, pop1 + "_VS_" + pop2 + "_chr" + chr).getAbsolutePath();
                        String group1S = new File(group2FileDirS, group2FileS[i].getName()).getAbsolutePath();
                        String group2S = new File(group2FileDirS, group2FileS[j].getName()).getAbsolutePath();
                        StringBuilder sb = new StringBuilder();
                        sb.append("vcftools --gzvcf ").append(infileS).append(" --weir-fst-pop ").append(group1S).append(" --weir-fst-pop ").append(group2S);
                        sb.append(" --fst-window-size ").append(window).append(" --fst-window-step ").append(step);
                        sb.append(" --out ").append(outfileS);
                        System.out.println(sb.toString());
                        bw.write(sb.toString());bw.newLine();
                    }
                }
            }

            bw.flush();bw.close();
            SplitScript.splitScript3(scriptS,numcmd); //脚本拆分

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }




    /**
     * 将fst window scan 计算的FST,多个文件合并起来成为一个文件，并添加一列分组Group信息
     * 暂时只将AE.tauschii合并起来吧
     *
     */
    public void mergeFSTwindow(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/004_fst_based2Mwindow_1Mstep/001";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/004_fst_based2Mwindow_1Mstep/002_merge/FST_Ae.tauschii_2Mwindow_1Mstep_20200227.txt";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/008_fst_based2Mwindow_1Mstep/001";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/101_Fst/008_fst_based2Mwindow_1Mstep/002_merge/FST_Ae.tauschii_subLandrace_2Mwindow_1Mstep_20200227.txt";

        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //read header
            bw.write(br.readLine() + "\tGroup");
            bw.newLine();

            int cnttotal = 0;
            //read context
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                String name = fs[i].getName();
                if (!name.contains("Ae.tauschii"))continue;
                String group = name.substring(0,name.indexOf("_chr"));
                br = AoFile.readFile(infileS);
                br.readLine();
                String temp = null; //read header
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
            BufferedReader br = AoFile.readFile(infileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            TDoubleArrayList vList = new TDoubleArrayList();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
//                String t = l.get(4); //********************** 需要修改
                String t = l.get(4);
                if (t.startsWith("-n") || t.startsWith("n") || t.startsWith("-") ) { //pos中 fst值为 -nan 的全部去掉
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



