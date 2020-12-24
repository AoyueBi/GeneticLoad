/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.SplitScript;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class SIFT {

    public SIFT() {
        //this.addGeneBiotype();
        //new FastaBit("/Users/Aoyue/Documents/chr036.fa.gz").getName(0);
        //this.parseGff3();
        //this.mkInputDirs();
        //this.mkParameterConfigChr1_42();
        //this.mkJavaCmdchr1_42_SIFTdb();
        //this.mvDatabase();
        //this.annotatorVCF("a","b");
        //new CountSites().extractHapPosAllele("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/test/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/test2/");
//        this.generateSIFT();
        //this.calVariantsType("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/001_DB_addSIFT/abd/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/002_calVariantsType/abd_variantsType.txt");
        //this.calVariantsType("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/001_DB_addSIFT/ab/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/002_calVariantsType/ab_variantsType.txt");
        //this.calVariantsType("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/001_DB_addSIFT/d/", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/002_calVariantsType/d_variantsType.txt");
        //this.calVariantsType("/data4/home/aoyue/vmap2/analysis/015_annoDB/003_addSIFT/", "/data4/home/aoyue/vmap2/analysis/015_annoDB/004_calSIFTNum/variantsType.txt");
        /**
         * ******** 正式数据的计算**************
         */
//        this.annotatorVCF2("/data4/home/aoyue/vmap2/genotype/mergedVCF/002_biMAF0.005VCF/", "/data4/home/aoyue/vmap2/analysis/008_sift/result/001_annotatorResult");
        //new SplitScript().splitBwaScript("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/004_newScript/sh_vcfAnnotator.sh", "vcfAnnotator", 7, 6);

        //this.calVariantsType2("/data4/home/aoyue/vmap2/analysis/015_annoDB/003_addSIFT/", "/data4/home/aoyue/vmap2/analysis/015_annoDB/004_calSIFTNum/cal_sift20190930.txt");
        //new CountSites().mergeChr1and2txt_double("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/005_result2/cal_sift20190930.txt", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/005_result2/cal_sift.byChrRef.txt");
        /**
         * ******** 正式数据的计算**************
         */
//        this.annotatorVCF2("/data4/home/aoyue/vmap2/genotype/mergedVCF/008_maf0.01SNPbyPop/", "/data4/home/aoyue/vmap2/analysis/008_sift/001_result/001_annotatorResult/");
//        new SplitScript().splitBwaScript("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/006_newScript_MAF0.01byPop/sh_vcfAnnotator20191029.sh", "vcfAnnotator", 7, 6);
//        this.mkdis();
//        this.ifDone("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/007_log_siftAnnotator", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/008_log_summary/log_summary.txt");
//this.move();

        /**
         * ******** for new test when rever the ref and alt **************
         *
         */
//        this.reverseRefAltallelebyExonVCF();
//        this.annotatorVCF2("/data4/home/aoyue/vmap2/analysis/008_sift/001_result/002_annotatorResult/001_reversedRefAltVCF","/data4/home/aoyue/vmap2/analysis/008_sift/001_result/002_annotatorResult");
//        new SplitScript().splitScript("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/010_exonVCF_reverseRefAlt/002_script/sh_vcfAnnotator20200103.sh","vcfAnnotator",7,6);
//        this.mkdis();
//        this.move();

//        this.annotatorVCF2("/data4/home/aoyue/vmap2/analysis/008_sift/002_result_exonVCF/001_annotatorResult/001_exonVCF","/data4/home/aoyue/vmap2/analysis/008_sift/002_result_exonVCF/001_annotatorResult");
//        new SplitScript().splitScript("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/010_exonVCF_reverseRefAlt/002_script/sh_vcfAnnotator20200105.sh","vcfAnnotator",7,6);

//        this.annotatorVCF2("/data4/home/aoyue/vmap2/analysis/008_sift/002_result_exonVCF/002_annotatorResult/001_reversedRefAltVCF","/data4/home/aoyue/vmap2/analysis/008_sift/002_result_exonVCF/002_annotatorResult");
//        new SplitScript().splitScript("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/010_exonVCF_reverseRefAlt/002_script/sh_vcfAnnotator20200105_2.sh","vcfAnnotator",7,6);

        //******** 2020 版本 VMapII 的计算 **************//
//        this.mkdis();
//        this.annotatorVCF3();
//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/034_sift/002_script/Annotator_20200719.sh",6,7);
//        this.move();
//        this.reverseRefAltallelebyExonVCF();

    }

//     find -name "*.vcf" | cut -c3- ; 本地获取运行脚本，输出在netbeans的output界面。
//     -c	To run on command line
//     -i	Path to your input variants file in VCF format
//     -d	Path to SIFT database directory
//     -r       Path to your output results folder
//     -t	To extract annotations for multiple transcripts (Optional)
//    输出的文件产生的日志也保存下来，便于计算时间。 注意输入文件必须是解压后的vcf

    public void annotatorVCF3(){ //一次性将所有的jar都运行上
//        String infileDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/002_exonSNPVCF";
//        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/008_sift/003_result_Vmap2.1-2020_exonVCF";

        String infileDirS = "/data4/home/aoyue/vmap2/analysis/028_sift/001_exon_refref/001_exonVCF_refref";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/028_sift/001_exon_refref";

        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            String infileS = new File(infileDirS,"chr" + chr + "_exon_vmap2.1.vcf").getAbsolutePath();
            String infileS = new File(infileDirS,"chr" + chr + "_exon_vmap2.1_RefRef.vcf").getAbsolutePath();

            String outfileDirS2 = new File(outfileDirS,"output" + chr).getAbsolutePath();
            File f = new File(infileS);
            String name = f.getName().split(".vcf")[0];
            String logFileS = new File(outfileDirS,"log_" + name + "_annotator.txt").getAbsolutePath();

            System.out.println("java -jar /data1/home/aoyue/biosoftware/SIFT4g_Annotator/SIFT4G_Annotator_v2.4.jar -c -i "
                    + infileS
                    + " -d /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/siftDatabase -r "
                    + outfileDirS2
                    + " -t > "
                    + logFileS
                    + "");
            //java -jar /data1/home/aoyue/biosoftware/SIFT4g_Annotator/SIFT4G_Annotator_v2.4.jar -c -i /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input001/dbSNP/chr036.Dgenome.bi.vcf.gz -d /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input001/abd_iwgscV1 -r /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/output/output036 -t &

        }
    }





    /**
     * 为了检查ref allele 和 alt allele 互换顺序后，进行SIFT计算，有无其他结果。实践证明，无结果。
     * 正确的思路是，Ref保持不变，将Alt替换成ref，只对ref进行sift的打分。
     *
     */
    public void reverseRefAltallelebyExonVCF(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/003_exonSNPVCF";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/010_exonVCF_reverseRefAlt/001_reverseRefAltallelebyExonVCF";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/010_exonVCF_reverseRefAlt/003_convertAltAlleletoRefAllelebyExonVCF";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/034_sift/001_exonVCF_RefRef";


        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().split(".vcf")[0] + "_RefRef.vcf").getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                int cnttotal = 0;
                int cntsubset = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        cnttotal++;
                        List<String> l = PStringUtils.fastSplit(temp);
                        if (l.get(3).contains(",")) {
                            continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                        }
                        // 写出的时候 3和4 互换位置
                        String ref = l.get(3);
                        String alt = l.get(4);

//                        bw.write(l.get(0)+ "\t"+l.get(1)+ "\t"+l.get(2)+ "\t"+l.get(4)+ "\t"+l.get(3));  //alt 和 ref 进行互换，然后用sift进行注释,此法不成立，被抛弃！！！
                        bw.write(l.get(0)+ "\t"+l.get(1)+ "\t"+l.get(2)+ "\t"+l.get(3)+ "\t"+l.get(3)); //只将alt进行替换，提成ref的碱基，然后用sift进行注释，从而判断ref是否是有害的。
                        for (int i = 5; i < l.size(); i++) {
                            bw.write("\t" + l.get(i));
                        }
                        bw.newLine();
                        cntsubset++;
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\twith " + cnttotal + " bp has a subset of\t" + cntsubset + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     * 将xls结果移动到同一个文件夹中，以备后续需要
     */
    public void move() {
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
//            System.out.println("mv output" + chr + "/chr" + chr + ".subgenome.maf0.01byPop.SNP_SIFTannotations.xls output/");
//            System.out.println("mv output" + chr + "/chr" + chr + "_exon_vmap2.1_reverseRefAlt_SIFTannotations.xls output/");
//            System.out.println("mv output" + chr + "/chr" + chr + "_exon_vmap2.1_SIFTannotations.xls output/");
//            System.out.println("mv output" + chr + "/chr" + chr + "_exon_vmap2.1_SIFTannotations.xls output/");
            System.out.println("mv output" + chr + "/chr" + chr + "_exon_vmap2.1_RefRef_SIFTannotations.xls output/");


        }

    }

    /**
     * 目的：检查/008_sift/007_log_siftAnnotator中所有log文件，查看是否全部运行完毕。车看第12行是否都包括completed,把第9行的内容提取出来，把第8行的表头提取一次
     *
     * @param infileDirS
     * @param outfileS
     */
    public void ifDone(String infileDirS, String outfileS) {
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        //System.out.println("Chr\tSNP_Num");

        try {
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }

            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 1; i < 15; i++) {
                String temp = br.readLine();
                if (i == 8) {
                    bw.write(temp);
                    bw.newLine();
                }
            }

            //读正文部分
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                String temp = null;
                //int chrint = Integer.parseInt(fs[i].getName().substring(3, 6));
                int cnt = 0;
                for (int j = 1; j < 15; j++) {  //每个log文件有14行，我们只要第9行的数据
                    temp = br.readLine();
                    if (j == 9) {
                        bw.write(temp);
                        bw.newLine();
                    }
                }
                System.out.println(fs[i].getName() + " is completed at" + outfileS);
            }
            br.close();
            bw.flush();
            bw.close();

            //下面一段程序，是用来去除文本中含有的tab键空格
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 建立42个空的文件夹，上传至集群
     */
    public void mkdis() {
        System.out.println("mkdir output");
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            System.out.println("mkdir output" + chr);
        }

    }

    /**
     * 本程序对merge后的vcf注释结果进行分类计算
     *
     * @param infileDirS
     * @param outfileS
     */
    public void calVariantsType2(String infileDirS, String outfileS) {

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);

        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chromosome\tNonSynonymousDeleteriousMutation\tNonSynonymousTolerentMutation\tSynonymousMutation\tNoncodingMutation\t"
                    + "StartLostDeleteriousMutation\tStartLostTolerentMutation\tStopGainMutation\tStopLossMutation\t"
                    + "FrameshiftInsertionMutation\n");
            for (int i = 0; i < fs.length; i++) {
                System.out.println("Now " + fs[i].getName() + " are being processed");
                String infileS = fs[i].getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }

                String temp = br.readLine();
                int cnt = 0;
                int NonSynonymousDeleteriousMutation = 0;
                int NonSynonymousTolerentMutation = 0;
                int SynonymousMutation = 0;
                int NoncodingMutation = 0;
                int StartLostMutation = 0;
                int StartLostDeleteriousMutation = 0;
                int StartLostTolerentMutation = 0;
                int StopGainMutation = 0;
                int StopLossMutation = 0;
                int FrameshiftInsertionMutation = 0;
                int naCount = 0; //
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    l = PStringUtils.fastSplit(temp);
                    String Variant_type = l.get(9);
                    String score = l.get(10);
                    if (Variant_type.equals("NONSYNONYMOUS")) { //如果score的值为NA，那么NONSYNONYMOUS的预测也是不确定的
                        if (score.equals("NA")) {
                            naCount++; //即是NONSYNONYMOUS但是没有sift值的个数
                        } else {
                            if (Double.valueOf(score) < 0.05) {
                                NonSynonymousDeleteriousMutation++;
                            } else {
                                NonSynonymousTolerentMutation++;
                            }

                        }

                    } else if (Variant_type.equals("SYNONYMOUS")) {
                        SynonymousMutation++;
                    } else if (Variant_type.equals("START-LOST")) {
                        StartLostMutation++;
                        if (score.equals("NA")) {
                            //StartLostTolerentMutation++;
                        } else {
                            if (Double.valueOf(score) < 0.05) {
                                StartLostDeleteriousMutation++;
                            } else {
                                StartLostTolerentMutation++;
                            }

                        }

                    } else if (Variant_type.equals("STOP-GAIN")) {
                        StopGainMutation++;
                    } else if (Variant_type.equals("STOP-LOSS")) {
                        StopLossMutation++;
                    } else if (Variant_type.equals("NONCODING")) {
                        NoncodingMutation++;
                    } else if (Variant_type.equals("FRAMESHIFT INSERTION")) {
                        FrameshiftInsertionMutation++;
                    }

                }
                br.close();

                bw.write(Integer.parseInt(fs[i].getName().substring(3, 6)) + "\t" + NonSynonymousDeleteriousMutation + "\t" + NonSynonymousTolerentMutation
                        + "\t" + SynonymousMutation + "\t" + NoncodingMutation + "\t" + StartLostDeleteriousMutation + "\t" + StartLostTolerentMutation
                        + "\t" + StopGainMutation + "\t" + StopLossMutation + "\t" + FrameshiftInsertionMutation + "\n");
            }
            bw.flush();
            bw.close();
            System.out.println("Calculating SIFT result is completed at " + outfileS);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

//     find -name "*.vcf" | cut -c3- ; 本地获取运行脚本，输出在netbeans的output界面。 
//     -c	To run on command line 
//     -i	Path to your input variants file in VCF format 
//     -d	Path to SIFT database directory 
//     -r       Path to your output results folder 
//     -t	To extract annotations for multiple transcripts (Optional)
//    输出的文件产生的日志也保存下来，便于计算时间。 注意输入文件必须是解压后的vcf
    /**
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void annotatorVCF2(String infileDirS, String outfileDirS) {
        //String infileDirS = "/data4/home/aoyue/vmap2/analysis/008_sift/dbSNP20190820/ab/";
        //String outfileDirS = "/data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/output/ab/";

        try {
//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_chrList/list_chrlineage.txt";
//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_chrList/list_chrMerged.txt";
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_chrList/list_chrMergedReverse.txt";
//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_chrList/list_Exon.txt";

            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                String chr = temp.substring(3, 6);
                System.out.println("java -jar /data1/home/aoyue/biosoftware/SIFT4g_Annotator/SIFT4G_Annotator_v2.4.jar -c -i "
                        + new File(infileDirS, temp).getAbsolutePath()
                        + " -d /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/siftDatabase -r "
                        + new File(outfileDirS, "output" + chr).getAbsolutePath()
                        + " -t > "
                        + new File(outfileDirS, "log_" + temp + ".annotator.txt").getAbsolutePath()
                        + "");
                //java -jar /data1/home/aoyue/biosoftware/SIFT4g_Annotator/SIFT4G_Annotator_v2.4.jar -c -i /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input001/dbSNP/chr036.Dgenome.bi.vcf.gz -d /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input001/abd_iwgscV1 -r /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/output/output036 -t &
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 目的：输出表格，计算出9类变异的数目。 NonSynonymousDeleteriousMutation
     * NonSynonymousTolerentMutation SynonymousMutation	NoncodingMutation
     * StartLostDeleteriousMutation StartLostTolerentMutation
     * FrameshiftInsertionMutation	StopGainMutation StopLossMutation
     *
     * @param infileDirS
     * @param outfileS
     */
    public void calVariantsType(String infileDirS, String outfileS) {

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);

        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chromosome\tNonSynonymousDeleteriousMutation\tNonSynonymousTolerentMutation\tSynonymousMutation\t"
                    + "StartLostDeleteriousMutation\tStartLostTolerentMutation\tStopGainMutation\tStopLossMutation\t"
                    + "NoncodingMutation\tFrameshiftInsertionMutation\n");
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }

                String temp = br.readLine();
                int cnt = 0;
                int NonSynonymousDeleteriousMutation = 0;
                int NonSynonymousTolerentMutation = 0;
                int SynonymousMutation = 0;
                int NoncodingMutation = 0;
                int StartLostMutation = 0;
                int StartLostDeleteriousMutation = 0;
                int StartLostTolerentMutation = 0;
                int StopGainMutation = 0;
                int StopLossMutation = 0;
                int FrameshiftInsertionMutation = 0;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    l = PStringUtils.fastSplit(temp);
                    String Variant_type = l.get(4);
                    String score = l.get(5);
                    if (Variant_type.equals("NONSYNONYMOUS")) { //如果score的值为NA，那么NONSYNONYMOUS的预测也是不确定的
                        if (score.equals("NA")) {
                            //NonSynonymousTolerentMutation++;
                        } else {
                            if (Double.valueOf(score) <= 0.05) {
                                NonSynonymousDeleteriousMutation++;
                            } else {
                                NonSynonymousTolerentMutation++;
                            }

                        }

                    } else if (Variant_type.equals("SYNONYMOUS")) {
                        SynonymousMutation++;
                    } else if (Variant_type.equals("START-LOST")) {
                        StartLostMutation++;
                        if (score.equals("NA")) {
                            StartLostTolerentMutation++;
                        } else {
                            if (Double.valueOf(score) <= 0.05) {
                                StartLostDeleteriousMutation++;
                            } else {
                                StartLostTolerentMutation++;
                            }

                        }

                    } else if (Variant_type.equals("STOP-GAIN")) {
                        StopGainMutation++;
                    } else if (Variant_type.equals("STOP-LOSS")) {
                        StopLossMutation++;
                    } else if (Variant_type.equals("NONCODING")) {
                        NoncodingMutation++;
                    } else if (Variant_type.equals("FRAMESHIFT INSERTION")) {
                        FrameshiftInsertionMutation++;
                    }

                }
                br.close();

                bw.write(Integer.parseInt(fs[i].getName().substring(3, 6)) + "\t" + NonSynonymousDeleteriousMutation + "\t" + NonSynonymousTolerentMutation
                        + "\t" + SynonymousMutation + "\t" + StartLostDeleteriousMutation + "\t" + StartLostTolerentMutation
                        + "\t" + StopGainMutation + "\t" + StopLossMutation + "\t" + NoncodingMutation + "\t" + FrameshiftInsertionMutation + "\n");
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * You need to modify 5 places:
     *
     * @num1: String siftDirS
     * @num2: String dbDirS
     * @num3: String outputDirS
     * @num4: String dbFileS
     * @num5: String outfileS
     */
    public void generateSIFT() {
        //分别建立染色体和对应亚基因组的关系
        int[] arra = {1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32, 37, 38};
        int[] arrb = {3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40};
        int[] arrd = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        HashMap<Integer, String> hml = new HashMap<>();
        Arrays.sort(arra);
        Arrays.sort(arrb);
        Arrays.sort(arrd);
        for (int i = 0; i < arra.length; i++) {
            hml.put(arra[i], "A");
            hml.put(arrb[i], "B");
            hml.put(arrd[i], "D");
        }
        /*==================================== 测试用 =============================================*/
//        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
//        //包含sift annotator注释结果的xls文件，可以用java文件输入读取， 非压缩格式
//        String siftDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/annotatorResult/";
//        //数据库的位点信息，包含 Chr	Pos	Ref	Alt    6	291	G	A，压缩格式
//        String dbDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/annoDB/";
//        //文件输出，在数据库添加SIFT信息  压缩格式
//        String outputDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/001_annoDB_addSIFT/";

//        /*==================================== abd =============================================*/
//        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
//        String siftDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/000_annotatorResult/ab/";
//        String dbDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/000_DB/ab/";
//        String outputDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/003_result/001_DB_addSIFT/ab/";
//                /*==================================== merge =============================================*/
//        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
//        String siftDirS = "/Users/Aoyue/Documents/siftOut";
//        String dbDirS = "/Users/Aoyue/Documents/annoDB";
//        String outputDirS = "/Users/Aoyue/Documents/annoDB2";
        /*==================================== HPC =============================================*/
//        String infileS = "/data1/publicData/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
//        String siftDirS = "/data4/home/aoyue/vmap2/analysis/008_sift/result/001_annotatorResult/output/";
//        String dbDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/001_step1/";
//        String outputDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/003_addSIFT/";

        String infileS = "/data1/publicData/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String siftDirS = "/data4/home/aoyue/vmap2/analysis/008_sift/001_result/001_annotatorResult/output";
        String dbDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/011_step2";
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/015_annoDB/012_addSIFT";
        //首先建立输出目录文件夹
        new File(outputDirS).mkdirs();
        //先处理pgf文件
        GeneFeature gf = new GeneFeature(infileS);
        int geneNum = gf.getGeneNumber(); //1.获取基因数
        String[] trans = new String[geneNum]; //2.定义trans数组,一个基因一个最长转录组
        for (int i = 0; i < geneNum; i++) { //循环基因
            int index = gf.getLongestTranscriptIndex(i); //获取最长转录本的索引
            trans[i] = gf.getTranscriptName(i, index); //根据索引，获取对应基因的最长转录本的名字
        }
        Arrays.sort(trans); //对转录本进行排序

        File[] fs = new File(siftDirS).listFiles(); //软件计算出的包含SIFT结果的xls文件和VCF文件；
        fs = IOUtils.listFilesEndsWith(fs, "xls");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> { //循环的是xls文件
            int chrIndex = Integer.parseInt(f.getName().substring(3, 6)) - 1; //获取文件的染色体号
            //chr006.Dgenome.bi.pos.txt.gz
            //chr001.ABDgenome.filterMiss.posAllele.txt.gz
            /*==================================== 此处需要修改 =============================================*/
            // dbFileS 即数据库的绝对路径
//            String dbFileS = new File(dbDirS, "chr" + PStringUtils.getNDigitNumber(3, chrIndex + 1) + ".ABgenome.filterMiss.posAllele.txt.gz").getAbsolutePath();
//            String outfileS = new File(outputDirS, "chr" + PStringUtils.getNDigitNumber(3, chrIndex + 1) + ".ABgenome.filterMiss.posAllele.SIFT.txt.gz").getAbsolutePath();
            //chr036.Dlineage.maf0.005.bi.AnnoDB.txt.gz
            //chr001_vmap2.1_AnnoDB_addDAF.txt.gz
//            String dbFileS = new File(dbDirS, "chr" + PStringUtils.getNDigitNumber(3, chrIndex + 1) + "." + hml.get(chrIndex + 1) + "lineage.maf0.005.bi.AnnoDB.txt.gz").getAbsolutePath();
//            String outfileS = new File(outputDirS, "chr" + PStringUtils.getNDigitNumber(3, chrIndex + 1) + ".lineage.maf0.005.bi.AnnoDB.addSIFT.txt.gz").getAbsolutePath();
            String dbFileS = new File(dbDirS, "chr" + PStringUtils.getNDigitNumber(3, chrIndex + 1) + "_vmap2.1_AnnoDB_addDAF.txt.gz").getAbsolutePath();
            String outfileS = new File(outputDirS, "chr" + PStringUtils.getNDigitNumber(3, chrIndex + 1) + "_vmap2.1_AnnoDB_addDAF_addSIFT.txt.gz").getAbsolutePath();

            String header = null;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath()); //是xls文件，故不是压缩的
                String temp = br.readLine(); //read header
                List<SIFT.SIFTRecord> sList = new ArrayList<>(); //建立sift记录的集合；建立sift记录的集合；建立sift记录的集合；建立sift记录的集合；建立sift记录的集合；建立sift记录的集合
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
//CHROM	POS	REF_ALLELE	ALT_ALLELE	TRANSCRIPT_ID	GENE_ID	GENE_NAME	REGION	VARIANT_TYPE	REF_AMINO	ALT_AMINO AMINO_POS	SIFT_SCORE	SIFT_MEDIAN	NUM_SEQS	dbSNP	SIFT_PREDICTION
//6	28419	G	A	TraesCS1D02G374800.1	TraesCS1D02G374800	TraesCS1D02G374800	CDS	NONSYNONYMOUS	M	I
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String alt = l.get(3);
                    String transcript = l.get(4);
                    String type = l.get(8);
                    String value = l.get(12);
                    SIFT.SIFTRecord s = new SIFT.SIFTRecord(pos, alt, transcript, type, value); //将xls文件中的sift结果放入集合中，建立库文件
                    sList.add(s);
                }
                br.close();
                Collections.sort(sList);

                //全部修改成根据文件后缀来识别压缩或者非压缩。
                if (dbFileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(dbFileS);
                } else if (dbFileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(dbFileS);
                }

                header = br.readLine() + "\tVariant_type\tSIFT_score\tTranscript";
                BufferedWriter bw = null;
                if(outfileS.endsWith(".txt")){
                    bw = IOUtils.getTextWriter(outfileS);
                }
                else if(outfileS.endsWith(".txt.gz")){
                    bw = IOUtils.getTextGzipWriter(outfileS);
                }
                bw.write(header);
                bw.newLine();
                int index = -1;
                int startIndex = -1;
                int endIndex = -1;
                StringBuilder sb = null;
                while ((temp = br.readLine()) != null) {
                    sb = new StringBuilder();
                    sb.append(temp).append("\t");
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String alt = l.get(3);
                    SIFT.SIFTRecord query = new SIFT.SIFTRecord(pos, alt, "", "", "");
                    index = Collections.binarySearch(sList, query);
                    if (index < 0) {

                        sb.append("NA").append("\t").append("NA").append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }
                    //意思是在库中，只有有sift类型的pos才会被列出来，但是在我们的数据库中，所有编码区非编码区都会被列出来，所有我们的数据库中pos一般比sift库中小
                    startIndex = index; //不是太懂  startIndex = -1;  startIndex - 1 = -2 < -1
                    endIndex = index; //不是太懂  endIndex = -1;  endIndex + 1 = 0 < sList.size
                    while ((startIndex - 1) > -1 && sList.get(startIndex - 1).isSimilar(pos, alt)) {
                        startIndex--;
                    }
                    while ((endIndex + 1) < sList.size() && sList.get(endIndex + 1).isSimilar(pos, alt)) {
                        endIndex++;
                    }
                    boolean status = false;
                    for (int i = startIndex; i < endIndex + 1; i++) { //对sList库的所有点进行循环，注意结果中可能有一个点很多个转录本结果，我们只选最长的转录本
                        if (Arrays.binarySearch(trans, sList.get(i).transcript) >= 0 && sList.get(i).type.equals("NONSYNONYMOUS") || sList.get(i).type.equals("SYNONYMOUS")
                                || sList.get(i).type.equals("START-LOST") || sList.get(i).type.equals("STOP-GAIN") || sList.get(i).type.equals("STOP-LOSS")
                                || sList.get(i).type.equals("NONCODING") || sList.get(i).type.equals("FRAMESHIFT INSERTION")) {
                            if (Arrays.binarySearch(trans, sList.get(i).transcript) < 0) {
                                sb.append("NA\tNA\tNA");
                            } else {
                                sb.append(sList.get(i).type).append("\t").append(sList.get(i).value).append("\t").append(sList.get(i).transcript);
                            }

                            bw.write(sb.toString());
                            bw.newLine();
                            status = true;
                            break;
                        }
                    }
                    if (status == false) {
                        sb.append("NA\tNA\tNA");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    class SIFTRecord implements Comparable<SIFTRecord> {

        public int pos;
        public String alt;
        public String transcript;
        public String type;
        public String value;

        public SIFTRecord(int pos, String alt, String transcript, String type, String value) {
            this.pos = pos;
            this.alt = alt;
            this.transcript = transcript;
            this.type = type;
            this.value = value;
        }

        public boolean isSimilar(int pos, String alt) {
            if (pos == this.pos && alt.equals(this.alt)) {
                return true;
            }
            return false;
        }

        @Override
        public int compareTo(SIFTRecord o) { //
            if (this.pos < o.pos) { 
                return -1;
            } else if (this.pos == o.pos) {
                return this.alt.compareTo(o.alt);
            } else {
                return 1;
            }

        }
    }

    public boolean ifMatch(List<String> l, int pos, String alt, String[] trans) {
        int posS = Integer.parseInt(l.get(1));
        if (posS != pos) {
            return false;
        }
        if (!alt.equals(l.get(3))) {
            return false;
        }
        String type = l.get(8);
        if (!(type.startsWith("NONSYNONYMOUS") || type.startsWith("SYNONYMOUS"))) {
            return false;
        }
        String query = l.get(4);
        if (Arrays.binarySearch(trans, query) < 0) {
            return false;
        }
        return true;
    }

    /**
     * find -name "*.vcf" | cut -c3- ; 本地获取运行脚本，输出在netbeans的output界面。 -c	To run
     * on command line -i	Path to your input variants file in VCF format -d	Path
     * to SIFT database directory -r Path to your output results folder -t	To
     * extract annotations for multiple transcripts (Optional)
     * 输出的文件产生的日志也保存下来，便于计算时间。 注意输入文件必须是解压后的vcf
     */
    public void annotatorVCF(String infileDirS, String outfileDirS) {
        //String infileDirS = "/data4/home/aoyue/vmap2/analysis/008_sift/dbSNP20190820/ab/";
        //String outfileDirS = "/data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/output/ab/";

        try {
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/002_annotatorVCF/001_fileList/list_chrABDgenomeTest.txt";
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                String chr = temp.substring(3, 6);
                System.out.println("java -jar /data1/home/aoyue/biosoftware/SIFT4g_Annotator/SIFT4G_Annotator_v2.4.jar -c -i "
                        + new File(infileDirS, temp).getAbsolutePath()
                        + " -d /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/siftDatabase -r "
                        + new File(outfileDirS, "output" + chr).getAbsolutePath()
                        + " -t > "
                        + new File(outfileDirS, "log_" + temp + ".annotator.txt").getAbsolutePath()
                        + " &");
                //java -jar /data1/home/aoyue/biosoftware/SIFT4g_Annotator/SIFT4G_Annotator_v2.4.jar -c -i /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input001/dbSNP/chr036.Dgenome.bi.vcf.gz -d /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input001/abd_iwgscV1 -r /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/output/output036 -t &
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mvDatabase() {
        String[] num = new String[42];
        String[] num2 = new String[42];
        for (int i = 0; i < num.length; i++) {
            num[i] = String.valueOf(i + 1);
            num2[i] = String.valueOf(PStringUtils.getNDigitNumber(3, i + 1));
            System.out.println("cp -f /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input" + num2[i] + "/abd_iwgscV1/" + num[i] + ".gz /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/siftDatabase/ ");
            System.out.println("cp -f /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input" + num2[i] + "/abd_iwgscV1/" + num[i] + ".regions /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/siftDatabase/ ");
            System.out.println("cp -f /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input" + num2[i] + "/abd_iwgscV1/" + num[i] + "_SIFTDB_stats.txt /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/siftDatabase/ ");
        }

    }

    public void mkJavaCmdchr1_42_SIFTdb() {
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/script/";
        String shfileS = new File("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/", "sh_mkSIFTdb.sh").getAbsolutePath();
        //nohup java -Xms200g -Xmx500g -jar FastCall.jar parameters_001_FastCall.txt > log_001_fastcall.txt &
        try {

            for (int i = 1; i < 43; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i);
                String outfileS = new File(outfileDirS, "chr" + chr + "_create_wheat_sift4g_db.sh").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("annotationHeader\n"
                        + "echo \"chr" + chr + "    :Create wheat sif4G database begins\"" + "\n"
                        + "perl make-SIFT-db-all.pl -config /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/config/abd_iwgscV1_SIFT4G_Create_Genomic_DB_config" + chr + ".txt > log_chr" + chr + ".wheat_sift4G_pipeline.step.txt\n"
                        + "annotationHeader\n"
                        + "echo \"chr" + chr + "    :Create wheat sif4G database ends\"" + "\n");

                bw.flush();
                bw.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try {
            File[] fs = new File(outfileDirS).listFiles();
            fs = IOUtils.listFilesEndsWith(fs, ".sh");
            Arrays.sort(fs);
            BufferedWriter bw = IOUtils.getTextWriter(shfileS);
            for (int i = 0; i < fs.length; i++) {
                bw.write("sh " + fs[i].getName() + " > alog_" + fs[i].getName() + " &");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkParameterConfigChr1_42() {
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/config/";
        try {
            for (int i = 1; i < 43; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i);
                String outfileS = new File(outfileDirS, "abd_iwgscV1_SIFT4G_Create_Genomic_DB_config" + chr + ".txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                sb.append("\n");
                sb.append("GENETIC_CODE_TABLE=1\n"
                        + "GENETIC_CODE_TABLENAME=Standard\n"
                        + "MITO_GENETIC_CODE_TABLE=11\n"
                        + "MITO_GENETIC_CODE_TABLENAME=Plant Plastids\n").append("\n");
                sb.append("PARENT_DIR=/data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input").append(chr).append("\n");
                sb.append("ORG=bread_wheat\n"
                        + "ORG_VERSION=abd_iwgscV1\n"
                        + "DBSNP_VCF_FILE=abd_iwgscV1.vcf.gz\n").append("\n");
                sb.append("#Running SIFT 4G\n"
                        + "SIFT4G_PATH=/data1/home/aoyue/biosoftware/SIFT4G/sift4g/bin/sift4g\n"
                        + "PROTEIN_DB=/data1/publicData/protein_db/uniprot90_Jul2018/uniref90.fasta\n"
                        + "COMPUTER=GIS-KATNISS\n"
                        + "\n"
                        + "\n"
                        + "# Sub-directories, don't need to change\n"
                        + "GENE_DOWNLOAD_DEST=gene-annotation-src\n"
                        + "CHR_DOWNLOAD_DEST=chr-src\n"
                        + "LOGFILE=Log.txt\n"
                        + "ZLOGFILE=Log2.txt\n"
                        + "FASTA_DIR=fasta\n"
                        + "SUBST_DIR=subst\n"
                        + "ALIGN_DIR=SIFT_alignments\n"
                        + "SIFT_SCORE_DIR=SIFT_predictions\n"
                        + "SINGLE_REC_BY_CHR_DIR=singleRecords\n"
                        + "SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores\n"
                        + "DBSNP_DIR=dbSNP\n"
                        + "\n"
                        + "# Doesn't need to change\n"
                        + "FASTA_LOG=fasta.log\n"
                        + "INVALID_LOG=invalid.log\n"
                        + "PEPTIDE_LOG=peptide.log\n"
                        + "ENS_PATTERN=ENS\n"
                        + "SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord\n"
                        + "");
                bw.write(sb.toString());
                bw.newLine();
                bw.flush();
                bw.close();

            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkInputDirs() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/003_inputFileList/";
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            File f = new File(infileDirS, "input" + chr);
            f.mkdirs();
            File chrsrc = new File(f, "chr-src");
            chrsrc.mkdirs();
            File geneAnnotationSrc = new File(f, "gene-annotation-src");
            geneAnnotationSrc.mkdirs();
            File dbSNP = new File(f, "dbSNP");
            dbSNP.mkdirs();

            //输出命令的脚本
            try {
                String scriptS = new File(infileDirS, chr + ".sh").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
//                bw.write("cp -f /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/002_splitGTF/chr" + chr +
//                        ".wheat_v1.1_Lulab.addBiotype.gtf.gz " 
//                + geneAnnotationSrc.getAbsolutePath());

                bw.write("cp -f /data1/publicData/wheat/reference/v1.0/byChr/chr" + chr
                        + ".fa.gz /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_d/input" + chr + "/chr-src/ &");
                bw.newLine();
                bw.flush();
                bw.close();

                StringBuilder sb = new StringBuilder("sh ");
                sb.append(scriptS);
                Runtime run = Runtime.getRuntime();
                Process p = run.exec(sb.toString());
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                StringBuilder ssb = new StringBuilder(); //有时候需要将屏幕信息输出到一个文件中，此时建立输出文本路径。
                String line = null;
                while ((line = br.readLine()) != null) {
                    System.out.println(line);
                    ssb.append(line + "\n");
                }
                String result = ssb.toString(); //有必要时输出出来
                p.waitFor();
                System.out.println(sb.toString());
                System.out.println(result);

                //new File(scriptS).delete();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

            //cat /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_chrList/chrABDgenomeList.txt | awk '{print $1}' | while read a;do echo "cp -f /data1/publicData/wheat/reference/v1.0/byChr/chr${a}.fa.gz /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input${a}/chr-src/ &";done
        }

    }

    /**
     * 将gtf文件拆分成单条染色体
     */
    public void parseGff3() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/001_gtf/wheat_v1.1_Lulab.addBiotype.gtf.gz";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/002_splitGTF/";
        String[] outfileS = new String[43];
        //String OutmtFileS = "";
        //"/Users/Aoyue/Documents/Zea_mays.AGPv4.38.gtf";
        for (int i = 0; i < outfileS.length; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            outfileS[i] = new File(outfileDirS, "chr" + chr + ".wheat_v1.1_Lulab.addBiotype.gtf.gz").getAbsolutePath();
            System.out.println(outfileS[i] + "\n");
        }
        try {
            for (int i = 0; i < outfileS.length; i++) {
                String chr = Integer.toString(i);
                BufferedReader br = IOUtils.getTextGzipReader(infileS);
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS[i]);
                String tem = null;
                List<String> l = null;
                while ((tem = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(tem);
                    if (l.get(0).equals(chr)) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(tem);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * modify gtf file, add gene_biotype to the 9th colum
     */
    public void addGeneBiotype() {
        String InputFileS = "/Users/Aoyue/Documents/wheat_v1.1_Lulab.gtf";
        String OutputFileS = "/Users/Aoyue/Documents/wheat_v1.1_Lulab.addBiotype.gtf";
        try {
            BufferedReader br = IOUtils.getTextReader(InputFileS);
            BufferedWriter bw = IOUtils.getTextWriter(OutputFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                //String plate = null;
                //plate = PStringUtils.fastSplit(temp).get(2);
                //if(plate.equalsIgnoreCase("exon") | plate.equalsIgnoreCase("CDS") | plate.equalsIgnoreCase("stop_codon") | plate.equalsIgnoreCase("start_codon")){
                List<String> l = null;
                l = PStringUtils.fastSplit(temp);
                StringBuilder sb = new StringBuilder();
                sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(l.get(3)).append("\t").
                        append(l.get(4)).append("\t").append(l.get(5)).append("\t").
                        append(l.get(6)).append("\t").append(l.get(7)).append("\t").append(l.get(8)).append(" ")
                        .append("gene_biotype").append(" ").append("\"protein_coding\";");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
