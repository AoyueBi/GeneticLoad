package PopulationAnalysis;

import AoUtils.AoFile;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import pgl.infra.window.SimpleWindow;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 * @author AoyueBi
 * @data 2020-10-22 08:39
 */
public class DeleteriousXPCLR2 {

    public DeleteriousXPCLR2(){
        this.pipeDeleteriousXPCLR();
//        this.pipeFilterTaxa();
//        this.getCDSlength();

    }

    public void pipeDeleteriousXPCLR(){
//        this.step0(); //第一次计算
//        this.step1(); //第二次计算
        this.step3(); //第三次计算:ref bias evaluation 后，load 变化



    }

    public void step3(){
        //************************************** new  2020-12-26 **************************************
//**************************** 不必修改 **************************//
        String variantType = null;
        String variantType1 = "001_synonymous";
        String variantType2 = "002_nonsynonymous";
        String variantType3 = "003_nonsynGERPandDerivedSIFT";
        String variantType4 = "004_nonsynDerivedSIFT";
        String variantType5 = "005_GERP";
        String variantType6 = "006_nonsynGERPandDerivedSIFT_correction";
        String variantType7 = "007_nonsynDerivedSIFT_correction";
        String variantType8 = "008_GERP_correction";

        String ratioType = null;
        String ratioType1 = "bySub";
        String ratioType2 = "bysub_mergeByTaxa";
        String group = null;

        String[] choice1 = {variantType1,variantType2,variantType3,variantType4,variantType5,variantType6,variantType7,variantType8};
        String[] choice3 = {ratioType1,ratioType2};

        //**************************** up 不必修改 **************************//
        String parentDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/001_delCount";
        new File(parentDirS).mkdirs();
        for (int k = 0; k < choice1.length; k++) {
            variantType = choice1[k];
            String addCountFileS = new File(parentDirS,variantType + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr.txt").getAbsolutePath();
            this.countDeleteriousVMapII_byChr_refBiasEvaluation(variantType,addCountFileS);
        }

//        int cnt = 0;
//        for (int k = 1; k < choice1.length; k++) {
//            variantType = choice1[k];
//            for (int l = 0; l < choice3.length; l++) {
//                ratioType = choice3[l];
//                String infileS1 = new File(parentDirS,variantType + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
//                String infileS2 = new File(parentDirS,"001_synonymous" + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType +".txt").getAbsolutePath();
//                this.DeltoSynonymousRatio(infileS1,ratioType,infileS2);
//                cnt++;
//                System.out.println("********" + cnt + " " + variantType + " " + ratioType);
//            }
//        }

//        this.filterLandrace();
    }


    /**
     * 将按照sub和taxa产生的有害突变的结果，进行过滤，使六倍体栽培品种都为欧洲的栽培品种，六倍体的农家种都为欧洲的农家种
     * 然后进行个体Load的判断
     */
    public void filterLandrace(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/001_delCount";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/047_referenceEvaluation/001_delCount_filterLR_CL";
        new File(outfileDirS).mkdirs();

        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        AoFile.readheader(taxaSummaryFileS);
        HashMap<String, String> taxaGenomeTypeMap = new HashMap();
        HashMap<String, String> taxaContinentMap = new HashMap();
        HashMap<String, String> taxaSubspeciesMap = new HashMap();
        RowTable<String> t = new RowTable (taxaSummaryFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGenomeTypeMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 3));
            taxaContinentMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 9));
            taxaSubspeciesMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 15));
        }

        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        for (int i = 0; i < fs.length; i++) {
            File f = fs[i];
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS, f.getName()).getAbsolutePath();
            try {
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String header = br.readLine();
                bw.write(header);bw.newLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    String taxa = l.get(0);
                    //过滤条件，如果是六倍体，如果是landrace，如果是欧洲的，就写出来；如果是cultivar,如果是欧洲的，也写出来； 如果不是六倍体，全都写出来。
                    //故，我需要知道倍性，大洲，树上的亚群信息
                    String genomeType = taxaGenomeTypeMap.get(taxa);
                    String continent = taxaContinentMap.get(taxa);
                    String subspecies = taxaSubspeciesMap.get(taxa);
                    if (genomeType.equals("AABBDD")){
                        if (subspecies.equals("Landrace")){
                            if (continent.equals("Europe")){
                                bw.write(temp);
                                bw.newLine();
                            }
                        }
                        if (subspecies.equals("Cultivar")){
                            bw.write(temp);
                            bw.newLine();
                        }
                    }
                    else {
                        bw.write(temp);
                        bw.newLine();
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }



    /**
     * 获取不同ref-obj对测试下的，100kb 有效cds长度，和非选择区域的cdslength
     */
    public void getCDSlength(){
        String outfileS = "/Users/Aoyue/Documents/cds.txt"; //写出所有CDS list 的 bed 格式，为下文计算 overlap 打基础
        String outfileS2 = "/Users/Aoyue/Documents/gene_cdsLength.txt"; //写出每个基因的cds 长度
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        String hcGeneFileS = "/Users/Aoyue/Documents/Data/wheat/gene/001_geneHC/geneHC.txt";
        RowTable<String> gt = new RowTable<>(hcGeneFileS);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName(); //通过名字排序

        try{
            BufferedWriter bw = AoFile.writeFile(outfileS);
            BufferedWriter bw2 = AoFile.writeFile(outfileS2);

            int geneIndex = -1;
            int tranIndex = -1;
            int cdsStart = -1;
            int cdsEnd = -1;
            for (int j = 0; j < gt.getRowNumber(); j++) { //gt已通过名字进行了排序
                geneIndex = gf.getGeneIndex(gt.getCell(j, 0)); //根据基因名字（不是转录本的名字，没有.后缀12）获取索引
                int chr = gf.getGeneChromosome(geneIndex);
                int geneCDSlength = 0;
                for (int k = 0; k < gf.getTranscriptNumber(geneIndex); k++) { //获取最长转录本的index
                    if (!gt.getCell(j,1).equals(gf.getTranscriptName(geneIndex, k)))continue; //如果 gt.getCell(j,1) 是该基因的最长转录本，已经提前总结出
                    tranIndex = k; //该循环的意思是：如果最长转录本不等于k，那么就跳出循环，一直到等于k为止，最后终止循环。
                    break;
                }
                List<Range> cdsList = gf.getCDSList(geneIndex, tranIndex); //根据 gene index 和 tranIndex 获取该转录本的 cdsList
                for (int k = 0; k < cdsList.size(); k++) {
                    cdsStart = cdsList.get(k).start;
                    cdsEnd = cdsList.get(k).end;
                    geneCDSlength = geneCDSlength + cdsEnd-cdsStart;
                    if (0<chr && chr <43){
                        bw.write(String.valueOf(chr) + "\t" + String.valueOf(cdsStart -1) + "\t" + String.valueOf(cdsEnd -1) + "\t" + String.valueOf( cdsEnd-cdsStart));
                        bw.newLine();
                    }
                }
                bw2.write(chr  + "\t" + gf.getGeneName(geneIndex)+ "\t" + geneCDSlength);
                bw2.newLine();

            }
            bw.flush();bw.close();
            bw2.flush();bw2.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void step0(){
        //**************************** 不必修改 **************************//
        String variantType = null;
        String variantType1 = "001_synonymous";
        String variantType2 = "002_nonsynonymous";
        String variantType3 = "003_nonsynGERPandDerivedSIFT";
        String ifselected = null;
        String ifselected1 = "1";
        String ifselected2 = "0";
        String ratioType = null;
        String ratioType1 = "bySub";
        String ratioType2 = "bysub_mergeByTaxa";
        String group = null;
        String group1 = "wede";
        String group2 = "dedurum";
        String group3 = "lrcul";

        String[] choice1 = {variantType1,variantType2,variantType3};
        String[] choice2 = {ifselected1,ifselected2};
        String[] choice3 = {ratioType1,ratioType2};
        String[] choice4 = {group1,group2,group3};

        //**************************** up 不必修改 **************************//


        //六倍体 LR-CL ，合计运行6次

//        for (int i = 0; i < choice2.length; i++) { //合计2*3=6次循环
//            ifselected = choice2[i];
//            for (int j = 0; j < choice1.length; j++) {
//                String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/008_topK/top0.01_ChrPosTrans.txt.gz"; //受选择的变异位点
//                String parentFileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid";
//                String addCountFileAddGroupS2 = new File(parentFileS1,"009_deleteriousXPCLR").getAbsolutePath(); new File(addCountFileAddGroupS2).mkdirs();
//                variantType = choice1[j];
//                addCountFileAddGroupS2 = new File(addCountFileAddGroupS2,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr.txt").getAbsolutePath();
//                this.countDeleteriousVMapII_byChr(infileS,ifselected,variantType,addCountFileAddGroupS2);
//            }
//        }

        //求六倍体ratio non/syn  del/syn
//        int cnt=0;
//        for (int i = 0; i < choice2.length; i++) { //合计2*2*2=8次循环
//            ifselected = choice2[i];
//            for (int j = 1; j < choice1.length; j++) {
//                variantType = choice1[j];
//                for (int k = 0; k < choice3.length; k++) {
//                    ratioType=choice3[k];
//                    String parentFileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/009_deleteriousXPCLR";
//                    String infileS1 = new File(parentFileS1,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
//                    String infileS2 = new File(parentFileS1, "001_synonymous" + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
//                    this.DeltoSynonymousRatio(infileS1,ratioType,infileS2);
//                    cnt++;
//                    System.out.println("********" + cnt + " " + "ifselected" + ifselected + " " + variantType + " " + ratioType);
//                }
//            }
//        }


        //四倍体 WE-DE
        ///////先单个运行测试，后写成循环进行一次性跑完
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid/008_topK/top0.01_ChrPosTrans.txt.gz"; //受选择的变异位点
//        String parentFileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid";
//        String addCountFileAddGroupS2 = new File(parentFileS1,"009_deleteriousXPCLR").getAbsolutePath(); new File(addCountFileAddGroupS2).mkdirs();
//        variantType = variantType3;
//        String ifselected = "0";
//        addCountFileAddGroupS2 = new File(addCountFileAddGroupS2,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr.txt").getAbsolutePath();
//        this.countDeleteriousVMapII_byChr(infileS,ifselected,variantType,addCountFileAddGroupS2);

        //        for (int i = 0; i < choice2.length; i++) { //合计2*3=6次循环
//            ifselected = choice2[i];
//            for (int j = 0; j < choice1.length; j++) {
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid/008_topK/top0.01_ChrPosTrans.txt.gz"; //受选择的变异位点
//        String parentFileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid";
//                String addCountFileAddGroupS2 = new File(parentFileS1,"009_deleteriousXPCLR").getAbsolutePath(); new File(addCountFileAddGroupS2).mkdirs();
//                variantType = choice1[j];
//                addCountFileAddGroupS2 = new File(addCountFileAddGroupS2,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr.txt").getAbsolutePath();
//                this.countDeleteriousVMapII_byChr(infileS,ifselected,variantType,addCountFileAddGroupS2);
//            }
//        }

        //求四倍体ratio non/syn  del/syn

//        int cnt=0;
//        for (int i = 0; i < choice2.length; i++) { //合计2*2*2=8次循环
//            ifselected = choice2[i];
//            for (int j = 1; j < choice1.length; j++) {
//                variantType = choice1[j];
//                for (int k = 0; k < choice3.length; k++) {
//                    ratioType=choice3[k];
//                    String parentFileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid/009_deleteriousXPCLR";
//                    String infileS1 = new File(parentFileS1,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
//                    String infileS2 = new File(parentFileS1, "001_synonymous" + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
//                    this.DeltoSynonymousRatio(infileS1,ratioType,infileS2);
//                    cnt++;
//                    System.out.println("********" + cnt + " " + "ifselected" + ifselected + " " + variantType + " " + ratioType);
//                }
//            }
//        }


        //四倍体 DE-FTT

//        for (int i = 0; i < choice2.length; i++) { //合计2*3=6次循环
//            ifselected = choice2[i];
//            for (int j = 0; j < choice1.length; j++) {
//                String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/006_tetraploid_FTT_DE/008_topK/top0.01_ChrPosTrans.txt.gz"; //受选择的变异位点
//                String parentFileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/006_tetraploid_FTT_DE";
//                String addCountFileAddGroupS2 = new File(parentFileS1,"009_deleteriousXPCLR").getAbsolutePath(); new File(addCountFileAddGroupS2).mkdirs();
//                variantType = choice1[j];
//                addCountFileAddGroupS2 = new File(addCountFileAddGroupS2,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr.txt").getAbsolutePath();
//                this.countDeleteriousVMapII_byChr(infileS,ifselected,variantType,addCountFileAddGroupS2);
//            }
//        }

        //求四倍体 DE-FTT ratio non/syn  del/syn

//        int cnt=0;
//        for (int i = 0; i < choice2.length; i++) { //合计2*2*2=8次循环
//            ifselected = choice2[i];
//            for (int j = 1; j < choice1.length; j++) {
//                variantType = choice1[j];
//                for (int k = 0; k < choice3.length; k++) {
//                    ratioType=choice3[k];
//                    String parentFileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/006_tetraploid_FTT_DE/009_deleteriousXPCLR";
//                    String infileS1 = new File(parentFileS1,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
//                    String infileS2 = new File(parentFileS1, "001_synonymous" + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
//                    this.DeltoSynonymousRatio(infileS1,ratioType,infileS2);
//                    cnt++;
//                    System.out.println("********" + cnt + " " + "ifselected" + ifselected + " " + variantType + " " + ratioType);
//                }
//            }
//        }




        //************************************** new result after asking the help of HuaChen  **************************************
        //************************************** new result after asking the help of HuaChen  **************************************
        //************************************** new result after asking the help of HuaChen  **************************************


        //四倍体 WE-DE
        ///////先单个运行测试，后写成循环进行一次性跑完
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid/008_topK/top0.01_ChrPosTrans.txt.gz"; //受选择的变异位点
//        String parentFileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid";
//        String addCountFileAddGroupS2 = new File(parentFileS1,"009_deleteriousXPCLR").getAbsolutePath(); new File(addCountFileAddGroupS2).mkdirs();
//        variantType = variantType3;
//        String ifselected = "0";
//        addCountFileAddGroupS2 = new File(addCountFileAddGroupS2,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr.txt").getAbsolutePath();
//        this.countDeleteriousVMapII_byChr(infileS,ifselected,variantType,addCountFileAddGroupS2);

        for (int i = 0; i < choice2.length; i++) { //合计2*3=6次循环  //String[] choice2 = {ifselected1,ifselected2};
            ifselected = choice2[i];
            for (int j = 0; j < choice1.length; j++) { //// String[] choice1 = {variantType1,variantType2,variantType3}; that is 001_synon 002_nonsyn 003_deleterios
                String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid/008_topK/top0.01_ChrPosTrans.txt.gz"; //受选择的变异位点
                String parentFileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid";
                String addCountFileAddGroupS2 = new File(parentFileS1,"009_deleteriousXPCLR").getAbsolutePath(); new File(addCountFileAddGroupS2).mkdirs();
                variantType = choice1[j];
                addCountFileAddGroupS2 = new File(addCountFileAddGroupS2,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr.txt").getAbsolutePath();
                this.countDeleteriousVMapII_byChr(infileS,ifselected,variantType,addCountFileAddGroupS2);
            }
        }

        //求四倍体ratio non/syn  del/syn

        int cnt=0;
        for (int i = 0; i < choice2.length; i++) { //合计2*2*2=8次循环
            ifselected = choice2[i];
            for (int j = 1; j < choice1.length; j++) {
                variantType = choice1[j];
                for (int k = 0; k < choice3.length; k++) {
                    ratioType=choice3[k];
                    String parentFileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid/009_deleteriousXPCLR";
                    String infileS1 = new File(parentFileS1,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
                    String infileS2 = new File(parentFileS1, "001_synonymous" + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
                    this.DeltoSynonymousRatio(infileS1,ratioType,infileS2);
                    cnt++;
                    System.out.println("********" + cnt + " " + "ifselected" + ifselected + " " + variantType + " " + ratioType);
                }
            }
        }
    }


    public void step1(){
        //************************************** new  2020-12-26 **************************************
//**************************** 不必修改 **************************//
        String variantType = null;
        String variantType1 = "001_synonymous";
        String variantType2 = "002_nonsynonymous";
        String variantType3 = "003_nonsynGERPandDerivedSIFT";
        String ifselected = null;
        String ifselected1 = "1";
        String ifselected2 = "0";
        String ratioType = null;
        String ratioType1 = "bySub";
        String ratioType2 = "bysub_mergeByTaxa";
        String group = null;
        String group1 = "wede";
        String group2 = "dedurum";
        String group3 = "lrcul";

        String[] choice1 = {variantType1,variantType2,variantType3};
        String[] choice2 = {ifselected1,ifselected2};
        String[] choice3 = {ratioType1,ratioType2};
        String[] choice4 = {group1,group2,group3};

        //**************************** up 不必修改 **************************//
        String parentDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/014_summary_XPCLR/003_deleteriousXPCLR";
//        new File(parentDirS).mkdirs();
//        for (int i = 0; i < choice4.length; i++) {
//            group = choice4[i];
//            String infileS = new File("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/014_summary_XPCLR/002_topK","top0.05_" + group + "_ChrPos_fromExonAnnotation.txt.gz").getAbsolutePath();
//            for (int j = 0; j < choice2.length; j++) {
//                ifselected = choice2[j];
//                for (int k = 0; k < choice1.length; k++) {
//                    variantType = choice1[k];
//                    String addCountFileS = new File(parentDirS,variantType + "_ifselected" + ifselected + "_" + group + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr.txt").getAbsolutePath();
//                    this.countDeleteriousVMapII_byChr(infileS,ifselected,variantType,addCountFileS);
//                }
//            }
//        }

        int cnt = 0;
        for (int i = 0; i < choice4.length; i++) {
            group = choice4[i];
            for (int j = 0; j < choice2.length; j++) {
                ifselected = choice2[j];
                for (int k = 1; k < choice1.length; k++) {
                    variantType = choice1[k];
                    for (int l = 0; l < choice3.length; l++) {
                        ratioType = choice3[l];
                        String infileS1 = new File(parentDirS,variantType + "_ifselected" + ifselected + "_" + group + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
                        String infileS2 = new File(parentDirS,"001_synonymous" + "_ifselected" + ifselected + "_" + group + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType +".txt").getAbsolutePath();
                        this.DeltoSynonymousRatio(infileS1,ratioType,infileS2);
                        cnt++;
                        System.out.println("********" + cnt + " " + "ifselected" + ifselected + " " + variantType + " " + ratioType + " " + group);

//addCountFileAddGroupS2 = new File(addCountFileAddGroupS2,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr.txt").getAbsolutePath();
//String infileS1 = new File(parentFileS1,variantType + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();
//String infileS2 = new File(parentFileS1, "001_synonymous" + "_ifselected" + ifselected + "_additiveDeleterious_ANCbarleyVSsecalePasimony_vmap2_bychr_" + ratioType + ".txt").getAbsolutePath();

                    }
                }
            }
        }
    }


    /**
     * @deprecated
     * this step is done by R language
     */
    public void pipeFilterTaxa(){

        String[] choice1 = {"WE_DE","DE_FTT","LR_CL"};
        String[] choice2 = {"005_tetraploid","006_tetraploid_FTT_DE","004_hexaploid"};
        String parentFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR";
        for (int i = 0; i < choice1.length; i++) {
            String infileDirS1 = new File(parentFileS,choice2[i]).getAbsolutePath();
            String infileDirS = new File(infileDirS1,"009_deleteriousXPCLR").getAbsolutePath();
            String outfileDirS = new File(infileDirS1,"010_deleteriousXPCLR_keep_" + choice1[i]).getAbsolutePath(); new File(outfileDirS).mkdirs();
            File[] fs = AoFile.getFileArrayInDir(infileDirS);
            for (int j = 0; j < fs.length; j++) {
                String infileS = fs[j].getAbsolutePath();
                String outfileS = new File(outfileDirS,fs[j].getName()).getAbsolutePath();
                this.filterTaxa(infileS,choice1[i],outfileS);
            }
        }
    }

    /**
     * 将按照sub和taxa产生的有害突变的结果，进行过滤，使六倍体栽培品种都为欧洲的栽培品种，六倍体的农家种都为欧洲的农家种
     * 然后进行个体Load的判断
     */
    public void filterTaxa(String infileS, String type,String outfileS){
        List<String> l1wede = new ArrayList<>();
        List<String> l2deftt = new ArrayList<>();
        List<String> l3lrcl = new ArrayList<>();


        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        AoFile.readheader(taxaSummaryFileS);


        RowTable<String> t = new RowTable (taxaSummaryFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String taxa = t.getCell(i,0);
            String subspecies = t.getCell(i,15);
            String query1 = t.getCell(i,20);
            if (type.equals("WE_DE")){
                if (subspecies.equals("Wild_emmer") || subspecies.equals("Domesticated_emmer")){
                    l1wede.add(taxa);
                    AoFile.filterTxtLines(infileS,0, l1wede,outfileS);
                }
            }
            if (type.equals("DE_FTT")){
                if (subspecies.equals("Domesticated_emmer") || subspecies.equals("Free_threshing_tetraploid")){
                    l2deftt.add(taxa);
                    AoFile.filterTxtLines(infileS,0, l2deftt,outfileS);
                }
            }
            if (type.equals("LR_CL")){
                if (query1.equals("LR_EU") || query1.equals("Cultivar")){
                    l3lrcl.add(taxa);
                    AoFile.filterTxtLines(infileS,0, l3lrcl,outfileS);
                }
            }
        }
    }



    public void countDeleteriousVMapII_byChr_refBiasEvaluation(String type, String addCountFileAddGroupS) {

        //model

        //######## 需要修改 ########//
//        String exonVCFDirS = ""; //外显子变异数据
//        String SNPAnnoFileS = ""; //注释信息库合并后的总文件
//        AoFile.readheader(SNPAnnoFileS);

        int cntNONSY = 0; //非同义突变的个数

        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################w
        String exonVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF"; //外显子变异数据
        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //注释信息库合并后的总文件
        AoFile.readheader(SNPAnnoFileS);
        //************* 无需修改的路径 ****************** //
        String addCountFileS = new File(addCountFileAddGroupS).getAbsolutePath().replaceFirst(".txt",".temp.txt"); //有害变异加性模型输出文件

        /**
         *  ################################### step1: 初始化染色体集合
         */
        int minDepth = 2;//inclusive
        int chrNum = 42;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrList.add(i + 1);
        }
        System.out.println("Finished step1: completing the initialization of chromosome.");

        /**
         *  ################################### step2: posList  charList 补充完整
         */

        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];

        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }

        String derivedAllele = null;

        try {
            BufferedReader br = AoFile.readFile(SNPAnnoFileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int index = Integer.parseInt(l.get(1)) - 1; //染色体号的索引 ################ 需要修改 需要修改 需要修改 ################
                int pos = Integer.parseInt(l.get(2)); //################ 需要修改 需要修改 需要修改 ################
                String variantType = l.get(12); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(16); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(20); //################ 需要修改 需要修改 需要修改 ################

                //********************* 过滤没有 ancestral allele 信息的位点
                //################### 需要修改 //###################//###################//###################//###################
                ////不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
                String ancestralAllele = l.get(15);
                //################### 需要修改 //###################//###################//###################//###################
                String ref = l.get(3);
                String alt = l.get(4);
                String majorAllele = l.get(5);
                String minorAllele = l.get(6);
                if(!(ancestralAllele.equals(ref)||ancestralAllele.equals(alt)))continue;


                /**
                 ******** 定义有害突变，不是有害突变，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */

                if (type.equals("001_synonymous")){
                    if (!variantType.equals("SYNONYMOUS"))continue;
                }
                if (type.equals("002_nonsynonymous")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue;
                }
                if (type.equals("003_nonsynGERPandDerivedSIFT")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                    if (sift.startsWith("N"))continue; //说明必须满足有sift值
                    double gerpd = Double.parseDouble(gerp);
                    double siftd = Double.parseDouble(sift);
                    if (gerpd < 1) continue; //说明必须满足gerp大于1
                    if (siftd > 0.05) continue; //说明必须满足sift小于等于0.05
                }

                if (type.equals("004_nonsynDerivedSIFT")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if (sift.startsWith("N"))continue; //说明必须满足有sift值
                    double siftd = Double.parseDouble(sift);
                    if (siftd > 0.04) continue; //说明必须满足sift小于等于0.05
                }

                if (type.equals("005_GERP")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                    double gerpd = Double.parseDouble(gerp);
                    if (gerpd < 1) continue; //说明必须满足gerp大于1
                }

                if (type.equals("006_nonsynGERPandDerivedSIFT_correction")){
                    if (ancestralAllele.equals(ref)){ //
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                        if (sift.startsWith("N"))continue; //说明必须满足有sift值
                        double gerpd = Double.parseDouble(gerp);
                        double siftd = Double.parseDouble(sift);
                        if (gerpd < 1.61) continue; //说明必须满足gerp大于1
                        if (siftd > 0.04) continue; //说明必须满足sift小于等于0.05
                    }
                    if (ancestralAllele.equals(alt)){
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                        if (sift.startsWith("N"))continue; //说明必须满足有sift值
                        double gerpd = Double.parseDouble(gerp);
                        double siftd = Double.parseDouble(sift);
                        if (gerpd < 0.0139) continue; //说明必须满足gerp大于1
                        if (siftd > 0.3) continue; //说明必须满足sift小于等于0.05
                    }
                }

                if (type.equals("007_nonsynDerivedSIFT_correction")){
                    if (ancestralAllele.equals(ref)){ //
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if (sift.startsWith("N"))continue; //说明必须满足有sift值
                        double siftd = Double.parseDouble(sift);
                        if (siftd > 0.04) continue; //说明必须满足sift小于等于0.05
                    }
                    if (ancestralAllele.equals(alt)){
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if (sift.startsWith("N"))continue; //说明必须满足有sift值
                        double siftd = Double.parseDouble(sift);
                        if (siftd > 0.3) continue; //说明必须满足sift小于等于0.05
                    }
                }

                if (type.equals("008_GERP_correction")){
                    if (ancestralAllele.equals(ref)){ //
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                        double gerpd = Double.parseDouble(gerp);
                        if (gerpd < 1.61) continue; //说明必须满足gerp大于1
                    }
                    if (ancestralAllele.equals(alt)){
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                        double gerpd = Double.parseDouble(gerp);
                        if (gerpd < 0.0139) continue; //说明必须满足gerp大于1
                    }
                }


                if (ancestralAllele.equals(majorAllele)) {
                    derivedAllele = minorAllele;
                    posList[index].add(pos); //将包含有derived allele的位点添加到Poslist
                    charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
                }
                if (ancestralAllele.equals(minorAllele)) {
                    derivedAllele = majorAllele;
                    posList[index].add(pos);
                    charList[index].add(derivedAllele.charAt(0));
                }
                else if (!(ancestralAllele.equals(majorAllele) || ancestralAllele.equals(minorAllele))){
                }
                cnt++;
                if (cnt%100000==0) System.out.println("cnt is " + cnt +" going on at step2");

            }
            br.close();
            System.out.println(cntNONSY + " nonsynonymous SNP num");

            for (int i = 0; i < chrNum; i++) { //将每一个list转化为数组
                delePos[i] = posList[i].toArray();
                deleChar[i] = charList[i].toArray();
                Arrays.sort(delePos[i]);
            }
            System.out.println("Finished step2: completing the posList  charList.");

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        /**
         *  ################################### step3: taxa 集合 642个taxa
         */
        String vmap2TaxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";

        String[] taxa = AoFile.getStringArraybyList(vmap2TaxaList,0);

        double[][] addCount = new double[chrNum][taxa.length];
        int[][] recCount = new int[chrNum][taxa.length];
        int[][] siteWithMinDepthCount = new int[chrNum][taxa.length]; //每个taxa在每条染色体中的有害突变位点

        chrList.parallelStream().forEach(chr -> {
            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_exon_vmap2.1.vcf.gz";
            //开始读写VCF文件
            delVmapFileS = new File(exonVCFDirS, delVmapFileS).getAbsolutePath();
            BufferedReader br = AoFile.readFile(delVmapFileS);
            int chrIndex = chr - 1;
            try {
                String temp = null;
                int cnt = 0;
                HashMap<String, Integer> hmtaxainVCFindex = new HashMap<>();
                HashMap<String, Integer> hmtaxainTaxaindex = new HashMap<>();
                //这里涉及2个数组概念， taxa and taxainVCFfile,即总的taxa和在VCF文件中的taxa
                String[] taxainVCFfile = null;
                List<String> taxainVCFlist = new ArrayList();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##")) {
                        continue;
                    }
                    if (temp.startsWith("#CHROM")) {//说明进入taxa信息列
                        List<String> l = PStringUtils.fastSplit(temp, "\t");
                        taxainVCFfile = new String[l.size() - 9];
                        for (int i = 9; i < l.size(); i++) {
                            taxainVCFfile[i - 9] = l.get(i);
                            hmtaxainVCFindex.put(taxainVCFfile[i - 9], i); //第0个taxa的genotype在第9行，第一个taxa的genotype在第10行，依次类推；
                        }
                        // ************** 在总的taxa中，搜索VCF文件中的taxa的index
                        for (int i = 0; i < taxainVCFfile.length; i++) { //找到taxa
                            int index = Arrays.binarySearch(taxa, taxainVCFfile[i]);
                            if (index > -1) {
                                hmtaxainTaxaindex.put(taxainVCFfile[i], index);
                            }
                        }
                        taxainVCFlist = Arrays.asList(taxainVCFfile);
                    }

                    if (!temp.startsWith("#")) {
                        cnt++;
                        if (cnt % 10000 == 0) {
                            System.out.println(String.valueOf(cnt/1000) + " kb lines on chr " + String.valueOf(chr));
                        }
                        List<String> l = PStringUtils.fastSplit(temp.substring(0, 100), "\t");
                        int pos = Integer.valueOf(l.get(1)); //根据pos找到index
                        int index = Arrays.binarySearch(delePos[chrIndex], pos);
                        if (index < 0) {
                            continue; //前面建立的 delePos[][] 和deleChar[][] 都是为现在在vcf文件中找位置贡献的，不是有害突变的位点，都过滤。*************************************************************
                        } //重要：假如该位点有基因型，也是非同义突变位点，但是没有祖先位点状态，于是就过滤。
                        l = PStringUtils.fastSplit(temp, "\t");
                        int[] idx = new int[2];
                        if (l.get(3).charAt(0) == deleChar[chrIndex][index]) { //如果ref allele = deleChar allele
                            idx[0] = 0; //等于0则代表是 derived allele
                            idx[1] = 1;
                        } else {
                            idx[0] = 1;
                            idx[1] = 0;
                        }


                        //合计642个taxa，在A Bgenome中只有606（419+187）个，在Dgenome中只有455（419+36）个，我们要找到每个VCF文件中的genotype所对应的taxa的index
                        //code:本段代码是每行SNP位点，每个taxa的有害位点的统计
                        for (int i = 0; i < taxainVCFlist.size(); i++) {
                            int genotypeIndex = hmtaxainVCFindex.get(taxainVCFlist.get(i)); //获取该taxa的基因型所在的列的索引
                            int taxaIndex = hmtaxainTaxaindex.get(taxainVCFlist.get(i)); //获取该taxa所在总的642个数组中的索引，为后续写文件进行统计
                            String genoS = l.get(genotypeIndex); //指的是 GT:AD:GL 信息
                            if (genoS.startsWith(".")) {
                                continue; //如果以.开头，说明没有基因型信息，此位点没有测到。
                            }
                            List<String> ll = PStringUtils.fastSplit(genoS, ":"); //分开为GT   AD   GL三类

                            List<String> lll = PStringUtils.fastSplit(ll.get(1), ","); //lll指将AD提取出来，并以"，"号分割。如 0/0:1,0:0,3,28中 ，1，0分别代表ref和alt的测序深度
                            int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1)); //总得测序深度等于 ref + alt
                            if (depth < minDepth) {
                                continue; //最小测序深度是2，如果小于2，则弃用
                            }
                            lll = PStringUtils.fastSplit(ll.get(0), "/"); //这里lll指的是基因型GT，lll被重新赋值，之前代表的是AD
                            int v1 = Integer.valueOf(lll.get(0)); //v1等于 ref
                            int v2 = Integer.valueOf(lll.get(1)); // v2 等于 alt
                            int sum = 0;
                            //如果ref是derived allele，那么idx[0]=0;idx[1]=1. 当是0/0时，sum=2; 0/1时，sum=1; 1/1时， sum=0.
                            //如果alt是derived allele，那么idx[0]=1;idx[1]=0. 当是0/0时，sum=0; 0/1时，sum=1; 1/1时， sum=2.
                            //
                            //
                            if (v1 == idx[0]) { //
                                sum++;
                            }
                            if (v2 == idx[0]) {
                                sum++;
                            }
                            if (sum == 0) {
                            } else if (sum == 1) {
                                addCount[chrIndex][taxaIndex] += 0.5;
                            } else {
                                addCount[chrIndex][taxaIndex] += 1;
                                recCount[chrIndex][taxaIndex] += 1;
                            }
                            siteWithMinDepthCount[chrIndex][taxaIndex]++; //taxa有多少个有害突变位点
                        }
                    }
                }
                br.close();
                System.out.println("Finished step3: calculating each taxa for each chromosome.");

            } catch (Exception e) {
                e.printStackTrace();
            }

            /**
             * 进行每个样品每条染色体的计算
             */
        });

        /**
         *  ################################### step4: 输出每个样品每条染色体的结果
         */

        try {
            BufferedWriter bw = AoFile.writeFile(addCountFileS);
            bw.write("Taxa\tChr\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) { //第一层是染色体号
                String chr = String.valueOf(i + 1);
                for (int j = 0; j < addCount[0].length; j++) { //第二层是该号染色体的有害变异个数
                    bw.write(taxa[j] + "\t" + chr + "\t" + String.valueOf(addCount[i][j]) + "\t" + String.valueOf(siteWithMinDepthCount[i][j]));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();

            System.out.println("Finished step4: writing taxa table.");

        } catch (Exception e) {
            e.printStackTrace();
        }


        /**
         *  ################################### step5: 将输出的文件添加分组信息：每个taxa的倍性，亚群分布，mutation burden Index
         */
        try{
            String taxaSummaryFileS = vmap2TaxaList;
//            String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
            AoFile.readheader(taxaSummaryFileS);
            double depthCut = 1; //只保留深度大于3的taxa样本
            HashMap<String, String> taxaGroupMap = new HashMap(); //TreeValidatedGroupbyPloidy
            HashMap<String, String> taxaSubMap = new HashMap(); //TreeValidatedGroupbySubspecies
            HashMap<String, String> taxaGroupIDMap = new HashMap(); //IndexforMutationBurden
            RowTable t = new RowTable (taxaSummaryFileS);
            ArrayList<String> taxaList = new ArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                taxaGroupMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 14));
                taxaSubMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 15));
                taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 12));
                if (t.getCellAsDouble(i, 2) < depthCut) continue;
                taxaList.add(t.getCellAsString(i, 0));
            }
            String[] taxaWithHighDepth = taxaList.toArray(new String[taxaList.size()]); //taxa是具有高深度的taxa列表
            Arrays.sort(taxaWithHighDepth);

            BufferedWriter bw = AoFile.writeFile(addCountFileAddGroupS);
            t=new RowTable(addCountFileS);
            List<String> headerl = t.getHeader();
            for (int i = 0; i < headerl.size(); i++) {
                bw.write(headerl.get(i) + "\t");
            }
            bw.write("Group\tSubspecies\tGroupID\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) { //这里的t是 addCountFileS
                int index = Arrays.binarySearch(taxa, t.getCellAsString(i, 0));
                if (index < 0) continue;
//                if(taxaGroupMap.get(taxa[index]).equals("ExclusionHexaploid") || taxaGroupMap.get(taxa[index]).equals("ExclusionTetraploid")) continue; //这里不去除其他 四倍体 六倍体
                double genotypesite = Double.valueOf(t.getCellAsDouble(i, 2)); //本列指的是在该条染色体中含有有害突变的计数，如果物种不含本条染色体，那么数值为0，过滤掉。
                if(genotypesite == 0) continue; //如山羊草在A B亚基因组没有值，故这里删去
                double ratio = Double.valueOf(t.getCellAsDouble(i, 2))/Double.valueOf(t.getCellAsDouble(i, 3));
                for (int j = 0; j < t.getColumnNumber(); j++) { //按列书写
                    bw.write(t.getCellAsString(i,j) + "\t");
                }
                bw.write(taxaGroupMap.get(taxa[index])
                        +"\t"+taxaSubMap.get(taxa[index])
                        +"\t"+taxaGroupIDMap.get(taxa[index])
                        +"\t"+String.format("%.4f",ratio));
                bw.newLine();
            }
            bw.close();
            System.out.println("Finished step5: adding group info for taxa table.");
            System.out.println("Finished in making the del count table for each taxa in each single chromsome");
            System.out.println("-----------------------------------------------------------------------------");
            System.out.println("Now begin to merge final file by subgenome.");
            this.mergeFinalfilebySub(addCountFileAddGroupS); //##################### 重点！！！！！！合并文件！！！！！
            new File(addCountFileS).delete();

        }
        catch(Exception e){
            e.printStackTrace();
        }
    }




    public void DeltoSynonymousRatio(String infileS1, String ratioType, String infileS2){
//        String infileS1 = ""; //有害突变文件
//        String infileS2 = ""; //同义突变文件

        String outfileS = new File(infileS1).getAbsolutePath().replaceFirst(".txt","_delVSsynonymous.txt");
        AoFile.readheader(infileS1);

        TDoubleArrayList del = new TDoubleArrayList();
        TDoubleArrayList syn = new TDoubleArrayList();
        TDoubleArrayList delcount = new TDoubleArrayList();
        TDoubleArrayList syncount = new TDoubleArrayList();

        if(ratioType.equals("bySub")){
            del = AoFile.getTDoubleList(infileS1,7); //bySub的情况
            syn = AoFile.getTDoubleList(infileS2,7); //bySub的情况
            delcount = AoFile.getTDoubleList(infileS1,2); //bySub的情况
            syncount = AoFile.getTDoubleList(infileS2,2); //bySub的情况


        }else if (ratioType.equals("bysub_mergeByTaxa")){
            del = AoFile.getTDoubleList(infileS1,6); //byTaxa的情况
            syn = AoFile.getTDoubleList(infileS2,6); //byTaxa的情况
            delcount = AoFile.getTDoubleList(infileS1,1); //byTaxa的情况
            syncount = AoFile.getTDoubleList(infileS2,1); //byTaxa的情况
        }

        TDoubleArrayList ratioList = new TDoubleArrayList(del.size());
        for (int i = 0; i < del.size(); i++) {
            double r = (double)del.get(i)/syn.get(i); //String.format("%.3f",r)
            ratioList.add(r);
        }

        //用count除以count
        TDoubleArrayList ratioList2 = new TDoubleArrayList(delcount.size());
        for (int i = 0; i < delcount.size(); i++) {
            double r = (double)delcount.get(i)/syncount.get(i); //String.format("%.3f",r)
            ratioList2.add(r);
        }

        try {
            BufferedReader br = AoFile.readFile(infileS1);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = br.readLine();
            bw.write(temp + "\tRatio_delVSsyn\tRatio_delVSsyn_bycount");
            bw.newLine();
            int i = 0;
            while ((temp = br.readLine()) != null) {
                bw.write(temp + "\t" + String.format("%.3f",ratioList.get(i))+ "\t" + String.format("%.3f",ratioList2.get(i)));
                bw.newLine();
                i++;
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



    /**
     * 计算指定选择区域的mutation burden
     * 一共有 5 步
     */
    public void countDeleteriousVMapII_byChr(String infileS, String ifselected,String type, String addCountFileAddGroupS) {

        //model

        //######## 需要修改 ########//
//        String exonVCFDirS = ""; //外显子变异数据
//        String SNPAnnoFileS = ""; //注释信息库合并后的总文件
//        AoFile.readheader(SNPAnnoFileS);

        int cntNONSY = 0; //非同义突变的个数

        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################w
        String exonVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF"; //外显子变异数据
        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //注释信息库合并后的总文件
        AoFile.readheader(SNPAnnoFileS);
        //************* 无需修改的路径 ****************** //
        String addCountFileS = new File(addCountFileAddGroupS).getAbsolutePath().replaceFirst(".txt",".temp.txt"); //有害变异加性模型输出文件

        /**
         * ################################### step0: 建立受选择区域的集合，并在下文进行 posList 和 ancestral charList 构建时进行适当的过滤。
         */

        // I only need to know the chr pos information about the selected region
        TIntArrayList[] selectedPosList = new TIntArrayList[42];
        for (int i = 0; i < selectedPosList.length; i++) {
            selectedPosList[i] = new TIntArrayList();
        }
        RowTable<String> selectedT = new RowTable<>(infileS);
        for (int i = 0; i < selectedT.getRowNumber(); i++) {
            int chr = selectedT.getCellAsInteger(i,0);
            int pos = selectedT.getCellAsInteger(i,1);
            selectedPosList[chr-1].add(pos);
        }
        for (int i = 0; i < selectedPosList.length; i++) {
            selectedPosList[i].sort();
        }
        System.out.println("Finished step0: completing the initialization of selected list.");


        /**
         *  ################################### step1: 初始化染色体集合
         */
        int minDepth = 2;//inclusive
        int chrNum = 42;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrList.add(i + 1);
        }
        System.out.println("Finished step1: completing the initialization of chromosome.");

        /**
         *  ################################### step2: posList  charList 补充完整
         */

        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];

        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }

        String derivedAllele = null;

        try {
            BufferedReader br = AoFile.readFile(SNPAnnoFileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int chrID = Integer.parseInt(l.get(1));
                int index = Integer.parseInt(l.get(1)) - 1; //染色体号的索引 ################ 需要修改 需要修改 需要修改 ################
                int pos = Integer.parseInt(l.get(2)); //################ 需要修改 需要修改 需要修改 ################

                /**
                 * 进行受选择区域的判断,如果受选择，那么判断该位点是否受选择，不受选择的进行过滤。
                 * 如果不受选择，那么判断该位点是否受选择，受选择的进行过滤。
                 */

                int indexselect = selectedPosList[index].binarySearch(pos);
                if (ifselected.equals("1")){ //只保留受选择的位点，把不受选择的进行过滤
                    if (indexselect <0) continue;
                }
                if (ifselected.equals("0")){ //只保留不受选择的位点，把受选择的进行过滤
                    if (indexselect>= 0) continue;
                }

                String variantType = l.get(12); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(16); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(20); //################ 需要修改 需要修改 需要修改 ################

                /**
                 ******** 定义有害突变，不是有害突变，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */

                if (type.equals("001_synonymous")){
                    if (!variantType.equals("SYNONYMOUS"))continue;
                }
                if (type.equals("002_nonsynonymous")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue;
                }
                if (type.equals("003_nonsynGERPandDerivedSIFT")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                    if (sift.startsWith("N"))continue; //说明必须满足有sift值
                    double gerpd = Double.parseDouble(gerp);
                    double siftd = Double.parseDouble(sift);
                    if (gerpd < 1) continue; //说明必须满足gerp大于1
                    if (siftd > 0.05) continue; //说明必须满足sift小于等于0.05
                }

                //*********** 有害突变情况一: 非同义突变，且 GERP > 1**************//
//                if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
//                if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
//                double gerpd = Double.parseDouble(gerp);
//                if (gerpd < 1) continue;

                //################### 需要修改 //###################//###################//###################//###################
                ////不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
                String ancestralAllele = l.get(15);
                //################### 需要修改 //###################//###################//###################//###################

                String majorAllele = l.get(5);
                String minorAllele = l.get(6);
                if (ancestralAllele.equals(majorAllele)) {
                    derivedAllele = minorAllele;
                    posList[index].add(pos); //将包含有derived allele的位点添加到Poslist
                    charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
                }
                if (ancestralAllele.equals(minorAllele)) {
                    derivedAllele = majorAllele;
                    posList[index].add(pos);
                    charList[index].add(derivedAllele.charAt(0));
                }
                else if (!(ancestralAllele.equals(majorAllele) || ancestralAllele.equals(minorAllele))){
                }
                cnt++;
                if (cnt%100000==0) System.out.println("cnt is " + cnt +" going on at step2");

            }
            br.close();
            System.out.println(cntNONSY + " nonsynonymous SNP num");

            for (int i = 0; i < chrNum; i++) { //将每一个list转化为数组
                delePos[i] = posList[i].toArray();
                deleChar[i] = charList[i].toArray();
                Arrays.sort(delePos[i]);
            }
            System.out.println("Finished step2: completing the posList  charList.");

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        /**
         *  ################################### step3: taxa 集合 642个taxa
         */
//        String vmap2TaxaList = "/data4/home/aoyue/vmap2/analysis/000_taxaList/102_taxaInfoDB/taxa_InfoDB.txt";
        String vmap2TaxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";

        String[] taxa = AoFile.getStringArraybyList(vmap2TaxaList,0);

        double[][] addCount = new double[chrNum][taxa.length];
        int[][] recCount = new int[chrNum][taxa.length];
        int[][] siteWithMinDepthCount = new int[chrNum][taxa.length]; //每个taxa在每条染色体中的有害突变位点

        chrList.parallelStream().forEach(chr -> {
            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_exon_vmap2.1.vcf.gz";
            //开始读写VCF文件
            delVmapFileS = new File(exonVCFDirS, delVmapFileS).getAbsolutePath();
            BufferedReader br = AoFile.readFile(delVmapFileS);
            int chrIndex = chr - 1;
            try {
                String temp = null;
                int cnt = 0;
                HashMap<String, Integer> hmtaxainVCFindex = new HashMap<>();
                HashMap<String, Integer> hmtaxainTaxaindex = new HashMap<>();
                //这里涉及2个数组概念， taxa and taxainVCFfile,即总的taxa和在VCF文件中的taxa
                String[] taxainVCFfile = null;
                List<String> taxainVCFlist = new ArrayList();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##")) {
                        continue;
                    }
                    if (temp.startsWith("#CHROM")) {//说明进入taxa信息列
                        List<String> l = PStringUtils.fastSplit(temp, "\t");
                        taxainVCFfile = new String[l.size() - 9];
                        for (int i = 9; i < l.size(); i++) {
                            taxainVCFfile[i - 9] = l.get(i);
                            hmtaxainVCFindex.put(taxainVCFfile[i - 9], i); //第0个taxa的genotype在第9行，第一个taxa的genotype在第10行，依次类推；
                        }
                        // ************** 在总的taxa中，搜索VCF文件中的taxa的index
                        for (int i = 0; i < taxainVCFfile.length; i++) { //找到taxa
                            int index = Arrays.binarySearch(taxa, taxainVCFfile[i]);
                            if (index > -1) {
                                hmtaxainTaxaindex.put(taxainVCFfile[i], index);
                            }
                        }
                        taxainVCFlist = Arrays.asList(taxainVCFfile);
                    }

                    if (!temp.startsWith("#")) {
                        cnt++;
                        if (cnt % 10000 == 0) {
                            System.out.println(String.valueOf(cnt/1000) + " kb lines on chr " + String.valueOf(chr));
                        }
                        List<String> l = PStringUtils.fastSplit(temp.substring(0, 100), "\t");
                        int pos = Integer.valueOf(l.get(1)); //根据pos找到index
                        int index = Arrays.binarySearch(delePos[chrIndex], pos);
                        if (index < 0) {
                            continue; //前面建立的 delePos[][] 和deleChar[][] 都是为现在在vcf文件中找位置贡献的，不是有害突变的位点，都过滤。*************************************************************
                        } //重要：假如该位点有基因型，也是非同义突变位点，但是没有祖先位点状态，于是就过滤。
                        l = PStringUtils.fastSplit(temp, "\t");
                        int[] idx = new int[2];
                        if (l.get(3).charAt(0) == deleChar[chrIndex][index]) { //如果ref allele = deleChar allele
                            idx[0] = 0; //等于0则代表是 derived allele
                            idx[1] = 1;
                        } else {
                            idx[0] = 1;
                            idx[1] = 0;
                        }


                        //合计642个taxa，在A Bgenome中只有606（419+187）个，在Dgenome中只有455（419+36）个，我们要找到每个VCF文件中的genotype所对应的taxa的index
                        //code:本段代码是每行SNP位点，每个taxa的有害位点的统计
                        for (int i = 0; i < taxainVCFlist.size(); i++) {
                            int genotypeIndex = hmtaxainVCFindex.get(taxainVCFlist.get(i)); //获取该taxa的基因型所在的列的索引
                            int taxaIndex = hmtaxainTaxaindex.get(taxainVCFlist.get(i)); //获取该taxa所在总的642个数组中的索引，为后续写文件进行统计
                            String genoS = l.get(genotypeIndex); //指的是 GT:AD:GL 信息
                            if (genoS.startsWith(".")) {
                                continue; //如果以.开头，说明没有基因型信息，此位点没有测到。
                            }
                            List<String> ll = PStringUtils.fastSplit(genoS, ":"); //分开为GT   AD   GL三类

                            List<String> lll = PStringUtils.fastSplit(ll.get(1), ","); //lll指将AD提取出来，并以"，"号分割。如 0/0:1,0:0,3,28中 ，1，0分别代表ref和alt的测序深度
                            int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1)); //总得测序深度等于 ref + alt
                            if (depth < minDepth) {
                                continue; //最小测序深度是2，如果小于2，则弃用
                            }
                            lll = PStringUtils.fastSplit(ll.get(0), "/"); //这里lll指的是基因型GT，lll被重新赋值，之前代表的是AD
                            int v1 = Integer.valueOf(lll.get(0)); //v1等于 ref
                            int v2 = Integer.valueOf(lll.get(1)); // v2 等于 alt
                            int sum = 0;
                            //如果ref是derived allele，那么idx[0]=0;idx[1]=1. 当是0/0时，sum=2; 0/1时，sum=1; 1/1时， sum=0.
                            //如果alt是derived allele，那么idx[0]=1;idx[1]=0. 当是0/0时，sum=0; 0/1时，sum=1; 1/1时， sum=2.
                            //
                            //
                            if (v1 == idx[0]) { //
                                sum++;
                            }
                            if (v2 == idx[0]) {
                                sum++;
                            }
                            if (sum == 0) {
                            } else if (sum == 1) {
                                addCount[chrIndex][taxaIndex] += 0.5;
                            } else {
                                addCount[chrIndex][taxaIndex] += 1;
                                recCount[chrIndex][taxaIndex] += 1;
                            }
                            siteWithMinDepthCount[chrIndex][taxaIndex]++; //taxa有多少个有害突变位点
                        }
                    }
                }
                br.close();
                System.out.println("Finished step3: calculating each taxa for each chromosome.");

            } catch (Exception e) {
                e.printStackTrace();
            }

            /**
             * 进行每个样品每条染色体的计算
             */
        });

        /**
         *  ################################### step4: 输出每个样品每条染色体的结果
         */

        try {
            BufferedWriter bw = AoFile.writeFile(addCountFileS);
            bw.write("Taxa\tChr\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) { //第一层是染色体号
                String chr = String.valueOf(i + 1);
                for (int j = 0; j < addCount[0].length; j++) { //第二层是该号染色体的有害变异个数
                    bw.write(taxa[j] + "\t" + chr + "\t" + String.valueOf(addCount[i][j]) + "\t" + String.valueOf(siteWithMinDepthCount[i][j]));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();

            System.out.println("Finished step4: writing taxa table.");

        } catch (Exception e) {
            e.printStackTrace();
        }


        /**
         *  ################################### step5: 将输出的文件添加分组信息：每个taxa的倍性，亚群分布，mutation burden Index
         */
        try{
            String taxaSummaryFileS = vmap2TaxaList;
//            String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
            AoFile.readheader(taxaSummaryFileS);
            double depthCut = 1; //只保留深度大于3的taxa样本
            HashMap<String, String> taxaGroupMap = new HashMap(); //TreeValidatedGroupbyPloidy
            HashMap<String, String> taxaSubMap = new HashMap(); //TreeValidatedGroupbySubspecies
            HashMap<String, String> taxaGroupIDMap = new HashMap(); //IndexforMutationBurden
            RowTable t = new RowTable (taxaSummaryFileS);
            ArrayList<String> taxaList = new ArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                taxaGroupMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 14));
                taxaSubMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 15));
                taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 12));
                if (t.getCellAsDouble(i, 2) < depthCut) continue;
                taxaList.add(t.getCellAsString(i, 0));
            }
            String[] taxaWithHighDepth = taxaList.toArray(new String[taxaList.size()]); //taxa是具有高深度的taxa列表
            Arrays.sort(taxaWithHighDepth);

            BufferedWriter bw = AoFile.writeFile(addCountFileAddGroupS);
            t=new RowTable(addCountFileS);
            List<String> headerl = t.getHeader();
            for (int i = 0; i < headerl.size(); i++) {
                bw.write(headerl.get(i) + "\t");
            }
            bw.write("Group\tSubspecies\tGroupID\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) { //这里的t是 addCountFileS
                int index = Arrays.binarySearch(taxa, t.getCellAsString(i, 0));
                if (index < 0) continue;
//                if(taxaGroupMap.get(taxa[index]).equals("ExclusionHexaploid") || taxaGroupMap.get(taxa[index]).equals("ExclusionTetraploid")) continue; //这里不去除其他 四倍体 六倍体
                double genotypesite = Double.valueOf(t.getCellAsDouble(i, 2)); //本列指的是在该条染色体中含有有害突变的计数，如果物种不含本条染色体，那么数值为0，过滤掉。
                if(genotypesite == 0) continue; //如山羊草在A B亚基因组没有值，故这里删去
                double ratio = Double.valueOf(t.getCellAsDouble(i, 2))/Double.valueOf(t.getCellAsDouble(i, 3));
                for (int j = 0; j < t.getColumnNumber(); j++) { //按列书写
                    bw.write(t.getCellAsString(i,j) + "\t");
                }
                bw.write(taxaGroupMap.get(taxa[index])
                        +"\t"+taxaSubMap.get(taxa[index])
                        +"\t"+taxaGroupIDMap.get(taxa[index])
                        +"\t"+String.format("%.4f",ratio));
                bw.newLine();
            }
            bw.close();
            System.out.println("Finished step5: adding group info for taxa table.");
            System.out.println("Finished in making the del count table for each taxa in each single chromsome");
            System.out.println("-----------------------------------------------------------------------------");
            System.out.println("Now begin to merge final file by subgenome.");
            this.mergeFinalfilebySub(addCountFileAddGroupS); //##################### 重点！！！！！！合并文件！！！！！
            new File(addCountFileS).delete();

        }
        catch(Exception e){
            e.printStackTrace();
        }
    }


    //根据最终生成的文件，进行 A B D sub的合并,使每个taxa具有Asub Bsub Dsub的结果
    /// hahhahaahha
    public void mergeFinalfilebySub(String infileS){

        //change
//        String infileS = "";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/008_deleteriousRegion/002_countDel/001_additiveDeleterious_vmap2_bychr_selectedRegion.txt";
        String splitDirS = new File(infileS).getParent()+"/split";
        new File(splitDirS).mkdirs();

        String outfileS = infileS.replaceFirst(".txt", "_bysub.txt");

        /**
         * ############################# step1: split del table by taxa
         */
        //查看有多少个taxa
        RowTable<String> t = new RowTable<>(infileS);
        Set<String> s = new HashSet<>(t.getColumn(0));
        String[] taxa = s.toArray(new String[s.size()]);
        Arrays.sort(taxa);
        /**
         * 将文件按照taxa进行分类写出
         */
        try {
            BufferedReader br = AoFile.readFile(infileS);
            String header = br.readLine();
            BufferedWriter[] bw = new BufferedWriter[taxa.length];
            for (int i = 0; i < taxa.length; i++) {
                String tempoutS = new File(splitDirS,taxa[i]+".txt").getAbsolutePath();
                bw[i] = AoFile.writeFile(tempoutS);
                bw[i].write(header);
                bw[i].newLine();
            }
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String ta = l.get(0);
                int index = Arrays.binarySearch(taxa,ta);
                bw[index].write(temp);
                bw[index].newLine();
                cnt++;

            }
            br.close();
            for (int i = 0; i < taxa.length; i++) {
                bw[i].flush();
                bw[i].close();
            }
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


        //不变的
//        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        AoFile.readheader(taxaSummaryFileS);
        HashMap<String, String> taxaGroupMap = new HashMap();
        HashMap<String, String> taxaSubMap = new HashMap();
        HashMap<String, String> taxaGroupIDMap = new HashMap();
        t = new RowTable (taxaSummaryFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 14));
            taxaSubMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 15));
            taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 12));
        }

        HashMap<String,Integer> h = new HashMap<>();
        h.put("A",0);
        h.put("B",1);
        h.put("D",2);

        /**
         * ############################# step2: calculating each taxa's deleterious mutation by A sub Bsub and D sub
         */
        String outDirS = null;
        for (int i = 0; i < taxa.length; i++) {
            try {
                outDirS = new File(splitDirS).getParent() + "/split_merge";
                new File(outDirS).mkdirs();
                String inS = new File(splitDirS,taxa[i]+".txt").getAbsolutePath();
                String outS = new File(outDirS,taxa[i]+"_bysub.txt").getAbsolutePath();
                BufferedReader br = AoFile.readFile(inS);
                BufferedWriter bw = AoFile.writeFile(outS);
                bw.write(br.readLine().replaceFirst("Chr","Sub"));
                bw.newLine();
                TDoubleArrayList[] derivedDelList = new TDoubleArrayList[3];
                TDoubleArrayList[] genotypeList = new TDoubleArrayList[3];
                for (int j = 0; j < derivedDelList.length; j++) {
                    derivedDelList[j] = new TDoubleArrayList();
                    genotypeList[j] = new TDoubleArrayList();

                }
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    int chr = Integer.parseInt(l.get(1));
                    double d = Double.parseDouble(l.get(2));
                    double k = Double.parseDouble(l.get(3));
                    String sub = this.convertChrtoSub(chr);
                    derivedDelList[h.get(sub)].add(d);
                    genotypeList[h.get(sub)].add(k);
                }
                br.close();


                //求和
                DescriptiveStatistics[] d = new DescriptiveStatistics[3];
                DescriptiveStatistics[] dd = new DescriptiveStatistics[3];

                Double[] ratio = new Double[3];
                for (int j = 0; j < derivedDelList.length; j++) {
                    d[j] = new DescriptiveStatistics(derivedDelList[j].toArray());
                    dd[j] = new DescriptiveStatistics(genotypeList[j].toArray());
                    ratio[j] = d[j].getSum()/dd[j].getSum();

                }

                HashMap<Integer,String> hhh = new HashMap<>();
                hhh.put(0,"A");
                hhh.put(1,"B");
                hhh.put(2,"D");


                for (int j = 0; j < derivedDelList.length; j++) {
                    if(dd[j].getSum()==0)continue;
                    bw.write(taxa[i] + "\t" + hhh.get(j) + "\t" + String.format("%.1f",d[j].getSum()) + "\t" + String.format("%.0f",dd[j].getSum())
                            + "\t" + taxaGroupMap.get(taxa[i]) + "\t" + taxaSubMap.get(taxa[i]) + "\t" + taxaGroupIDMap.get(taxa[i]) + "\t" + String.format("%.4f",ratio[j]));
                    bw.newLine();
                }

                bw.flush();
                bw.close();

            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        /**
         * ############################# step3: merge each taxa's deleterious mutation by A sub Bsub and D sub
         */

        AoFile.mergeTxtbysuffix(outDirS,outfileS,".txt");
        new File(splitDirS).delete();
        new File(outDirS).delete();
        //######################## 依据taxa进行合并
        this.mergeByTaxa(outfileS);


    }

    /**
     * 根据染色体号，返回所在亚基因组
     * @param a
     * @return
     */
    public String convertChrtoSub(int a){
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
        String subgenome = hml.get(a);
        return subgenome;
    }


    /**
     * 将A B D的结果合并，使一个taxa拥有一个结果，不分亚基因组
     * 即一个taxa的所有有害等位位点数，基因型数，比率。
     */

    public void mergeByTaxa(String infileS){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/008_deleteriousRegion/002_countDel/001_additiveDeleterious_vmap2_bychr_selectedRegion_bysub.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/008_deleteriousRegion/002_countDel/001_WEvsDE_additiveDeleterious_vmap2_bychr_selectedRegion_bysub.txt";
        String outfileS = new File(infileS).getAbsolutePath().replaceFirst(".txt","_mergeByTaxa.txt");
        String[] taxa = AoFile.getStringArraybySet(infileS,0);

        TDoubleArrayList[] derivedDelList = new TDoubleArrayList[taxa.length];
        TDoubleArrayList[] genotypeList = new TDoubleArrayList[taxa.length];
        for (int j = 0; j < derivedDelList.length; j++) {
            derivedDelList[j] = new TDoubleArrayList();
            genotypeList[j] = new TDoubleArrayList();

        }


//        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";

        AoFile.readheader(taxaSummaryFileS);
        HashMap<String, String> taxaGroupMap = new HashMap();
        HashMap<String, String> taxaSubMap = new HashMap();
        HashMap<String, String> taxaGroupIDMap = new HashMap();
        RowTable<String> t = new RowTable (taxaSummaryFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 14));
            taxaSubMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 15));
            taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 12));
        }

        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = br.readLine(); //read header
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String taxon = l.get(0);
                int index = Arrays.binarySearch(taxa,taxon);
                if (index<0)continue;
                double d = Double.parseDouble(l.get(2));
                double k = Double.parseDouble(l.get(3));
                derivedDelList[index].add(d);
                genotypeList[index].add(k);
            }
            br.close();

            DescriptiveStatistics[] d = new DescriptiveStatistics[taxa.length];
            DescriptiveStatistics[] dd = new DescriptiveStatistics[taxa.length];

            Double[] ratio = new Double[taxa.length];
            for (int i = 0; i < derivedDelList.length; i++) {
                d[i] = new DescriptiveStatistics(derivedDelList[i].toArray());
                dd[i] = new DescriptiveStatistics(genotypeList[i].toArray());
                ratio[i] = d[i].getSum()/dd[i].getSum();
            }

            bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth\tGroup\tSubspecies\tGroupID\tRatio");
            bw.newLine();
            for (int i = 0; i < taxa.length; i++) {
                if(dd[i].getSum()==0)continue;
                bw.write(taxa[i]  + "\t" + String.format("%.1f",d[i].getSum()) + "\t" + String.format("%.0f",dd[i].getSum())
                        + "\t" + taxaGroupMap.get(taxa[i]) + "\t" + taxaSubMap.get(taxa[i]) + "\t" + taxaGroupIDMap.get(taxa[i]) + "\t" + String.format("%.4f",ratio[i]));
                bw.newLine();
            }

            bw.flush();
            bw.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

}
