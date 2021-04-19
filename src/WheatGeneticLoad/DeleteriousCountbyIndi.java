package WheatGeneticLoad;

import AoUtils.AoColor;
import AoUtils.AoFile;
import AoUtils.CountSites;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 * @author AoyueBi
 * @data 2020-06-11 22:38
 */
public class DeleteriousCountbyIndi {

    ////// 最新版结果
    public DeleteriousCountbyIndi(){

        this.countDeleteriousVMapII_byChr();
//        this.DeltoSynonymousRatio();
//        this.filterLandrace();

        //************* residual analysis ************//
//        this.addIBSdistancetoCS2017();
//        this.filterLandrace();
        //**Obstacle:由于我把642VCF全部合并求一个Dxy值，分析时出现偏差，故需要分开计算
//        this.mergeExonVCFbySub();
//        this.extractIBSdistance();
//        this.addIBSdistancebySub();
//        this.addGrouptoCorrelationFile();
//        this.addNewGrouptoLoadFile();

//        this.countDeleteriousVMapII_byChr_onlyHomo(); //统计分子是纯合，分母也是纯合基因型的3种分类； 统计分子是纯合，分母是纯合加杂合的3种分类； 统计分子是杂合，分母
//        this.DeltoSynonymousRatio();
//        this.filterLandrace();

        /**
         * ******** 计算个体load的时候，只看纯合子。 遇到杂合子就跳过
         */


    }

    //************************************************************************************ BEGIN ********************************************************************************************

    /**
     * 计算所有区域的mutation burden
     * 一共有 5 步
     *
     * 程序运行说明：每次运行程序，需更改的内容有：
     * 1.输出文件类型  是同义 非同义 有害
     * 2.判断条件
     * 3.数据库中对应的index要核查
     *
     * 4.程序运行完毕后，需将结果进行核验，确保万无一失，再进行画图。
     */
    public void countDeleteriousVMapII_byChr_onlyHomo() {
        //######## 需要修改 ########//
//        String exonVCFDirS = ""; //外显子变异数据
//        String SNPAnnoFileS = ""; //注释信息库合并后的总文件
//        AoFile.readheader(SNPAnnoFileS);

        int cntNONSY = 0; //非同义突变的个数
//**************** VMap2.0-2020 加上 derived SIFT 的数据库 ************************* //
        String exonVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF"; //外显子变异数据
        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //注释信息库合并后的总文件
        AoFile.readheader(SNPAnnoFileS);
        int test =3;

//        new AoFile().readheader(SNPAnnoFileS);
//        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt");
//        System.out.println("************************************************************");
//        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt");
//        int a =3;

        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################
        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################
        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################
        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################

        //###### synonymous
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/005_homo_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/005_homo_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";

        //###### non-synonymous
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/006_homo_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/006_homo_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";

        //###### deleterious by GERP and derived SIFT (分子是有害突变的纯合子计数， 分母是该位点在全基因组是有害突变，并且有基因型，并且是纯合子)
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/007_homo_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr.txt"; //指的是非同义突变，并且GERP大于1,SIFT<0.05的条件
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/007_homo_recessive_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr.txt"; //指的是非同义突变，并且GERP大于1,SIFT<0.05的条件



        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################ 分母是该类型的突变有基因型，分子是纯合

        //###### synonymous
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/001_homoTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/001_homoTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";

        //###### non-synonymous
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/001_homoTotal_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/001_homoTotal_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";

        //###### deleterious by GERP and derived SIFT (分子是有害突变的纯合子计数， 分母是该位点在全基因组是有害突变，并且有基因型，并且是纯合子)
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/003_homoTotal_nonsynGERPandDerivedSIFT_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/003_homoTotal_nonsynGERPandDerivedSIFT_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";


        ////########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################ 分母是该类型的突变有基因型，分子是纯合

        ////###### synonymous
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/005_heterTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";

//        //###### non-synonymous
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/006_heterTotal_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";
//
//        //###### deleterious by GERP and derived SIFT (分子是有害突变的纯合子计数， 分母是该位点在全基因组是有害突变，并且有基因型，并且是纯合子)
        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/007_heterTotal_nonsynGERPandDerivedSIFT_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";





//************* 无需修改的路径 ****************** //
        String addCountFileS = new File(addCountFileAddGroupS).getAbsolutePath().replaceFirst(".txt",".temp.txt"); //有害变异加性模型输出文件
//

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
//                String sift = l.get(13); //################ 需要修改 需要修改 需要修改 ################
//                String gerp = l.get(18); //################ 需要修改 需要修改 需要修改 ################

                String sift = l.get(16); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(20); //################ 需要修改 需要修改 需要修改 ################

                /**
                 ******** 定义有害突变，不是有害突变，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */

                //*********** 有害突变情况一: 非同义突变，且 GERP > 1**************//
//                if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
//                if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
//                double gerpd = Double.parseDouble(gerp);
//                if (gerpd < 1) continue;

                //*********** 有害突变情况二: 非同义突变，且 GERP >= 1，  SIFT <= 0.05 **************//
                if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                if (sift.startsWith("N"))continue; //说明必须满足有sift值
                double gerpd = Double.parseDouble(gerp);
                double siftd = Double.parseDouble(sift);
                if (gerpd < 1) continue; //说明必须满足gerp大于1
                if (siftd > 0.05) continue; //说明必须满足sift小于等于0.05

                /**
                 ********* 定义 synonymous, 不是 synonymous,就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */
//                if (!variantType.equals("SYNONYMOUS"))continue;

                /**
                 ********** 定义 non-synonymous ，不是 non-synonymous，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */
//                if (!variantType.equals("NONSYNONYMOUS"))continue;


                //################### 需要修改 //###################//###################//###################//###################
                ////不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
//                String ancestralAllele = l.get(14);
                String ancestralAllele = l.get(15);
                //################### 需要修改 //###################//###################//###################//###################

                String majorAllele = l.get(5);
                String minorAllele = l.get(6);
                if (ancestralAllele.equals(majorAllele)) {
                    derivedAllele = minorAllele;
                    posList[index].add(Integer.parseInt(l.get(2))); //将包含有derived allele的位点添加到Poslist
                    charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
                }
                if (ancestralAllele.equals(minorAllele)) {
                    derivedAllele = majorAllele;
                    posList[index].add(Integer.parseInt(l.get(2)));
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
//            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_exon_vmap2.1_filterbyHeter0.05.vcf.gz";
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

                            /**
                             * 添加一步杂合子过滤
                             *
                             */
//                            if (ll.get(0).equals("0/1")) continue;

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
//                                siteWithMinDepthCount[chrIndex][taxaIndex]++; //taxa有多少个有害突变位点

                            } else if (sum == 1) {
                                addCount[chrIndex][taxaIndex] += 0.5;
                            } else if (sum ==2){
//                                addCount[chrIndex][taxaIndex] += 1;
//                                siteWithMinDepthCount[chrIndex][taxaIndex]++; //taxa有多少个有害突变位点
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
                    bw.write(taxa[j] + "\t" + chr + "\t" + addCount[i][j] + "\t" + siteWithMinDepthCount[i][j]);
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

    //************************************************************************************ END ********************************************************************************************



    /**
     * 在load结果中添加landrace的分组7，和cultivar的分组7
     */
    public void addNewGrouptoLoadFile(){
//        String hmfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
//        AoFile.readheader(hmfileS);
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/005/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous.txt";
//        int[] columnIndexes = {20,21};
//        HashMap<String,String>[] hm = new AoFile().getHashMapsStringKey(hmfileS,0,columnIndexes);
//        AoFile.addColumsbyString(outfileS,0,hm,"Landrace7Cultivar1_byContinent\tLandrace7Cultivar7_byContinent");

        String hmfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        AoFile.readheader(hmfileS);
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/nonsyn/002/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_delVSsynonymous.txt";
        int[] columnIndexes = {20,21};
        HashMap<String,String>[] hm = new AoFile().getHashMapsStringKey(hmfileS,0,columnIndexes);
        AoFile.addColumsbyString(outfileS,0,hm,"Landrace7Cultivar1_byContinent\tLandrace7Cultivar7_byContinent");

    }

    /**
     * 根据R语言得到的20个小群中 load 和 ibs 的相关性 value1 以及ibs值 value2，再添加具体倍性信息，可以在画图时进行分组
     */
    public void addGrouptoCorrelationFile(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/006/001_loadCorreandIBS.txt";
        String taxainfoS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt";
        AoFile.readheader(taxainfoS);
        HashMap<String,String> hm = AoFile.getHashMapStringKey(taxainfoS,12,3);
        AoFile.addColumbyString(infileS,1,hm,"GenomeType");

    }

    //，再添加whole genome VCF文件的IBS结果，然后再进行分开残差校正
    public void addIBSdistancebySub(){
//        this.step1();
//        this.step2();
        this.step3();

    }

    private void step3(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/014_VMap2.1DelCount_derivedSIFT_filterLR_CL/splitbySub/004";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/014_VMap2.1DelCount_derivedSIFT_filterLR_CL/splitbySub/005/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous.txt";
//        AoFile.mergeTxt(infileDirS,outfileS);

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/004";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/005/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous.txt";
//        AoFile.mergeTxt(infileDirS,outfileS);

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/nonsyn/001";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/nonsyn/002/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_delVSsynonymous.txt";
        AoFile.mergeTxt(infileDirS,outfileS);

    }

    /**
     * add IBS
     */
    private void step2(){
//        String infileS = "";
//        String ibsFileS = "";

        //过滤 LR 只剩欧洲的
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/014_VMap2.1DelCount_derivedSIFT_filterLR_CL/splitbySub/004/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous_A.txt";
//        String ibsFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/001_IBSdistance2CS2017/IBSdistance_byAsub.txt";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/014_VMap2.1DelCount_derivedSIFT_filterLR_CL/splitbySub/004/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous_B.txt";
//        String ibsFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/001_IBSdistance2CS2017/IBSdistance_byBsub.txt";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/014_VMap2.1DelCount_derivedSIFT_filterLR_CL/splitbySub/004/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous_D.txt";
//        String ibsFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/001_IBSdistance2CS2017/IBSdistance_byDsub.txt";

        //所有taxa不进行欧洲landrace过滤
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/004/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous_A.txt";
//        String ibsFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/001_IBSdistance2CS2017/IBSdistance_byAsub.txt";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/004/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous_B.txt";
//        String ibsFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/001_IBSdistance2CS2017/IBSdistance_byBsub.txt";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/004/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous_D.txt";
//        String ibsFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/001_IBSdistance2CS2017/IBSdistance_byDsub.txt";


        //所有taxa不进行欧洲landrace过滤 nonsynonymous

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/nonsyn/001/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_delVSsynonymous_A.txt";
//        String ibsFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/001_IBSdistance2CS2017/IBSdistance_byAsub.txt";
//
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/nonsyn/001/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_delVSsynonymous_B.txt";
//        String ibsFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/001_IBSdistance2CS2017/IBSdistance_byBsub.txt";

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/nonsyn/001/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_delVSsynonymous_D.txt";
        String ibsFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/001_IBSdistance2CS2017/IBSdistance_byDsub.txt";

        HashMap<String,String> hm = AoFile.getHashMapStringKey(ibsFileS,0,1);
        AoFile.addColumbyString(infileS,0,hm,"IBSdistancetoCS2017");
    }

    /**
     * 拆分文件
     */
    public void step1(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous.txt";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/004";
//        AoFile.splitFilebyGroup(infileS,1,outfileDirS);

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/014_VMap2.1DelCount_derivedSIFT_filterLR_CL/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_delVSsynonymous.txt";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/014_VMap2.1DelCount_derivedSIFT_filterLR_CL/splitbySub/004";
//        AoFile.splitFilebyGroup(infileS,1,outfileDirS);

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_delVSsynonymous.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/splitbySub/nonsyn/001";
        AoFile.splitFilebyGroup(infileS,1,outfileDirS);


    }


    public void extractIBSdistance(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/004_ExonVCF_IBSdistance/source";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/004_ExonVCF_IBSdistance/001_IBSdistance2CS2017";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/source";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/005_subsetVCF_IBSdistance/001_IBSdistance2CS2017";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f ->{
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
            CountSites.extractIBSdistanceFromMatrix(infileS,outfileS,"CS-2017");
        });
    }

    /**
     * 计算exonVCF中不同亚基因组下的遗传距离，看看差别有多少
     */
    public void mergeExonVCFbySub(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/008_exonVCF/001_exonVCF/000_exonVCF_bySub";
        CountSites.mergeVCFbysubgenome(infileDirS,outfileDirS);
    }



    public void getLargeDelAnnotation(){
        int cntNONSY = 0;
//        String exonVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/002_exonSNPVCF"; //有害变异的VCF文件路径
//        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/004_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //有害变异信息库

        String exonVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/010_exonSNPVCF_filterHeter0.05"; //有害变异的VCF文件路径
        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/012_exonSNPAnnotation_merge_filterHeter0.05/001_exonSNP_anno_filterHeter0.05.txt.gz"; //有害变异信息库

        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/008_largeDelAnnotation/exonSNP_anno_filterHeter0.05_largeDelgerp3.txt";
        new AoFile().readheader(SNPAnnoFileS);

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
            BufferedWriter bw = AoFile.writeFile(outfileS);

            String temp = null;
            String header = br.readLine();
            bw.write(header+"\tSub");bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int index = Integer.parseInt(l.get(1)) - 1; //染色体号的索引
                int chr = Integer.parseInt(l.get(1));
                int pos = Integer.parseInt(l.get(2));
                String sub = RefV1Utils.getChromosome(chr,pos).substring(1);
                String trans = l.get(12);
                String variantType = l.get(12);
                String sift = l.get(13);
                String gerp = l.get(20);
                /**
                 * 定义有害突变，不是有害突变，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */

                if (!variantType.equals("NONSYNONYMOUS") || sift.equals("NA") || gerp.equals("NA"))continue;

//                if (!variantType.equals("NONSYNONYMOUS") || sift.equals("NA"))continue;
//                if (!variantType.equals("NONSYNONYMOUS") || gerp.equals("NA"))continue;
                double siftd = Double.parseDouble(sift);
                double gerpd = Double.parseDouble(gerp);

                if (siftd >= 0.05 || gerpd < 3)continue;



                /**
                 * 定义有害突变，不是有害突变，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */
//                if (!variantType.equals("SYNONYMOUS"))continue;

                /**
                 * 定义有害突变，不是有害突变，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */
//                if (!variantType.equals("NONSYNONYMOUS"))continue;

                /**
                 * 只是用来计数用 ################ 需要修改 需要修改 需要修改 ################
                 */

//                if (variantType.equals("NONSYNONYMOUS") ){
//                    cntNONSY++;
//                }

//                if (variantType.equals("NONSYNONYMOUS") && (siftd < 0.05) && (gerpd > 1) ){
//                    cntNONSY++;
//                }

//                if (variantType.equals("NONSYNONYMOUS") && (siftd < 0.05)){
//                    cntNONSY++;
//                }

//                if (variantType.equals("NONSYNONYMOUS") && (gerpd > 1)){
//                    cntNONSY++;
//                }

                //################### 需要修改 //###################//###################//###################//###################
//                String ancestralAllele = l.get(22); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
//                String ancestralAllele = l.get(15); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
                String ancestralAllele = l.get(31);
                //################### 需要修改 //###################//###################//###################//###################

                String majorAllele = l.get(5);
                String minorAllele = l.get(6);
                if (ancestralAllele.equals(majorAllele)) {
                    derivedAllele = minorAllele;
                    posList[index].add(Integer.parseInt(l.get(2))); //将包含有derived allele的位点添加到Poslist
                    charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
                    bw.write(temp + "\t" + sub);
                    bw.newLine();
                }
                if (ancestralAllele.equals(minorAllele)) {
                    derivedAllele = majorAllele;
                    posList[index].add(Integer.parseInt(l.get(2)));
                    charList[index].add(derivedAllele.charAt(0));
                    bw.write(temp + "\t" + sub);
                    bw.newLine();
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

    }

    public void mergeExonSNPAnnotation(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/011_exonSNPAnnotation_filterHeter0.05";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/012_exonSNPAnnotation_merge_filterHeter0.05/001_exonSNP_anno_filterHeter0.05.txt.gz";
        AoFile.mergeTxt(infileDirS,outfileS);
    }

    public void getFiltedExonSNPAnnotation() {
        String exonVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/010_exonSNPVCF_filterHeter0.05";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/011_exonSNPAnnotation_filterHeter0.05";
        String exonAnnotationDBS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony";
        List<File> fsList = AoFile.getFileListInDir(exonAnnotationDBS);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_filterHeter0.05.txt.gz").getAbsolutePath();
                String exonVCFS = new File(exonVCFDirS,f.getName().split("_")[0]+"_exon_vmap2.1_filterbyHeter0.05.vcf.gz").getAbsolutePath();
                //获取vcf文件中的pos信息库
                TIntArrayList posl = AoFile.getNumListfromVCF(exonVCFS);

                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                String header = br.readLine();
                bw.write(header);
                bw.newLine();
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(2));
                    int index = posl.binarySearch(pos);
                    if (index > -1){
                        bw.write(temp);
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
    }


    public void DeltoSynonymousRatio(){
//        String infileS1 = ""; //有害突变文件
//        String infileS2 = ""; //同义突变文件

        //########### 大麦和黑麦简约法 VMap2.0-2020 结果
        //by Taxa
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件

//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/003_additiveDeleterious_nonsynGERP_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件

//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/004_additiveDeleterious_nonsynGERPandSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件


        //by Sub
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件

//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/003_additiveDeleterious_nonsynGERP_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件

//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/004_additiveDeleterious_nonsynGERPandSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件

        //########### 大麦和黑麦简约法 VMap2.0-2020 2020-07-21 增加 derived sift 结果
        //****** by taxa ******
        //nonsyn/syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件
        //nonsyn del/ syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件

        //****** by sub ******
        //nonsyn/syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件
        //nonsyn del/ syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件


        ////////////////////////// ********************** 纯杂合分子 ********************* /////////////////////////////////
        /**
         * 把所有杂合子都去掉，只分析纯合状态下，所有是纯合子，分母也是纯合子的load的高低情况
         */
        //****** by taxa ******
        //nonsyn/syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/006_homo_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件
//                String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/005_homo_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件

        //del/syn
                String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/007_homo_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/005_homo_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件

        //****** by sub ******
        //nonsyn/syn
//                String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/006_homo_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/005_homo_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件

        //del/syn
//                String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/007_homo_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/005_homo_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件


        /**
         * 分析纯纯纯纯纯合状态下，load的高低情况,分母是所有该类突变基因型在这个个体中的计数
         */
        //########### 大麦和黑麦简约法 VMap2.0-2020 2020-07-21 增加 derived sift 结果
        //****** by taxa ******
        //nonsyn/syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/002_homoTotal_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/001_homoTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件
        // del/syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/003_homoTotal_nonsynGERPandDerivedSIFT_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/001_homoTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件


        //****** by sub ******
        //nonsyn/syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/002_homoTotal_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/001_homoTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件
        // del/syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/003_homoTotal_nonsynGERPandDerivedSIFT_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/001_homoTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件


        /**
         * 分析杂杂杂杂杂杂杂杂合状态下，load的高低情况,分母是所有该类突变基因型在这个个体中的计数
         */

        //****** by taxa ******
        //nonsyn/syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/006_heterTotal_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/005_heterTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件
        //del/syn
//                String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/007_heterTotal_nonsynGERPandDerivedSIFT_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/005_heterTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub_mergeByTaxa.txt"; //同义突变文件


        //****** by sub ******
        //nonsyn/syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/006_heterTotal_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/005_heterTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件
        //del/syn
//        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/007_heterTotal_nonsynGERPandDerivedSIFT_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //有害突变文件
//        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter/005_heterTotal_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr_bysub.txt"; //同义突变文件





        String outfileS = new File(infileS1).getAbsolutePath().replaceFirst(".txt","_delVSsynonymous.txt");
        AoFile.readheader(infileS1);

//        TDoubleArrayList del = AoFile.getTDoubleList(infileS1,7); //bySub的情况
//        TDoubleArrayList syn = AoFile.getTDoubleList(infileS2,7); //bySub的情况
//        TDoubleArrayList delcount = AoFile.getTDoubleList(infileS1,2); //bySub的情况
//        TDoubleArrayList syncount = AoFile.getTDoubleList(infileS2,2); //bySub的情况

        TDoubleArrayList del = AoFile.getTDoubleList(infileS1,6); //byTaxa的情况
        TDoubleArrayList syn = AoFile.getTDoubleList(infileS2,6); //byTaxa的情况
        TDoubleArrayList delcount = AoFile.getTDoubleList(infileS1,1); //byTaxa的情况
        TDoubleArrayList syncount = AoFile.getTDoubleList(infileS2,1); //byTaxa的情况

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
     * Goal: 从 IBS distance matrix 中获取所有taxa 到 CS-2017 的距离值,然后添加到load表单中，进行残差分析
     * ##IBS_Distance_Matrix.AverageTotalSites=191737.66392009504
     * ##IBS_Distance_Matrix.NumAlleles=3
     * ##IBS_Distance_Matrix.TrueIBS=false
     * ##Matrix_Type=IBS_Distance_Matrix
     * 644
     */
    public void addIBSdistancetoCS2017(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/001_IBSdistance/001_IBSdistance2CS2017.txt";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT";
//        List<File> fsList = IOUtils.getFileListInDirEndsWith(infileDirS,"mergeByTaxa_delVSsynonymous.txt");
        List<File> fsList = IOUtils.getFileListInDirEndsWith(infileDirS,"bysub_delVSsynonymous.txt");

        Collections.sort(fsList);
        HashMap<String,String> hm = AoFile.getHashMapStringKey(infileS,0,1);
        for (int i = 0; i < fsList.size(); i++) {
            String outfileS = fsList.get(i).getAbsolutePath();
            AoFile.addColumbyString(outfileS,0,hm,"IBSdistancetoCS2017");
        }
    }

    /**
     * 将按照sub和taxa产生的有害突变的结果，进行过滤，使六倍体栽培品种都为欧洲的栽培品种，六倍体的农家种都为欧洲的农家种
     * 然后进行个体Load的判断
     */
    public void filterLandrace(){
//        String infileDirS = "";
//        String outfileDirS = "";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/004_VMap2.1DelCount_filterLR_CL";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/014_VMap2.1DelCount_derivedSIFT_filterLR_CL";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/002_addIBS";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/015_IBSdistanceAdjust/003_filterLR";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/homoAndHeter";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/014_VMap2.1DelCount_derivedSIFT_filterLR_CL/homoAndHeter_filterLR";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/014_VMap2.1DelCount_derivedSIFT_filterLR_CL";


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
     * 下一阶段：通过上一步骤得到的结果，进行指定群体的结果提取，进行Mutation burden 测试
     * need 4 file path
     */
    public void getPopmutationBurden(){
        String objectPop = "CL"; //pop1
        String refPop = "EU"; //pop2

        //输入输出文件  可变
//        String infileS = "";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/003_VMap2.1DelCount/002_additiveDeleterious_ANCbarleyVSsecale_vmap2_bychr_bysub_mergeByTaxa.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/003_VMap2.1DelCount/002_additiveDeleterious_ANCbarleyVSsecale_vmap2_bychr_bysub.txt";

        String outfileS = new File(infileS).getAbsolutePath().replaceFirst(".txt","_"+ objectPop + "_vs_" + refPop+".txt");
        //分组文件 可变
        String pop1fileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/008_treeValidatedGroup_bySubspecies/Cultivar.txt";
        String pop2fileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/009_treeValidatedFroup_byRegion/002_Landrace_European/Landrace_Europe.txt";

//        String pop1fileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/008_treeValidatedGroup_bySubspecies/Domesticated_emmer.txt";
//        String pop2fileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/008_treeValidatedGroup_bySubspecies/Wild_emmer.txt";
        String[] pop1 = AoFile.getStringArraybyList_withoutHeader(pop1fileS,0);
        String[] pop2 = AoFile.getStringArraybyList_withoutHeader(pop2fileS,0);
        System.out.println(pop1.length);
        System.out.println(pop2.length);
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = br.readLine();
            String pop = null;
            bw.write(temp+"\tGroup_bySubsubspecies");
            bw.newLine();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String taxon = l.get(0);
                int index1 = Arrays.binarySearch(pop1,taxon);
                int index2 = Arrays.binarySearch(pop2,taxon);
                if (index1 <0 && index2<0)continue;
                if (index1 > -1){
//                    pop = "Cultivar";
                    pop = objectPop;
                }
                if (index2 > -1){
//                    pop = "Landrace_Europe";
                    pop = refPop;
                }

                bw.write(temp+"\t"+pop);
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
     * 计算所有区域的mutation burden
     * 一共有 5 步
     */
    public void countDeleteriousVMapII_byChr() {
        //######## 需要修改 ########//
//        String exonVCFDirS = ""; //外显子变异数据
//        String SNPAnnoFileS = ""; //注释信息库合并后的总文件
//        AoFile.readheader(SNPAnnoFileS);

        int cntNONSY = 0; //非同义突变的个数
//**************** VMap2.0-2020 ************************* //
//        String exonVCFDirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/002_exonSNPVCF";
//        String SNPAnnoFileS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/004_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //有害变异信息库
        //本地信息库
//        String exonVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
//        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/015_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //有害变异信息库

//**************** VMap2.0-2020 加上 derived SIFT 的数据库 ************************* //
        String exonVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF"; //外显子变异数据
        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/019_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //注释信息库合并后的总文件
        AoFile.readheader(SNPAnnoFileS);
        int test =3;

//        new AoFile().readheader(SNPAnnoFileS);
//        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/011_taxaInfoDB/taxa_InfoDB.txt");
//        System.out.println("************************************************************");
//        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt");
//        int a =3;
        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** new data test ################
        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** new data test ################
        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** new data test ################
        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** new data test ################

        //###### synonymous
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/001_recessiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";

        //###### non-synonymous
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/002_recessiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";

        //###### deleterious by GERP
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/003_additiveDeleterious_nonsynGERP_ANCbarleyVSsecalePasimony_vmap2_bychr.txt"; //指的是非同义突变，并且GERP大于1的条件
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/003_recessiveDeleterious_nonsynGERP_ANCbarleyVSsecalePasimony_vmap2_bychr.txt"; //指的是非同义突变，并且GERP大于1的条件

        //###### deleterious by GERP and SIFT
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/004_additiveDeleterious_nonsynGERPandSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr.txt"; //指的是非同义突变，并且GERP大于1,SIFT<0.05的条件
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/003_VMap2.1DelCount/004_recessiveDeleterious_nonsynGERPandSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr.txt"; //指的是非同义突变，并且GERP大于1,SIFT<0.05的条件

        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################
        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################
        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################
        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################

        //###### synonymous
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/001_additiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/001_recessiveDeleterious_synonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";

        //###### non-synonymous
//        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/002_additiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";
//        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/002_recessiveDeleterious_nonsynonymous_ANCbarleyVSsecaleParsimony_vmap2_bychr.txt";

        //###### deleterious by GERP and derived SIFT
        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/004_additiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr.txt"; //指的是非同义突变，并且GERP大于1,SIFT<0.05的条件
        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/033_annoDB/013_VMap2.1DelCount_derivedSIFT/004_recessiveDeleterious_nonsynGERPandDerivedSIFT_ANCbarleyVSsecalePasimony_vmap2_bychr.txt"; //指的是非同义突变，并且GERP大于1,SIFT<0.05的条件


//************* 无需修改的路径 ****************** //
        String addCountFileS = new File(addCountFileAddGroupS).getAbsolutePath().replaceFirst(".txt",".temp.txt"); //有害变异加性模型输出文件
        String recCountFileS = new File(recCountFileAddGroupS).getAbsolutePath().replaceFirst(".txt",".temp.txt"); //有害变异隐形模型输出文件


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
//                String sift = l.get(13); //################ 需要修改 需要修改 需要修改 ################
//                String gerp = l.get(18); //################ 需要修改 需要修改 需要修改 ################

                String sift = l.get(16); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(20); //################ 需要修改 需要修改 需要修改 ################

                /**
                 ******** 定义有害突变，不是有害突变，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */

                //*********** 有害突变情况一: 非同义突变，且 GERP > 1**************//
//                if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
//                if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
//                double gerpd = Double.parseDouble(gerp);
//                if (gerpd < 1) continue;

                //*********** 有害突变情况二: 非同义突变，且 GERP >= 1，  SIFT <= 0.05 **************//
                if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                if (sift.startsWith("N"))continue; //说明必须满足有sift值
                double gerpd = Double.parseDouble(gerp);
                double siftd = Double.parseDouble(sift);
                if (gerpd < 1) continue; //说明必须满足gerp大于1
                if (siftd > 0.05) continue; //说明必须满足sift小于等于0.05

                /**
                 ********* 定义 synonymous, 不是 synonymous,就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */
//                if (!variantType.equals("SYNONYMOUS"))continue;

                /**
                 ********** 定义 non-synonymous ，不是 non-synonymous，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */
//                if (!variantType.equals("NONSYNONYMOUS"))continue;


                //################### 需要修改 //###################//###################//###################//###################
                ////不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
//                String ancestralAllele = l.get(14);
                String ancestralAllele = l.get(15);
                //################### 需要修改 //###################//###################//###################//###################

                String majorAllele = l.get(5);
                String minorAllele = l.get(6);
                if (ancestralAllele.equals(majorAllele)) {
                    derivedAllele = minorAllele;
                    posList[index].add(Integer.parseInt(l.get(2))); //将包含有derived allele的位点添加到Poslist
                    charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
                }
                if (ancestralAllele.equals(minorAllele)) {
                    derivedAllele = majorAllele;
                    posList[index].add(Integer.parseInt(l.get(2)));
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
//            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_exon_vmap2.1_filterbyHeter0.05.vcf.gz";
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

            bw = IOUtils.getTextWriter(recCountFileS);
            bw.write("Taxa\tChr\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) {
                String chr = String.valueOf(i + 1);
                for (int j = 0; j < addCount[0].length; j++) {
                    bw.write(taxa[j] + "\t" + chr + "\t" + String.valueOf(recCount[i][j]) + "\t" + String.valueOf(siteWithMinDepthCount[i][j]));
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
            new File(recCountFileS).delete();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
}
